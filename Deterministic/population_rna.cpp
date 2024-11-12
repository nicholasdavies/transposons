#include "population_rna.h"

#include <iomanip>
#include <vector>
#include <algorithm>
#include <numeric>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/weighted_mean.hpp>
#include <boost/accumulators/statistics/weighted_variance.hpp>

// Construct population
Population::Population()
 : A(P.n_max),
   B(P.n_max),
   f(ZeroTriangularish(P.n_max, P.n_max)),
   w(ZeroTriangularish(P.n_max, P.n_max)),
   N_A(0), N_B(0), V_A(0), V_B(0),
   max_freq(0), still(0), slow(0), f_other(0),
   saved_sWU(0), saved_b(0),

   DuplicationP_A(A, Vector(A, 0)),
   DuplicationP_B(B, Vector(B, 0)),
   ExcisionP_A(A, Vector(A, 0)),
   ExcisionP_B(B, Vector(B, 0)),
   TransitionP_A(A, Vector(A, 0)),
   TransitionP_B(B, Vector(B, 0)),
   GameteP_A(A, Vector(A, 0)),
   GameteP_B(B, Vector(B, 0)),

   W_lookup(),

   f_temp(ZeroTriangularish(P.n_max, P.n_max)),
   gametes(ZeroTriangularish(P.n_max, P.n_max)),
   offspring(ZeroTriangularish(P.n_max, P.n_max)),
   gamete_is()
{
}

// Population destructor
Population::~Population()
{
}

// Return the number of columns in row r of the frequency matrix.
int Population::Minor(int r) const
{
    if (B == 1) return 1;
    else        return min(B, A - r);
}

// Initialize population
void Population::Init()
{
    MakeLookupTables();

    // Genotype frequencies -- initially Poisson distributed, all wild-type
    for (int a = 0; a < A; ++a)
    {
        for (int b = 0; b < Minor(a); ++b)
        {
            if (b == 0)
                f[a][b] = PoissonPMF(P.n_0, a);
            else
                f[a][b] = 0;
        }
    }

    Normalize(f);
    CalculateStatistics();
}

// Do lookup tables for fitness
void Population::ConstructFitness()
{
    W_lookup.clear();
    W_lookup.reserve(A * B);

    for (int a = 0; a < A; ++a)
    {
        for (int b = 0; b < Minor(a); ++b)
        {
            w[a][b] = Fitness(a, P.u, b, P.u_B);
            W_lookup.push_back({a, b, w[a][b]});
        }
    }

    sort(W_lookup.begin(), W_lookup.end(), [](Lookup& a, Lookup& b) { return a.w == b.w ? (a.a == b.a ? a.b < b.b : a.a < b.a) : a.w < b.w; });

    W_list.resize(1, W_lookup.front().w);
    for (unsigned int i = 0; i < W_lookup.size(); ++i)
        if (W_lookup[i].w != W_list.back())
            W_list.push_back(W_lookup[i].w);
    W_list.push_back(1);
}

// Populate lookup tables
void Population::MakeLookupTables()
{
    ConstructFitness();

    // Calculate probabilities of duplications, excisions, and gamete formation.
    // The model assumes the following sequence for transposon events:
    // 1) Each transposon produces transposases through a Poisson process at average rate u.
    // 2) Each transposon has a fixed probability v of being excised/lost.
    // 3) The transposases do their work, creating a new transposon each regardless of whether their producing transposon has been lost.
    // In this way, the number of duplications is D ~ Poisson(n * u) and the number of excisions is E ~ Binomial(n, v).
    auto basic = [&](Matrix& dup, Matrix& exc, Matrix& gam, int G, double u, double v)
    {
        for (int g = 0; g < G; ++g)
        {
            for (int e = 0; e < G; ++e)
            {
                dup[g][e] = PoissonPMF(g * u, e);
                exc[g][e] = BinomialPMF(g, v, e);
                gam[g][e] = BinomialPMF(g, 0.5, e);
            }
        }
    };
    basic(DuplicationP_A, ExcisionP_A, GameteP_A, A, P.u, P.v);
    basic(DuplicationP_B, ExcisionP_B, GameteP_B, B, P.u_B, P.v);

    // Calculate transition probabilities, which are combinations of duplications and excisions
    auto transition = [&](Matrix& tr, Matrix& dup, Matrix& exc, int G)
    {
        // Initialise, as MakeLookupTables may get called multiple times
        tr = Matrix(G, Vector(G, 0));

        // Calculate transitions from genotype i -> j
        for (int j = 0; j < G; ++j)
        {
            for (int i = 0; i < G; ++i)
            {
                // An individual with i transposons can transition into an individual with j transposons
                // in many ways, where each way can be defined by the number of excisions that happen.
                // The number of duplications is then fixed, for that particular transition.
                // If i > j then we must have at least i - j excisions.
                for (int excisions = max(0, i - j), duplications = j - i + excisions; excisions <= i; ++excisions, ++duplications)
                {
                    tr[i][j] += exc[i][excisions] * dup[i][duplications];
                }
            }
        }
    };
    transition(TransitionP_A, DuplicationP_A, ExcisionP_A, A);
    transition(TransitionP_B, DuplicationP_B, ExcisionP_B, B);
}

// Swap lookup tables with buffered version (for analysis2)
void Population::SwapLookupTables()
{
    if (Swap.empty())
    {
        Swap.push_back(DuplicationP_A);
        Swap.push_back(ExcisionP_A);
        Swap.push_back(GameteP_A);
        Swap.push_back(TransitionP_A);
        Swap.push_back(DuplicationP_B);
        Swap.push_back(ExcisionP_B);
        Swap.push_back(GameteP_B);
        Swap.push_back(TransitionP_B);
    }
    else
    {
        std::swap(DuplicationP_A, Swap[0]);
        std::swap(ExcisionP_A, Swap[1]);
        std::swap(GameteP_A, Swap[2]);
        std::swap(TransitionP_A, Swap[3]);
        std::swap(DuplicationP_B, Swap[4]);
        std::swap(ExcisionP_B, Swap[5]);
        std::swap(GameteP_B, Swap[6]);
        std::swap(TransitionP_B, Swap[7]);
    }
}

// Calculate statistics
void Population::CalculateStatistics()
{
    typedef boost::accumulators::accumulator_set<double,
        boost::accumulators::stats<
            boost::accumulators::tag::weighted_mean,
            boost::accumulators::tag::lazy_weighted_variance>,
        double> accumulator_t;
    accumulator_t A_acc;
    accumulator_t B_acc;
    accumulator_t W_acc;

    // Calculate running weighted mean and variance, find the maximum genotype frequency,
    // and get marginal frequency distributions over A, B, and fitness.
    max_freq = 0;
    vector<double> A_marginal(A, 0.0);
    vector<double> B_marginal(A, 0.0);
    vector<double> W_marginal(W_list.size(), 0.0);
    int W_list_i = 0;
    for (unsigned int i = 0; i < W_lookup.size(); ++i)
    {
        int a = W_lookup[i].a, b = W_lookup[i].b;

        if (f[a][b] > 0)
        {
            A_acc(a, boost::accumulators::weight = f[a][b]);
            B_acc(b, boost::accumulators::weight = f[a][b]);
            W_acc(W_lookup[i].w, boost::accumulators::weight = f[a][b]);
        }

        if (f[a][b] > max_freq)
            max_freq = f[a][b];

        if (i > 0 && W_lookup[i].w > W_lookup[i - 1].w)
            ++W_list_i;

        A_marginal[a] += f[a][b];
        B_marginal[b] += f[a][b];
        W_marginal[W_list_i] += f[a][b];
    }

    N_A = boost::accumulators::weighted_mean(A_acc);
    V_A = boost::accumulators::lazy_weighted_variance(A_acc);
    N_B = boost::accumulators::weighted_mean(B_acc);
    V_B = boost::accumulators::lazy_weighted_variance(B_acc);
    Av_W = boost::accumulators::weighted_mean(W_acc);
    V_W = boost::accumulators::lazy_weighted_variance(W_acc);

    auto interval2 = [&](double prob, vector<double> marginal, int N, double& lower, double& upper) -> void
    {
        int li, ui;
        double lt = 0, ut = 0;

        for (li = 0; li < N; ++li)
        {
            if (lt + marginal[li] > prob)
            {
                lower = li + (prob - lt) / marginal[li];
                break;
            }
            lt += marginal[li];
        }

        for (ui = N - 1; ui >= 0; --ui)
        {
            if (ut + marginal[ui] > prob)
            {
                upper = 1 + ui - (prob - ut) / marginal[ui];
                break;
            }
            ut += marginal[ui];
        }
    };

    double W05i, W95i;
    interval2(0.05, A_marginal, A, A05, A95);
    interval2(0.05, B_marginal, B, B05, B95);
    interval2(0.05, W_marginal, W_list.size(), W05i, W95i);
    W05 = W_list[(int)W05i] + (W05i - (int)W05i) * (W_list[(int)W05i + 1] - W_list[(int)W05i]);
    W95 = W_list[(int)W95i] + (W95i - (int)W95i) * (W_list[(int)W95i + 1] - W_list[(int)W95i]);
}

// Get duplication rate of wild-type and mutant transposons
double Population::DuplicationRateWildType(int a, int b) const
{
    (void)a; (void)b;
    return P.u;
}

double Population::DuplicationRateMutant(int a, int b) const
{
    (void)a; (void)b;
    return P.u_B;
}

// Calculate average fitness
double Population::AvgW() const
{
    double s = 0;
    for (int a = 0; a < A; ++a)
        for (int b = 0; b < Minor(a); ++b)
            s += f[a][b] * w[a][b];
    return s;
}

// Calculate proportion of individuals with in-bound genotypes who are inviable
double Population::Inviables()
{
    double inviables = 0;
    for (int a = 0; a < A; ++a)
        for (int b = 0; b < Minor(a); ++b)
            if (w[a][b] <= 0)
                inviables += f[a][b];
    return inviables;
}

// Update genotypes through duplication and loss of transposons.
void Population::Update()
{
    // Transition genotypes, doing duplication and loss of transposons.
    TransitionTri(A, B, f, f, f_temp, TransitionP_A, TransitionP_B);

    // The above step may have produced some "out of bounds" genotypes; keep track of their number.
    f_other = 1 - Sum(f);

    // Determine what proportion of the population has died as a result of transposon duplication in their lifetime.
    // We don't want to include the fraction of individuals who were born with zero fitness, just those who became
    // inviable as a result of transposon duplication.
    slow = f_other + Inviables() - still;
}

// Do viability selection, multiplying the frequency of each type by its relative fitness.
// This is equivalent to killing off a fraction of each type, where the survival probability is proportional to fitness.
void Population::Selection()
{
    // Multiply each type's frequency by its relative fitness.
    double avg_w = AvgW();

    for (int a = 0; a < A; ++a)
        for (int b = 0; b < Minor(a); ++b)
            f[a][b] *= w[a][b] / avg_w;

    // There should be no out-of-bound genotypes left.
    f_other = 0;
}

// Reproduction
void Population::Reproduce()
{
    // Get gamete frequencies
    TransitionTri(A, B, f, gametes, f_temp, GameteP_A, GameteP_B);

    // Zero offspring matrix to same size as frequency vector
    offspring = ZeroTriangularish(f.size(), f.size());

    // Gamete fusion
    gamete_is.Build(gametes, A, B, P.gamete_threshold);
    for (auto K = gamete_is.begin(); K.good; ++K)
      for (auto k = gamete_is.begin_sub(K, A); k.good; ++k)
        {
            // Determine offspring's genotype
            int a = K.i + k.i, b = K.j + k.j;

            // Spawn offspring from this gametic fusion
            offspring[a][b] += gametes[K.i][K.j] * gametes[k.i][k.j];
        }

    // Enable the next generation
    swap(f, offspring);

    // Reproduction may produce some "out of bounds" genotypes; keep track of their number.
    f_other = 1 - Sum(f);

    // Determine what proportion of offspring is inviable due to out-of-bounds genotype or zero fitness.
    still = f_other + Inviables();
}

// "Special" update step for analysis mode
void Population::SpecialUpdate()
{
    const int max_special_duplications = 5;

    if (P.v > 0)
        throw std::runtime_error("SpecialUpdate cannot do excision (v > 0).\n");

    // Zero out f_temp
    for (int a = 0; a < A; ++a)
        f_temp[a].assign(Minor(a), 0.0);

    double ua = P.u;
    double ub = P.u_B;
    double ubp = P.analysis_special_prob;

    // Iterate through each genotype
    for (int ai = 0; ai < A; ++ai)
    {
        for (int bi = 0; bi < Minor(ai); ++bi)
        {
            for (int aj = ai; aj < A; ++aj)
            {
                for (int bj = bi; bj < Minor(aj); ++bj)
                {
                    for (int bs = 0; bs < min(max_special_duplications, Minor(aj) - bj); ++bs)
                    {
                        if (P.analysis_special_focus == "none")
                        {
                            // Here, the incremental B-type TE goes into [aj][bj + bs], as normal
                            f_temp[aj][bj + bs] += f[ai][bi] * PoissonPMF(ua * ai, aj - ai) * PoissonPMF(ub * bi, bj - bi) * BinomialPMF(bi, ubp, bs);
                        }
                        else if (P.analysis_special_focus == "old")
                        {
                            // Here, the incremental B-type TE is cast as an A-type TE
                            f_temp[aj + bs][bj] += f[ai][bi] * PoissonPMF(ua * ai, aj - ai) * PoissonPMF(ub * bi, bj - bi) * BinomialPMF(bi, ubp, bs);
                        }
                        else if (P.analysis_special_focus == "new")
                        {
                            // Here, all the old B-type TEs get cast as A-type TEs, and only the incremental B-type is cast as a B-type
                            f_temp[aj + bj][bs] += f[ai][bi] * PoissonPMF(ua * ai, aj - ai) * PoissonPMF(ub * bi, bj - bi) * BinomialPMF(bi, ubp, bs);
                        }
                        else if (P.analysis_special_focus == "par")
                        {
                            // Here, the mutant genome-mates of any incremental B-type transposons get cast as B-types; all others cast as A-types
                            if (bs == 0)
                                f_temp[aj + bj][0] += f[ai][bi] * PoissonPMF(ua * ai, aj - ai) * PoissonPMF(ub * bi, bj - bi) * BinomialPMF(bi, ubp, bs);
                            else
                                f_temp[aj + bs][bj] += f[ai][bi] * PoissonPMF(ua * ai, aj - ai) * PoissonPMF(ub * bi, bj - bi) * BinomialPMF(bi, ubp, bs);
                        }
                        else
                        {
                            throw std::runtime_error("Unrecognized analysis_special_focus string.\n");
                        }

                    }
                }
            }
        }
    }

    // Swap f and f_temp
    swap(f, f_temp);
}


// Print some basic statistics about the population
void Population::PrintStatistics()
{
    CalculateStatistics();

    double a_freq = N_A / (N_A + N_B), b_freq = N_B / (N_A + N_B);

    cout << "     WT             Mutant         Overall\n" << setprecision(10);
    if (N_B > 0)
        cout << "u    " << setw(15) << P.u    << setw(15) << P.u_B  << setw(15) << P.u * a_freq + P.u_B * b_freq << "\n";
    else
        cout << "u    " << setw(15) << P.u    << setw(15) << "-"    << setw(15) << P.u * a_freq + P.u_B * b_freq << "\n";
    cout << "n    " << setw(15) << N_A     << setw(15) << N_B     << setw(15) << N_A + N_B << "(A " << A05 << ", " << A95 << "); (B " << B05 << ", " << B95 << ")\n";
    cout << "freq " << setw(15) << a_freq  << setw(15) << b_freq  << setw(15) << a_freq + b_freq << "\n";
    cout << "Vn   " << setw(15) << V_A     << setw(15) << V_B     << setw(15) << V_A + V_B << "\n";
    cout << "Vn/n " << setw(15) << V_A/N_A << setw(15) << V_B/N_B << setw(15) << (V_A + V_B)/(N_A + N_B) << "\n";
    cout << "Avg. w = " << Av_W << " (var. " << V_W << ")\n";
    cout << "Total representable " << Sum(f) << ", " << f_other << " out of bounds and " << Inviables() << " inviable.\n";
    cout << "Still " << still << ", slow " << slow << ", dead " << still + slow << ".\n" << flush;
}

// Graphically show the distribution of genotypes in the population
void Population::ShowGenotypeDistribution()
{
    auto symbol = [&](double f) -> char {
        if (f == max_freq)
            return '*';
        else if (f / max_freq <= 1e-9)
            return '0';
        else
            return '0' + min(9, max(0, 9 + (int)log10(f / max_freq)));
    };

    CalculateStatistics();
    if (B == 1)
    {
        for (int a = 0; a < A; ++a)
            cout << symbol(f[a][0]);
        cout << "\n";
    }
    else
    {
        cout << "  A →\n";
        for (int b = 0; b < B; ++b)
        {
            cout << (b == 0 ? "B" : (b == 1 ? "↓" : " ")) << ' ';
            for (int a = 0; a < A - b; ++a)
                cout << symbol(f[a][b]);
            cout << "\n";
        }
    }
}

// Report progress through PrintStatistics(), optionally at a certain frequency
void Population::Report(int g, string stage, bool always)
{
    if (g % P.report_freq == 0 || always)
    {
        cout << "Generation " << g << "\n";
        if (!stage.empty()) cout << stage << "\n";
        PrintStatistics();
        cout << "\nGenotypes:\n";
        ShowGenotypeDistribution();
        cout << "\n";
    }
}

// Record a generation
void Population::RecordGeneration(Datafile& df, int g)
{
    CalculateStatistics();
    df.Set(g,  0, g);
    df.Set(g,  1, N_A);
    df.Set(g,  2, V_A);
    df.Set(g,  3, A05);
    df.Set(g,  4, A95);
    df.Set(g,  5, N_B);
    df.Set(g,  6, V_B);
    df.Set(g,  7, B05);
    df.Set(g,  8, B95);
    df.Set(g,  9, P.u);
    df.Set(g, 10, 0);
    df.Set(g, 11, P.u_B);
    df.Set(g, 12, 0);
    df.Set(g, 13, N_A * P.u + N_B * P.u_B);
    df.Set(g, 14, 0);
    df.Set(g, 15, Av_W);
    df.Set(g, 16, V_W);
    df.Set(g, 17, W05);
    df.Set(g, 18, W95);
    df.Set(g, 19, still + slow);
}

// Record a generation, analysis2 mode
void Population::Analysis2RecordGeneration(Datafile& df, int row, double g)
{
    CalculateStatistics();

    typedef boost::accumulators::accumulator_set<
        double,
        boost::accumulators::stats<boost::accumulators::tag::weighted_mean, boost::accumulators::tag::lazy_weighted_variance>,
        double> basic_acc;

    basic_acc bas_w, bas_u, bas_n, bas_a, bas_b;
    basic_acc bas_W, bas_U, bas_N, bas_A, bas_B;

    basic_acc bas_WA, bas_UA, bas_NA, bas_AA, bas_BA;
    basic_acc bas_WB, bas_UB, bas_NB, bas_AB, bas_BB;

    Covariance Cov_WU, Cov_UU, Cov_NU, Cov_AU, Cov_BU;
    Covariance Cov_WX, Cov_UX, Cov_NX, Cov_AX, Cov_BX, Cov_XX;

    basic_acc bas_R;
    Covariance Cov_nn, Cov_NN, Cov_NN_A, Cov_NN_B, Cov_aa, Cov_bb;
    Covariance Cov_an, Cov_bn, Cov_ab;

    Covariance Cov_wn, Cov_WN, Cov_WA, Cov_WB;
    Covariance Cov_AA, Cov_BB;


    for (int a = 0; a < A; ++a)
    {
        for (int b = 0; b < Minor(a); ++b)
        {
            double ww = w[a][b] / Av_W;

            // Host-level statistics
            bas_w(ww,    boost::accumulators::weight = f[a][b]);
            if (a + b > 0)
                bas_u((P.u * a + P.u_B * b) / (a + b), boost::accumulators::weight = f[a][b]);
            bas_n(a + b, boost::accumulators::weight = f[a][b]);
            bas_a(a,     boost::accumulators::weight = f[a][b]);
            bas_b(b,     boost::accumulators::weight = f[a][b]);

            double uA = DuplicationRateWildType(a, b);
            double uB = DuplicationRateMutant(a, b);

            // Transposon-level statistics...
            bas_W(ww,    boost::accumulators::weight = (a + b) * f[a][b]);
            bas_U(uA,    boost::accumulators::weight = a * f[a][b]);
            bas_U(uB,    boost::accumulators::weight = b * f[a][b]);
            bas_N(a + b, boost::accumulators::weight = (a + b) * f[a][b]);
            bas_A(a,     boost::accumulators::weight = (a + b) * f[a][b]);
            bas_B(b,     boost::accumulators::weight = (a + b) * f[a][b]);

            bas_WA(ww,    boost::accumulators::weight = a * f[a][b]);
            bas_UA(uA,    boost::accumulators::weight = a * f[a][b]);
            bas_NA(a + b, boost::accumulators::weight = a * f[a][b]);
            bas_AA(a,     boost::accumulators::weight = a * f[a][b]);
            bas_BA(b,     boost::accumulators::weight = a * f[a][b]);

            bas_WB(ww,    boost::accumulators::weight = b * f[a][b]);
            bas_UB(uB,    boost::accumulators::weight = b * f[a][b]);
            bas_NB(a + b, boost::accumulators::weight = b * f[a][b]);
            bas_AB(a,     boost::accumulators::weight = b * f[a][b]);
            bas_BB(b,     boost::accumulators::weight = b * f[a][b]);

            Cov_WU(uA,  ww,    a * f[a][b]);
            Cov_WU(uB,  ww,    b * f[a][b]);
            Cov_UU(uA,  uA,    a * f[a][b]);
            Cov_UU(uB,  uB,    b * f[a][b]);
            Cov_NU(uA,  a + b, a * f[a][b]);
            Cov_NU(uB,  a + b, b * f[a][b]);
            Cov_AU(uA,  a,     a * f[a][b]);
            Cov_AU(uB,  a,     b * f[a][b]);
            Cov_BU(uA,  b,     a * f[a][b]);
            Cov_BU(uB,  b,     b * f[a][b]);

            Cov_WX(P.u,   ww,    a * f[a][b]);
            Cov_WX(P.u_B, ww,    b * f[a][b]);
            Cov_UX(P.u,   uA,    a * f[a][b]);
            Cov_UX(P.u_B, uB,    b * f[a][b]);
            Cov_NX(P.u,   a + b, a * f[a][b]);
            Cov_NX(P.u_B, a + b, b * f[a][b]);
            Cov_AX(P.u,   a,     a * f[a][b]);
            Cov_AX(P.u_B, a,     b * f[a][b]);
            Cov_BX(P.u,   b,     a * f[a][b]);
            Cov_BX(P.u_B, b,     b * f[a][b]);
            Cov_XX(P.u,   P.u,   a * f[a][b]);
            Cov_XX(P.u_B, P.u_B, b * f[a][b]);

            if (a + b > 0) 
                bas_R(double(b) / double(a + b), boost::accumulators::weight = b * f[a][b]);

            Cov_nn  (a + b, a + b, f[a][b]);
            Cov_NN  (a + b, a + b, (a + b) * f[a][b]);
            Cov_NN_A(a + b, a + b, a       * f[a][b]);
            Cov_NN_B(a + b, a + b, b       * f[a][b]);
            Cov_aa  (a,     a,     f[a][b]);
            Cov_bb  (b,     b,     f[a][b]);
            Cov_an  (a,     a + b, f[a][b]);
            Cov_bn  (b,     a + b, f[a][b]);
            Cov_ab  (a,     b,     f[a][b]);

            Cov_wn  (ww,    a + b, f[a][b]);
            Cov_WN  (ww,    a + b, (a + b) * f[a][b]);
            Cov_WA  (ww,    a,     a * f[a][b]);
            Cov_WB  (ww,    b,     b * f[a][b]);

            Cov_AA  (a,     a,     a * f[a][b]);
            Cov_BB  (b,     b,     b * f[a][b]);
        }
    }

    df.Set(row, 0, g);

    // avg w, u, n, a, b
    saved_b = boost::accumulators::weighted_mean(bas_b);
    df.Set(row, 1, boost::accumulators::weighted_mean(bas_w));
    df.Set(row, 2, boost::accumulators::weighted_mean(bas_u));
    df.Set(row, 3, boost::accumulators::weighted_mean(bas_n));
    df.Set(row, 4, boost::accumulators::weighted_mean(bas_a));
    df.Set(row, 5, boost::accumulators::weighted_mean(bas_b));

    // avg W, U, N, A, B
    df.Set(row, 6, boost::accumulators::weighted_mean(bas_W));
    df.Set(row, 7, boost::accumulators::weighted_mean(bas_U));
    df.Set(row, 8, boost::accumulators::weighted_mean(bas_N));
    df.Set(row, 9, boost::accumulators::weighted_mean(bas_A));
    df.Set(row, 10, boost::accumulators::weighted_mean(bas_B));

    // avg W, U, N, A, B for A only
    df.Set(row, 11, boost::accumulators::weighted_mean(bas_WA));
    df.Set(row, 12, boost::accumulators::weighted_mean(bas_UA));
    df.Set(row, 13, boost::accumulators::weighted_mean(bas_NA));
    df.Set(row, 14, boost::accumulators::weighted_mean(bas_AA));
    df.Set(row, 15, boost::accumulators::weighted_mean(bas_BA));

    // avg W, U, N, A, B for B only                
    df.Set(row, 16, boost::accumulators::weighted_mean(bas_WB));
    df.Set(row, 17, boost::accumulators::weighted_mean(bas_UB));
    df.Set(row, 18, boost::accumulators::weighted_mean(bas_NB));
    df.Set(row, 19, boost::accumulators::weighted_mean(bas_AB));
    df.Set(row, 20, boost::accumulators::weighted_mean(bas_BB));

    // slope of W on U, U on U, N on U, A on U, B on U
    saved_sWU = Cov_WU.Cov() / Cov_UU.Cov();
    df.Set(row, 21, saved_sWU);
    df.Set(row, 22, Cov_UU.Cov() / Cov_UU.Cov());
    df.Set(row, 23, Cov_NU.Cov() / Cov_UU.Cov());
    df.Set(row, 24, Cov_AU.Cov() / Cov_UU.Cov());
    df.Set(row, 25, Cov_BU.Cov() / Cov_UU.Cov());

    // As above on X
    df.Set(row, 26, Cov_WX.Cov() / Cov_XX.Cov());
    df.Set(row, 27, Cov_UX.Cov() / Cov_XX.Cov());
    df.Set(row, 28, Cov_NX.Cov() / Cov_XX.Cov());
    df.Set(row, 29, Cov_AX.Cov() / Cov_XX.Cov());
    df.Set(row, 30, Cov_BX.Cov() / Cov_XX.Cov());
    
    // relatedness and variance of copy number from overall/WT/mutant TE perspective
    df.Set(row, 31, boost::accumulators::weighted_mean(bas_R));
    df.Set(row, 32, Cov_nn.Cov());
    df.Set(row, 33, Cov_NN.Cov());
    df.Set(row, 34, Cov_NN_A.Cov());
    df.Set(row, 35, Cov_NN_B.Cov());

    // more variances / covariances
    df.Set(row, 36, Cov_aa.Cov());
    df.Set(row, 37, Cov_bb.Cov());
    df.Set(row, 38, Cov_an.Cov());
    df.Set(row, 39, Cov_bn.Cov());
    df.Set(row, 40, Cov_ab.Cov());

    // fitness
    df.Set(row, 41, Av_W);
    df.Set(row, 42, Cov_wn.Cov() / Cov_nn.Cov());
    df.Set(row, 43, Cov_WN.Cov() / Cov_NN.Cov());
    df.Set(row, 44, Cov_WA.Cov() / Cov_AA.Cov());
    df.Set(row, 45, Cov_WB.Cov() / Cov_BB.Cov());

    // reporting duplication rates
    df.Set(row, 46, P.u);
    df.Set(row, 47, P.u_B);
}

// Initialise and run simulation with one transposon type until a certain condition is met.
// Parameters:
//  n_generations   maximum number of generations for which to run the simulation.
//  epsilon         simulation terminates early if distance between two consecutive population
//                  frequency matrices is less than epsilon. can pass negative number to override.
//  verbose         if true, reports on progress every P.report_freq generations.
//  df              if non-null, records statistics to df every generation.
//  record_start    which datafile row to start recording in.
void Population::SingleRun(int n_generations, double epsilon, bool verbose, Datafile* df, int record_start)
{
    // Initialize population with wild-type transposons
    B = 1; Init();

    // Record prior to first generation
    if (df)
        RecordGeneration(*df, record_start);
    Report(0, "Start of single-transposon run.", true);

    // Simulation
    bool quit = false;
    for (int g = 0; !quit; ++g)
    {
        // 1. UPDATE INDIVIDUALS
        Update();
        if (verbose)
            Report(g + 1, "After updating.");

        // 2. RECORDING
        if (df)
            RecordGeneration(*df, record_start + g + 1);

        // 3. SELECTION
        Selection();
        if (verbose)
            Report(g + 1, "After selection.");

        // Check for exit due to epsilon: at this point, offspring matrix contains f from previous generation.
        if (epsilon > 0 && Distance(f, offspring) < epsilon)
        {
            Report(g + 1, "Epsilon reached after " + to_string(g + 1) + " generations.", true);
            quit = true;
        }

        // Check for exit due to end of generations
        if (n_generations > 0 && g == n_generations - 1)
        {
            Report(g + 1, "End reached after " + to_string(g + 1) + " generations.", true);
            quit = true;
        }

        // 4. REPRODUCTION
        Reproduce();
        if (verbose)
            Report(g + 1, "After reproduction.");
    }
}

// Introduce mutants to a wild-type-only population.
void Population::IntroduceMutants(double mutant_freq, int n_B)
{
    if (B > 1)
    {
        throw std::runtime_error("Cannot reintroduce mutants.\n");
    }

    // Enable genotypes with mutant transposons
    if (n_B > 0)
        B = n_B;
    else
        B = A;
    MakeLookupTables();

    // Introduce mutants
    for (int a = 1; a < A; ++a)
    {
        double frequency = f[a][0];
        for (int m = 0; m <= min(a, B - 1); ++m)
        {
            f[a - m][m] = frequency * BinomialPMF(a, mutant_freq, m);
        }
    }
}

// Consolidate all transposons to wild-type.
void Population::Consolidate()
{
    for (int a = 0; a < A; ++a)
    {
        for (int b = 1; b < Minor(a); ++b)
        {
            f[a + b][0] += f[a][b];
            f[a][b] = 0;
        }
    }
    B = 1;
    MakeLookupTables();
}

// Run simulation
void Population::Run()
{
    vector<string> datafile_columns = { "g", "n", "Vn", "n05", "n95", "m", "Vm", "m05", "m95", "u", "Vu", "U", "VU", "a", "Va", "w", "Vw", "w05", "w95", "d" };

    // Single run, for a certain number of generations or until a desired tolerance is reached.
    if (P.mode == "static")
    {
        Datafile df(1, P.generations + 1, datafile_columns.size(), datafile_columns);
        SingleRun(P.generations, P.epsilon, true, &df, 0);
        df.ShrinkToFit();
        df.Save(P.fileout);
    }

    // Introduce a mutant genotype into a population at equilibrium.
    else if (P.mode == "invasion")
    {
        Datafile df(1, P.generations + 1, datafile_columns.size(), datafile_columns);

        // Run with the wild-type transposon until epsilon is reached
        SingleRun(-1, P.epsilon, false);

        // Enable genotypes with mutant transposons
        IntroduceMutants(P.invasion_seed);
        cout << "Introduced mutants.\n";
        PrintStatistics();

        // Record prior to first generation
        RecordGeneration(df, 0);

        // Proceed until invasion_threshold or generations is reached
        bool invasion = false;
        cout << "Testing invasion...\n";
        int i;
        for (i = 0; i < P.generations; ++i)
        {
            Update();

            // Record and check for end
            RecordGeneration(df, i + 1);
            if (N_B / (N_A + N_B) > P.invasion_threshold)
            {
                invasion = true;
                break;
            }

            Selection();
            Reproduce();
            Report(i + 1, "After reproduction.");
        }

        // Print results . . .
        if (invasion)
            cout << "Mutant allele invaded after " << i << " generations.\n";
        else
            cout << "No invasion after " << i << " generations.\n";
        PrintStatistics();

        df.ShrinkToFit();
        df.Save(P.fileout);

        // Set n_0 to current n, as a guess for any subsequent runs
        P.n_0 = N_A + N_B;
    }

    // Sequentially introduce mutants into a population at regular intervals.
    else if (P.mode == "sequential")
    {
        // Check parameters
        if (P.sequential_g.size() != P.sequential_u.size())
            throw std::runtime_error("sequential_g and sequential_u must be same size!");
        if (P.sequential_g.empty())
            throw std::runtime_error("sequential_g and sequential_u must be populated!");
        if (!P.sequential_n_max.empty() && P.sequential_n_max.size() != P.sequential_g.size())
            throw std::runtime_error("sequential_n_max must be size zero or equal to size of sequential_g and sequential_u.");
        P.u_B = 0;

        int total_generations = 1 + accumulate(P.sequential_g.begin(), P.sequential_g.end(), 0);
        Datafile df(1, total_generations, datafile_columns.size(), datafile_columns);

        // Run with the wild-type transposon for specified number of generations
        P.u = P.sequential_u[0];
        if (!P.sequential_n_max.empty())
            A = P.sequential_n_max[0];
        SingleRun(P.sequential_g[0], 0, true, &df, 0);
        int g = P.sequential_g[0] + 1;

        // Start invasion phase
        for (unsigned int s = 1; s < P.sequential_u.size(); ++s)
        {
            cout << "Invasion " << s << ".\n";

            P.u_B = P.sequential_u[s];
            if (!P.sequential_n_max.empty())
                A = P.sequential_n_max[s];

            cout << "Introducing mutants...\n" << flush;
            IntroduceMutants(P.invasion_seed);

            cout << "Running invasion... " << flush;
            for (int i = 0; i < P.sequential_g[s]; ++i)
            {
                Update();
                RecordGeneration(df, g);
                Selection();
                Reproduce();
                if (i % P.report_freq == 0)
                    cout << "." << flush;
                ++g;
            }
            cout << " OK.\n";
            Report(g, "Sequential stage.", true);

            cout << "Consolidating genotypes.\n";
            Consolidate();
            ShowGenotypeDistribution();
            cout << "\n";

            // Reset for next invasion
            P.u = P.u_B;
        }

        df.ShrinkToFit();
        df.Save(P.fileout);
    }

    // Investigate, for a very rare mutant, various "gene's-eye" statistics.
    else if (P.mode == "analysis")
    {
        vector<string> columns = { "g",
                                   "n", "n_x", "n_xa", "n_xb",
                                   "na", "na_x", "na_xa", "na_xb",
                                   "nb", "nb_x", "nb_xa", "nb_xb",
                                   "w", "w_x", "w_xa", "w_xb",
                                   "W", "W_x", "W_xa", "W_xb",
                                   "nb_withb" };

        Datafile df(1, P.analysis_gens + 2, columns.size(), columns);

        // First run to equilibrium with wild-type.
        SingleRun(-1, P.epsilon, false);

        // Lambda for collecting one generation's worth of information
        auto record = [&](Datafile& df, int row, int g)
        {
            CalculateStatistics();

            double n = 0, n_x = 0, n_xa = 0, n_xb = 0;      // Average copy number + experienced by all / wildtype / mutant transposons
            double na = 0, na_x = 0, na_xa = 0, na_xb = 0;  // Average wildtype number + experienced by all / wildtype / mutant transposons
            double nb = 0, nb_x = 0, nb_xa = 0, nb_xb = 0;  // Average mutant number + experienced by all / wildtype / mutant transposons
            double w_ = 0, w_x = 0, w_xa = 0, w_xb = 0;      // Average host fitness + experienced by all / wildtype / mutant transposons

            double s = 0, sa = 0, sb = 0;     // Frequency of experiences
            double f_withb = 0;

            for (int a = 0; a < A; ++a)
            {
                for (int b = 0; b < Minor(a); ++b)
                {
                    // Copy number
                    n    += (a + b) * f[a][b];
                    n_x  += (a + b) * f[a][b] * (a + b);
                    n_xa += (a + b) * f[a][b] * a;
                    n_xb += (a + b) * f[a][b] * b;

                    // Wildtype number
                    na    += a * f[a][b];
                    na_x  += a * f[a][b] * (a + b);
                    na_xa += a * f[a][b] * a;
                    na_xb += a * f[a][b] * b;

                    // Mutant number
                    nb    += b * f[a][b];
                    nb_x  += b * f[a][b] * (a + b);
                    nb_xa += b * f[a][b] * a;
                    nb_xb += b * f[a][b] * b;

                    // Host fitness
                    w_   += w[a][b] * f[a][b];
                    w_x  += w[a][b] * f[a][b] * (a + b);
                    w_xa += w[a][b] * f[a][b] * a;
                    w_xb += w[a][b] * f[a][b] * b;

                    // Frequency of experiences
                    s  += f[a][b] * (a + b);
                    sa += f[a][b] * a;
                    sb += f[a][b] * b;

                    if (b > 0)
                    {
                        f_withb += f[a][b];
                    }
                }
            }

            df.Set(row, 0, g);

            df.Set(row, 1, n);
            df.Set(row, 2, n_x / s);
            df.Set(row, 3, n_xa / sa);
            df.Set(row, 4, n_xb / sb);

            df.Set(row, 5, na);
            df.Set(row, 6, na_x / s);
            df.Set(row, 7, na_xa / sa);
            df.Set(row, 8, na_xb / sb);

            df.Set(row, 9, nb);
            df.Set(row, 10, nb_x / s);
            df.Set(row, 11, nb_xa / sa);
            df.Set(row, 12, nb_xb / sb);

            df.Set(row, 13, w_);
            df.Set(row, 14, w_x / s);
            df.Set(row, 15, w_xa / sa);
            df.Set(row, 16, w_xb / sb);

            df.Set(row, 17, (w_ / Av_W));
            df.Set(row, 18, (w_x / Av_W) / s);
            df.Set(row, 19, (w_xa / Av_W) / sa);
            df.Set(row, 20, (w_xb / Av_W) / sb);

            df.Set(row, 21, nb / f_withb);
        };

        // Record first observation in absence of mutants.
        record(df, 0, -1);

        // Now, introduce mutant at low frequency.
        cout << "Introducing mutants:\n" << flush;
        IntroduceMutants(P.analysis_seed, P.analysis_nb_max);
        PrintStatistics();

        cout << "\nRunning generations...\n" << flush;

        // Allow to run for a small number of generations, recording associations, etc.
        for (int g = 0; g <= P.analysis_gens; ++g)
        {
            if (g % 10 == 0)
                cout << "." << flush;

            //record(df, g + 1, g);                

            Update();
            record(df, g + 1, g);
            Selection();
            Reproduce();
        }

        cout << "\n" << flush;

        df.ShrinkToFit();
        df.Save(P.fileout);
    }

    // Investigate various other "gene's-eye" statistics.
    else if (P.mode == "analysis2")
    {
        // throughout, n = a+b; a, b are transposon numbers of wild-type and mutant type; s is slope; uppercase means over TEs, lowercase over hosts
        // note that w is relative fitness below -- abs_w is absolute fitness
        vector<string> columns = { "g",                                     // generation
                                   "w", "u", "n", "a", "b",                 // mean of w, u, n, a, b over hosts
                                   "W", "U", "N", "A", "B",                 // mean of w, u, n, a, b over TEs
                                   "W_A", "U_A", "N_A", "A_A", "B_A",       // mean of w, u, n, a, b over wild-type TEs
                                   "W_B", "U_B", "N_B", "A_B", "B_B",       // mean of w, u, n, a, b over mutant TEs
                                   "sWU", "sUU", "sNU", "sAU", "sBU",       // slope of w, u, n, a, b on realized duplication rate over TEs
                                   "sWX", "sUX", "sNX", "sAX", "sBX",       // slope of w, u, n, a, b on activity level over TEs
                                   "r", "Vn", "VN", "VN_A", "VN_B",         // mean relatedness over mutant TEs; variance in n over hosts, TEs, wild-type TEs, mutant TEs
                                   "Va", "Vb", "Van", "Vbn", "Vab",         // variance and covarince in number of wild-type, mutant TEs over hosts
                                   "abs_w", "swn", "sWN", "sWN_A", "sWN_B", // mean absolute fitness over hosts; slope of w on n over hosts; slope of W on N over TEs
                                   "uA", "uB"                               // u for wild-type, mutant TEs
                                   };

        Datafile df(1, 1 + 3 * (P.analysis_gens + 1), columns.size(), columns);

        // First, run to equilibrium with wild-type.
        SingleRun(P.analysis_gens0, P.epsilon, false);

        // Save the last values
        double comp_b = 0;

        // Record first observation in absence of mutants.
        Analysis2RecordGeneration(df, 0, 0);

        // Now, introduce mutant at low frequency.
        cout << "Introducing mutants:\n" << flush;
        IntroduceMutants(P.analysis_seed, P.analysis_nb_max);
        PrintStatistics();

        // Saved genotypes
        vector<string> genotype_columns = { "g", "A", "B", "f" };
        unsigned int n_genotype_row = 0;
        if (!P.fileout_genotypes.empty())
        {
            for (int a = 0; a < A; ++a)
                n_genotype_row += Minor(a);
            n_genotype_row *= 1 + P.analysis_gens / P.analysis_genotype_report_freq;
        }
        Datafile dfg(1, n_genotype_row, genotype_columns.size(), genotype_columns);

        cout << "\nRunning generations...\n" << flush;

        // Allow to run for a small number of generations, recording associations, etc.
        unsigned int genotype_row = 0;
        unsigned int swap_number = 0; // for analysis_u_change_gens
        for (int g = 0; g < P.analysis_gens; ++g)
        {
            if (g % 10 == 0)
                cout << "." << flush;

            if (!P.fileout_genotypes.empty() && g % P.analysis_genotype_report_freq == 0)
            {
                for (int a = 0; a < A; ++a)
                {
                    for (int b = 0; b < Minor(a); ++b)
                    {
                        dfg.Set(genotype_row, 0, g + 1);
                        dfg.Set(genotype_row, 1, a);
                        dfg.Set(genotype_row, 2, b);
                        dfg.Set(genotype_row, 3, f[a][b]);
                        ++genotype_row;
                    }
                }
            }

            if (find(P.analysis_u_change_gens.begin(), P.analysis_u_change_gens.end(), g + 1) != P.analysis_u_change_gens.end())
            {
                ++swap_number;
                if (swap_number == 1) {
                    // Set u_B to u_B + amt
                    P.u_B += P.analysis_u_change_amt;
                    SwapLookupTables();
                    MakeLookupTables();
                } else if (swap_number % 2 == 0) {
                    // Restore u_B to original u_B
                    P.u_B -= P.analysis_u_change_amt;
                    ConstructFitness();
                    SwapLookupTables();
                } else {
                    // Restore u_B to original u_B + amt
                    P.u_B += P.analysis_u_change_amt;
                    ConstructFitness();
                    SwapLookupTables();
                }
            }

            // if (g == 1) -- used to compare to introduction, but now just compare to previous generation.
            comp_b = saved_b;

            Analysis2RecordGeneration(df, (g * 3) + 1, g + 1.1);
            if (g + 1 != P.analysis_special_gen)
                Update();
            else
                SpecialUpdate();
            Analysis2RecordGeneration(df, (g * 3) + 2, g + 1.2);
            Selection();
            Analysis2RecordGeneration(df, (g * 3) + 3, g + 1.3);
            Reproduce();
        }

        cout << (comp_b > saved_b ? "\nNo invasion.\n" : "\nInvasion.\n");
        cout << "final sWU: " << saved_sWU << "\n" << flush;

        df.ShrinkToFit();
        df.Save(P.fileout);
        if (!P.fileout_genotypes.empty())
            dfg.Save(P.fileout_genotypes);
    }


    // How much do equilibrium frequencies differ from a Poisson distribution over a range of duplication rates?
    else if (P.mode == "poisson")
    {
        vector<double> u_levels;
        int n_u_levels = int((P.u_B - P.u) / P.poisson_u_step + 0.5);
        for (int u = 0; u <= n_u_levels; ++u)
            u_levels.push_back(P.u + u * P.poisson_u_step);

        vector<string> columns = { "i", "u", "n", "diff" };

        Datafile df(1, u_levels.size(), columns.size(), columns);

        // Initialize population with wild-type transposons
        B = 1; Init();

        for (unsigned int i = 0; i < u_levels.size(); ++i)
        {
            P.u = u_levels[i];
            MakeLookupTables();

            // Run until epsilon is reached
            bool quit = false;
            while (!quit)
            {
                Update();
                Selection();
                if (Distance(f, offspring) < P.epsilon)
                    quit = true;
                Reproduce();
            }

            // note: findings are that: quitting after Selection makes the results least like a Poisson distribution;
            // quitting after Update is the best match to the Poisson distribution (at least over some ranges),
            // and makes the fit better for lower ranges of u; quitting after Reproduce is in the middle.
            // Update();
            // Selection();
            CalculateStatistics();
            
            vector<double> freq(P.n_max, 0);
            for (int i = 0; i < P.n_max; ++i)
                freq[i] = f[i][0];

            vector<double> poisson_distribution(P.n_max, 0);
            for (int i = 0; i < P.n_max; ++i)
                poisson_distribution[i] = PoissonPMF(N_A, i);
            
            double diff = Distance(freq, poisson_distribution);

            df.Set(i, 0, i);
            df.Set(i, 1, u_levels[i]);
            df.Set(i, 2, N_A);
            df.Set(i, 3, diff);

            cout << i << ". For u = " << u_levels[i] << ", diff = " << diff << "\n" << flush;
            cout << "Sums: ∑ freq = " << Sum(freq) << ", ∑ poisson = " << Sum(poisson_distribution) << "\n";
        }

        df.ShrinkToFit();
        df.Save(P.fileout);
    }

    else
    {
        cout << "Unrecognised simulation mode " << P.mode << ".\n";
    }
}
