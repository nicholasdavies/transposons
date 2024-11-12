#include "population_dna.h"
#include <csignal>

PopulationDNA::PopulationDNA()
 : Population(),
   TransitionP()
{ }

// Populate lookup tables
void PopulationDNA::MakeLookupTables()
{
    ConstructFitness();

    // Calculate probabilities of gamete formation.
    auto gamete = [&](Matrix& gam, int G)
    {
        for (int g = 0; g < G; ++g)
            for (int e = 0; e <= g; ++e)
                gam[g][e] = BinomialPMF(g, 0.5, e);
    };
    gamete(GameteP_A, A);
    gamete(GameteP_B, B);

    // The model assumes the following sequence for transposon events:
    // 1) Each transposon produces transposases through a Poisson process at rate equal to the average duplication rate.
    // 2) Each transposon has a fixed probability v of being excised/lost.
    // 3) The transposases do their work, creating a new transposon each regardless of whether their producing transposon has been lost.
    // In this way, the number of duplications is D ~ Poisson(n * u) and the number of excisions is E ~ Binomial(n, v).

    // Precalculate excision probabilities
    int G = max(A, B);
    Matrix exc(G, Vector(G, 0));
    for (int g = 0; g < G; ++g)
        for (int e = 0; e <= g; ++e)
            exc[g][e] = BinomialPMF(g, P.v, e);

    Matrix dup(G, Vector(G, 0));

    // Calculate transition probabilities, which are combinations of duplications and excisions, for each genotype
    cout << "Calculating transition probabilities .... \n" << string(A, '_') << "\n";
    TransitionP.clear();
    for (int a = 0; a < A; ++a)
    {
        cout << "*" << flush;
        TransitionP.push_back(vector<SparseMatrix>());
        for (int b = 0; b < Minor(a); ++b)
        {
            // Work out average duplication rate for this genotype
            double u = a + b > 0 ? (a * P.u + b * P.u_B) / (a + b) : 0;

            // Initialise transition matrix
            Matrix tr;
            if (B == 1) tr = Matrix(A, Vector(1, 0.0));
            else        tr = ZeroTriangularish(A, B);

            // Precalculate duplication probabilities
            for (int g = 0; g < G; ++g)
                for (int e = 0; e < G; ++e)
                    dup[g][e] = PoissonPMF(g * u, e);

            // Calculate transitions from genotype (a,b) -> (i,j)
            double marginalA = 1.0, marginalB = 1.0;    // marginal probabilities of transition from a -> i, and b -> j, respectively
            for (int i = 0; i < A; ++i)
            {
                for (int j = 0; j < Minor(i); ++j)
                {
                    // An individual with a (b) transposons can transition into an individual with i (j) transposons
                    // in many ways, where each way can be defined by the number of excisions that happen.
                    // The number of duplications is then fixed, for that particular transition.
                    // If a > i (b > j) then we must have at least a - i (b - j) excisions.
                    for (int excisionsA = max(0, a - i), duplicationsA = i - a + excisionsA; excisionsA <= a; ++excisionsA, ++duplicationsA)
                    {
                        marginalA = exc[a][excisionsA] * dup[a][duplicationsA];

                        if (marginalA < P.marginal_threshold) break;
                        
                        for (int excisionsB = max(0, b - j), duplicationsB = j - b + excisionsB; excisionsB <= b; ++excisionsB, ++duplicationsB)
                        {
                            marginalB = exc[b][excisionsB] * dup[b][duplicationsB];
                            tr[i][j] += marginalA * marginalB;

                            if (marginalA * marginalB < P.marginal_threshold) break;
                        }
                    }
                }
            }

            // Insert result into collection of transition probabilities
            TransitionP.back().push_back(SparseMatrix(tr, P.transition_threshold));
        }
    }
    cout << "\n" << flush;
}

// Swap lookup tables with buffered version (for analysis2)
void PopulationDNA::SwapLookupTables()
{
    Population::SwapLookupTables();
    if (TransitionP_Swap.empty())
    {
        TransitionP_Swap = TransitionP;
    }
    else
    {
        std::swap(TransitionP, TransitionP_Swap);
    }
}

// Get duplication rate of wild-type and mutant transposons
double PopulationDNA::DuplicationRateWildType(int a, int b) const
{
    return (a + b > 0) ? (P.u * a + P.u_B * b) / (a + b) : 0;
}

double PopulationDNA::DuplicationRateMutant(int a, int b) const
{
    return (a + b > 0) ? (P.u * a + P.u_B * b) / (a + b) : 0;
}

// Update genotypes through duplication and loss of transposons.
void PopulationDNA::Update()
{
    // Transition genotypes, doing duplication and loss of transposons.
    f_temp = ZeroTriangularish(A, A);
    for (int a = 0; a < A; ++a)
    {
        for (int b = 0; b < Minor(a); ++b)
        {
            // Move frequency from this genotype to others
            if (f[a][b] > P.transition_threshold)
                for (auto& e : TransitionP[a][b].d)
                    f_temp[e.i][e.j] += e.x * f[a][b];
        }
    }
    f = f_temp;

    // The above step may have produced some "out of bounds" genotypes; keep track of their number.
    f_other = 1 - Sum(f);

    // Determine what proportion of the population has died as a result of transposon duplication in their lifetime.
    // We don't want to include the fraction of individuals who were born with zero fitness, just those who became
    // inviable as a result of transposon duplication.
    slow = f_other + Inviables() - still;
}

// // "Special" low-probability incremental duplication step for analysis mode
// void PopulationDNA::SpecialDuplication()
// {
//     // Zero out f_temp
//     for (int a = 0; a < A; ++a)
//         f_temp[a].assign(Minor(a), 0.0);

//     // Iterate through each genotype
//     for (int a = 0; a < A; ++a)
//     {
//         for (int b = 0; b < Minor(a); ++b)
//         {
//             int a0, b0, a1A, b1A, a1B, b1B; // Genotypes to add frequency of (0) no duplication (1A) one duplication of wild type (1B) one duplication of mutant type

//             if (P.analysis_special_focus == "none")
//             {
//                 // Here, the duplicated TE goes into [a + 1][b] if wildtype and [a][b + 1] if mutant, as normal
//                 a0 = a;
//                 b0 = b;
//                 a1A = a + 1;
//                 b1A = b;
//                 a1B = a;
//                 b1B = b + 1;
//             }
//             else if (P.analysis_special_focus == "old")
//             {
//                 // Here, the duplicated TE is cast as an A-type TE (irrespective of whether it is wild-type or mutant)
//                 a0 = a;
//                 b0 = b;
//                 a1A = a + 1;
//                 b1A = b;
//                 a1B = a + 1;
//                 b1B = b;
//             }
//             else if (P.analysis_special_focus == "new")
//             {
//                 // Here, all the old B-type TEs get cast as A-type TEs, and only the duplicated B-type is cast as a B-type, only if it is a mutant
//                 a0 = a + b;
//                 b0 = 0;
//                 a1A = a + b + 1;
//                 b1A = 0;
//                 a1B = a + b;
//                 b1B = 1;
//             }
//             else
//             {
//                 throw std::runtime_error("Unrecognized analysis_special_focus string.\n");
//             }

//             // Enforce bounds
//             a0  = min(A - 1, a0);
//             b0  = min(Minor(a0) - 1, b0);
//             a1A = min(A - 1, a1A);
//             b1A = min(Minor(a1A) - 1, b1A);
//             a1B = min(A - 1, a1B);
//             b1B = min(Minor(a1B) - 1, b1B);

//             // Update genotypes
//             double fA = 0.0;
//             if (a + b > 0) 
//                 fA = a / double(a + b);

//             f_temp[a0][b0]   += f[a][b] * (1.0 - P.analysis_special_prob * b);      // Frequency of f[a][b] genotypes not changing
//             f_temp[a1A][b1A] += f[a][b] * P.analysis_special_prob * b * fA;         // Frequency of f[a][b] genotypes adding one extra wild-type TE
//             f_temp[a1B][b1B] += f[a][b] * P.analysis_special_prob * b * (1. - fA);  // Frequency of f[a][b] genotypes adding one extra mutant TE
//         }
//     }

//     // Swap f and f_temp
//     swap(f, f_temp);
// }

// "Special" update step for analysis mode
void PopulationDNA::SpecialUpdate()
{
    const int max_special_duplications = 5;

    if (P.v > 0)
        throw std::runtime_error("SpecialUpdate cannot do excision (v > 0).\n");

    if (P.u != P.u_B)
        throw std::runtime_error("SpecialUpdate requires u == u_B.\n");

    // Zero out f_temp
    for (int a = 0; a < A; ++a)
        f_temp[a].assign(Minor(a), 0.0);

    double u = P.u;
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
                    int an = aj - ai;
                    int bn = bj - bi;
                    double pnd = PoissonPMF((ai + bi) * u, an + bn); // Probability of an + bn normal duplications
                    double pa = 0;
                    if (ai + bi > 0)
                        pa = double(ai) / (ai + bi);
                    double pna = BinomialPMF(an + bn, pa, an); // Probability that an of the an + bn normal duplications are type A

                    for (int sn = 0; sn < min(max_special_duplications, Minor(aj) - bj); ++sn)
                    {
                        double psd = BinomialPMF(bi, ubp, sn); // Probability of sn special duplications
                        for (int as = 0, bs = sn; as <= sn; ++as, --bs)
                        {
                            if (P.analysis_special_focus == "none")
                            {
                                // Here, the incremental TE goes into A or B according to its type
                                f_temp[aj + as][bj + bs] += f[ai][bi] * pnd * pna * psd * BinomialPMF(sn, pa, as);
                            }
                            else if (P.analysis_special_focus == "old")
                            {
                                // Here, the incremental TE is cast as an A-type TE
                                f_temp[aj + as + bs][bj] += f[ai][bi] * pnd * pna * psd * BinomialPMF(sn, pa, as);
                            }
                            else if (P.analysis_special_focus == "new")
                            {
                                // Here, all the old B-type TEs get cast as A-type TEs, and only an incremental B-type TE is cast as a B-type
                                f_temp[aj + as + bj][bs] += f[ai][bi] * pnd * pna * psd * BinomialPMF(sn, pa, as);
                            }
                            else if (P.analysis_special_focus == "par")
                            {
                                // Here, the mutant genome-mates of any incremental transposons get cast as B-types; all others cast as A-types
                                if (sn == 0)
                                    f_temp[aj + bj][0] += f[ai][bi] * pnd * pna * psd * BinomialPMF(sn, pa, as);
                                else
                                    f_temp[aj + sn][bj] += f[ai][bi] * pnd * pna * psd * BinomialPMF(sn, pa, as);
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
    }

    // Swap f and f_temp
    swap(f, f_temp);
}

