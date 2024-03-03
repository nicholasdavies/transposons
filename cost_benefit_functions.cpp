#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <cmath>
#include <Rcpp.h>

// Helper functions
double dpois(int x, double lambda)
{
    if (x < 0)
        return 0.0;
    if (lambda == 0)
        return x == 0 ? 1.0 : 0;
    double y = x * log(lambda) - lambda - lgamma(x + 1);
    return exp(y);
}

double dbinom(int x, int n, double p)
{
    if (x > n)
        return 0.0;
    double y = lgamma(n+1) + x*log(p) + (n-x)*log(1-p) - lgamma(x+1) - lgamma(n-x+1);
    if (std::isnan(exp(y))) {
        Rcpp::Rcout << y << "!!!\n";
    }
    return exp(y);
}

double dquasi(int x, double u, double phi = 1.0)
{
    if (x == 0)
        return 1 - phi;
    return phi * (1 - u) * pow(u, x - 1);
}

double fitness(double n, double cost)
{
    return pow(1.0 - cost, 0.5 * n * n);
}

// Wrapper for a 3D matrix
class Matrix3d
{
public:
    int n_rows;
    int n_cols;
    int n_slices;
    SEXP S;
    
    Matrix3d(SEXP s)
     : S(s)
    {
        SEXP dims = Rf_getAttrib(S, R_DimSymbol);
        if (TYPEOF(dims) != INTSXP)
            Rf_error("dims should be type INTSXP");
        n_rows = INTEGER(dims)[0];
        n_cols = INTEGER(dims)[1];
        n_slices = INTEGER(dims)[2];
    }
    
    double& operator()(int x, int y, int z)
    {
        --x; --y; --z;
        int i = x + y * n_rows + z * n_rows * n_cols;
        return REAL(S)[i];
    }
};

// Run generations of a focal-mutant-lineage simulation with no selection
// xd       an AxBxG matrix to hold genotype frequencies after duplication; should be zeros for generation g
// xr       an AxBxG matrix to hold genotype frequencies after recombination; should be zeros for generation g
// g        the generation to run; 0-based, but must be > 0 (because information from g-1 is read from xr) and < G
// u        duplication rate
// p        extra duplication probability
// omega    mean number of type-A TEs in outcross gametes
// mode     "pub" or "pri"
// focus    "none", "new", or "old" -- as in deterministic simulation
// gametes  if provided, use this for distribution of TE number in outcross 
//              gametes rather than Poisson(omega); must be either size 0 or be
//              size max(A, B)
// [[Rcpp::export]]
Rcpp::List mutant_lineage(
    SEXP xdS, SEXP xrS, unsigned int g, 
    double u, double p, double omega, 
    std::string mode, std::string focus,
    Rcpp::NumericVector gametes)
{
    Matrix3d xd(xdS);
    Matrix3d xr(xrS);
    
    // Read A and B dimensions from xr
    const int A = xr.n_rows - 1;
    const int B = xr.n_cols - 1;

    // Build gamete vector
    if (gametes.size() == 0) {
        gametes = Rcpp::NumericVector(std::max(A, B), 0.0);
        for (unsigned int k = 0; k < gametes.size(); ++k)
            gametes[k] = dpois(k, omega);
    }
    
    // 1. Do duplication...
    if (p <= 0) {
        // Version with no additional duplication
        // a1, b1 is destination genotype and a0, b0 is source genotype
        for (int a0 = 0; a0 <= A; ++a0) {
            for (int b0 = 0; b0 <= B - a0; ++b0) {
                for (int a1 = a0; a1 <= A; ++a1) {
                    for (int b1 = b0; b1 <= B - a1; ++b1) {
                        xd(a1, b1, g) +=
                            xr(a0, b0, g - 1) * dpois(a1 - a0, u * a0) * dpois(b1 - b0, u * b0);
                    }
                }
            }
        }
    } else if (mode == "pri") {
        // Version with additional duplication (private replication)
        const int max_special_duplications = 5;
    
        // Iterate through each genotype
        for (int a0 = 0; a0 <= A; ++a0) {
            for (int b0 = 0; b0 <= B - a0; ++b0) {
                for (int a1 = a0; a1 <= A; ++a1) {
                    for (int b1 = b0; b1 <= B - a1; ++b1) {
                        if (focus == "none") {
                            // Here, the incremental B-type TE goes into [a1][b1 + bs], as normal
                            for (int bs = 0; bs < max_special_duplications && b1 + bs <= B - a1; ++bs) {
                                xd(a1, b1 + bs, g) += xr(a0, b0, g - 1) * dpois(a1 - a0, u * a0) * dpois(b1 - b0, u * b0) * dbinom(bs, b0, p);
                            }
                        } else if (focus == "old") {
                            // Here, the incremental B-type TE is cast as an A-type TE
                            for (int bs = 0; bs < max_special_duplications && a1 + bs <= A; ++bs) {
                                xd(a1 + bs, b1, g) += xr(a0, b0, g - 1) * dpois(a1 - a0, u * a0) * dpois(b1 - b0, u * b0) * dbinom(bs, b0, p);
                            }
                        } else if (focus == "new") {
                            // Here, all the old B-type TEs get cast as A-type TEs, and only the incremental B-type is cast as a B-type
                            for (int bs = 0; bs < max_special_duplications && a1 + b1 <= A; ++bs) {
                                xd(a1 + b1, bs, g) += xr(a0, b0, g - 1) * dpois(a1 - a0, u * a0) * dpois(b1 - b0, u * b0) * dbinom(bs, b0, p);
                            }
                        } else {
                            throw "Unrecognized focus.";
                        }
                    }
                }
            }
        }
    } else if (mode == "pub") {
        // Version with additional duplication (public replication)
        const int max_special_duplications = 5;
    
        // Iterate through each genotype
        for (int a0 = 0; a0 <= A; ++a0) {
            for (int b0 = 0; b0 <= B - a0; ++b0) {
                for (int a1 = a0; a1 <= A; ++a1) {
                    for (int b1 = b0; b1 <= B - a1; ++b1) {
                        int an = a1 - a0;
                        int bn = b1 - b0;
                        double pnd = dpois(an + bn, (a0 + b0) * u); // Probability of an + bn normal duplications
                        double pa = 0;
                        if (a0 + b0 > 0)
                            pa = double(a0) / (a0 + b0);
                        double pna = dbinom(an, an + bn, pa); // Probability that an of the an + bn normal duplications are type A
    
                        for (int sn = 0; sn < std::min(max_special_duplications, B - a1 - b1 + 1); ++sn) {
                            double psd = dbinom(sn, b0, p); // Probability of sn special duplications
                            for (int as = 0, bs = sn; as <= sn; ++as, --bs) {
                                if (focus == "none") {
                                    // Here, the incremental TE goes into A or B according to its type
                                    xd(a1 + as, b1 + bs, g) += xr(a0, b0, g - 1) * pnd * pna * psd * dbinom(as, sn, pa);
                                } else if (focus == "old") {
                                    // Here, the incremental TE is cast as an A-type TE
                                    xd(a1 + as + bs, b1, g) += xr(a0, b0, g - 1) * pnd * pna * psd * dbinom(as, sn, pa);
                                } else if (focus == "new") {
                                    // Here, all the old B-type TEs get cast as A-type TEs, and only an incremental B-type TE is cast as a B-type
                                    xd(a1 + as + b1, bs, g) += xr(a0, b0, g - 1) * pnd * pna * psd * dbinom(as, sn, pa);
                                } else {
                                    throw "Unrecognized focus.";
                                }
                            }
                        }
                    }
                }
            }
        }
    } else {
        throw "Unrecognized mode.";
    }
    
    // 2. Do recombination...
    // a2 + k, b2 is destination genotype and a1, b1 is source genotype
    for (unsigned int a1 = 0; a1 <= A; ++a1) {
        for (unsigned int b1 = 0; b1 <= B - a1; ++b1) {
            for (unsigned int a2 = 0; a2 <= a1; ++a2) {
                for (unsigned int b2 = 0; b2 <= b1; ++b2) {
                    for (unsigned int k = 0; a2 + k <= A && b2 <= B - a2 - k; ++k) {
                        xr(a2 + k, b2, g) +=
                            xd(a1, b1, g) * dbinom(a2, a1, 0.5) * dbinom(b2, b1, 0.5) * gametes[k];
                    }
                }
            }
        }
    }

    // Normalize xd
    double norm = 0;
    for (unsigned int a = 0; a <= A; ++a)
    {
        xd(a, 0, g) = 0; // Do not track genotypes with no mutant transposons
        for (unsigned int b = 1; b <= B; ++b)
            norm += xd(a, b, g);
    }

    for (unsigned int a = 0; a <= A; ++a)
        for (unsigned int b = 0; b <= B; ++b)
            xd(a, b, g) /= norm;
    
    // Normalize xr
    norm = 0;
    for (unsigned int a = 0; a <= A; ++a)
    {
        xr(a, 0, g) = 0; // Do not track genotypes with no mutant transposons
        for (unsigned int b = 1; b <= B; ++b)
            norm += xr(a, b, g);
    }

    for (unsigned int a = 0; a <= A; ++a)
        for (unsigned int b = 0; b <= B; ++b)
            xr(a, b, g) /= norm;
            
    Rcpp::List x = Rcpp::List::create(xdS, xrS);
    return x;
}

// Get equilibrium distribution of TE number in a one-type model
// n_max    number of host genotypes to track
// g_max    number of generations to run
// u        duplication rate
// cost     for fitness function (1-cost)^(0.5*n^2)
// n_0      starting mean copy number
// epsilon  when sum of absolute differences from one generation to the next 
//     is less than this, quit before g_max
// after_selection 
//          if false, quit after duplication; if true, quit after selection.
// [[Rcpp::export]]
Rcpp::NumericVector eq_te(
    unsigned int n_max, unsigned int g_max, 
    double u, double cost, double n_0,
    double epsilon = 0.0,
    bool after_selection = false)
{
    // Create initial distribution
    Rcpp::NumericVector x(n_max, 0.0);
    for (unsigned int i = 0; i < n_max; ++i) {
        x[i] = dpois(i, n_0);
    }
    
    // Working memory
    Rcpp::NumericVector z(n_max, 0.0);
    Rcpp::NumericVector xprev(n_max, 0.0);
    bool stationary = false;
    
    // Cycle
    for (unsigned int g = 0; true; ++g) {
        // Create xprev
        for (unsigned int i = 0; i < n_max; ++i) 
            xprev[i] = x[i];
            
        // Zero out z
        z.fill(0.0);
        
        // Duplication i -> j
        for (unsigned int i = 0; i < n_max; ++i) {
            for (unsigned int j = i; j < n_max; ++j) {
                z[j] += x[i] * dpois(j - i, u * i);
            }
        }
        
        // Quit after duplication
        if ((g == g_max || stationary) && !after_selection) {
            break;
        }

        // Selection
        double ff = 0.0;
        for (unsigned int i = 0; i < n_max; ++i) {
            z[i] *= fitness(i, cost);
            ff += z[i];
        }
        
        // Normalize
        for (unsigned int i = 0; i < n_max; ++i) {
            z[i] /= ff;
        }
        
        // Quit after selection
        if ((g == g_max || stationary) && after_selection) {
            break;
        }

        // Make gametes i -> j
        x.fill(0);
        for (unsigned int i = 0; i < n_max; ++i) {
            for (unsigned int j = 0; j <= i; ++j) {
                x[j] += z[i] * dbinom(j, i, 0.5);
            }
        }
        
        // Fuse gametes
        z.fill(0);
        for (unsigned int i = 0; i < n_max; ++i) {
            for (unsigned int j = 0; j < n_max; ++j) {
                z[i + j >= n_max ? n_max - 1 : i + j] += x[i] * x[j];
            }
        }
        
        std::copy(z.begin(), z.end(), x.begin()); // x = z;
        
        // Check for stationarity
        double absdiff = 0.0;
        for (unsigned int i = 0; i < n_max; ++i) 
            absdiff += fabs(x[i] - xprev[i]);
        if (absdiff < epsilon)
            stationary = true;
    }
    
    return z;
}


// [[Rcpp::export]]
Rcpp::NumericVector eq_mut(
    unsigned int n_max, unsigned int g_max, 
    double u, double epsilon = 0.0)
{
    // Create initial distribution
    Rcpp::NumericVector x(n_max, 1.0 / n_max);
    // for (unsigned int i = 0; i < n_max; ++i) {
    //     x[i] = ...;
    // }
    
    // Working memory
    Rcpp::NumericVector z(n_max, 0.0);
    bool stationary = false;
    
    // Cycle
    for (unsigned int g = 0; g < g_max && !stationary; ++g)
    {
        double sum = 0.0;
        for (unsigned int i = 1; i <= n_max; ++i)
        {
            z[i - 1] = 0;
            for (unsigned int j = i; j <= n_max; ++j)
            {
                for (unsigned int k = 1; k <= j; ++k)
                {
                    z[i - 1] += x[k - 1] * dpois(j - k, k * u) * dbinom(i, j, 0.5);
                }
            }
            sum += z[i - 1];
        }
        
        // Normalize
        for (unsigned int i = 0; i < n_max; ++i)
        {
            z[i] /= sum;
        }
        
        // Check for stationarity
        double absdiff = 0.0;
        for (unsigned int i = 0; i < n_max; ++i) 
            absdiff += fabs(z[i] - x[i]);
        if (absdiff < epsilon)
            stationary = true;
        
        std::copy(z.begin(), z.end(), x.begin()); // x = z;
    }
    
    return x;
}



// [[Rcpp::export]]
Rcpp::NumericVector Ns(Rcpp::NumericMatrix f, double u, double delta = 1e-4)
{
    double N0_n = 0.0,      N0_d = 0.0;
    double NPar1_n = 0.0,   NPar1_d = 0.0;
    double NOffPr1_n = 0.0, NOffPr1_d = 0.0;
    double NOffPu1_n = 0.0, NOffPu1_d = 0.0;
    
    for (unsigned int n = 0; n < 100; ++n)
    {
        for (unsigned int m = 1; m < 100 - n; ++m)
        {
            for (unsigned int i = 0; i < 25; ++i)
            {
                for (unsigned int j = 0; j < 25; ++j)
                {
                    N0_n += f(n, m) * dpois(i, n * u) * dpois(j, m * u) * (m + j) * (n + m + i + j);
                    N0_d += f(n, m) * dpois(i, n * u) * dpois(j, m * u) * (m + j);

                    NPar1_n += f(n, m) * dpois(i, n * u) * dpois(j, m * u) * (m + j) * (n + m + i + j + m * delta);
                    NPar1_d += f(n, m) * dpois(i, n * u) * dpois(j, m * u) * (m + j);

                    NOffPr1_n += f(n, m) * dpois(i, n * u) * dpois(j, m * u) * m * (n + m + i + j + 1);
                    NOffPr1_d += f(n, m) * dpois(i, n * u) * dpois(j, m * u) * m;

                    NOffPu1_n += f(n, m) * dpois(i, n * u) * dpois(j, m * u) * m * m / double(n + m) * (n + m + i + j + 1);
                    NOffPu1_d += f(n, m) * dpois(i, n * u) * dpois(j, m * u) * m * m / double(n + m);
                }
            }
        }
    }
    
    return Rcpp::NumericVector::create(
        Rcpp::Named("N0",      N0_n      / N0_d),
        Rcpp::Named("NPar1",   NPar1_n   / NPar1_d),
        Rcpp::Named("NOffPr1", NOffPr1_n / NOffPr1_d),
        Rcpp::Named("NOffPu1", NOffPu1_n / NOffPu1_d)
    );
}
