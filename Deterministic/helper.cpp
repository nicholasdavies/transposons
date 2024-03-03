// helper.cpp - helper functions and classes.

#include "helper.h"
#include "string.hpp"
#include <numeric>
#include <ctime>
#include <iostream>
#include <boost/math/distributions/poisson.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <stdexcept>

// Simulation parameters
Parameters P;

// Determine mode of fitness function.
static int FitnessMode = 0;
void PrepareFitnessFunction()
{
    string args = P.w.CompactForm();
    args = replace(args.substr(0, args.find("->")), " ", "");
    
    // can take: n; na,nb; n,u; na,ua,nb,ub
    if (args == "n")
        FitnessMode = 1;
    else if (args == "na,nb")
        FitnessMode = 2;
    else if (args == "n,u")
        FitnessMode = 3;
    else if (args == "na,ua,nb,ub")
        FitnessMode = 4;
    else
        throw runtime_error("Unrecognised arguments to fitness function (" + args + "). Must be one of n; na,nb; n,u; na,ua,nb,ub.");
}

// The fitness function
double Fitness(int na, double ua, int nb, double ub)
{
    double w = 0;

    switch (FitnessMode)
    {
    case 1:
        w = P.w(na + nb);
        break;
    case 2:
        w = P.w(na, nb);
        break;
    case 3:
        w = P.w(na + nb, (na * ua + nb * ub) / (na + nb > 0 ? na + nb : 1));
        break;
    case 4:
        w = P.w(na, ua, nb, ub);
        break;
    default:
        w = 0;
        break;
    }

    return max(0.0, w);
}

// Increment to the next i, j, setting good to false if past-the-end.
void Index::operator++()
{
    if (is.use_indices)
    {
        n_internal += 3;

        while (true)
        {
            if (n_internal >= (int)is.indices.size())
            {
                good = false;
                return;
            }
            i = is.indices[n_internal];
            j = is.indices[n_internal + 1];

            if (limit >= 0 && i + j >= limit)
                n_internal = is.indices[n_internal - 1];
            else
                break;
        }
    }
    else
    {
        int max = limit < 0 ? is.A : limit;
        if (++j >= (is.B == 1 ? 1 : max - i))
        {
            j = 0;
            if (++i >= max)
                good = false;
        }
    }
}

// Return an upper-left "triangularish" matrix of dimension rows x cols, initialised to zeros.
// e.g. for row = col = 3, gives
//  0 0 0
//  0 0
//  0
// for rows = 3, cols = 2, gives
//  0 0
//  0 0
//  0
Matrix ZeroTriangularish(int rows, int cols)
{
    Matrix m;
    for (int r = 0; r < rows; ++r)
        m.push_back(Vector(min(rows - r, cols), 0.0));

    return m;
}

// For an AxB matrix f, an AxA matrix P_A and a BxB matrix P_B, performs
//
//    r [I,J]    =     SUM     { f[i,j] * P_A[i,I] * P_B[j,J] }
//  0≤I<A,0≤J<B    0≤i<A,0≤j<B
//
// using the temporary AxB matrix temp. r can be the same matrix as f.
//
// (In other words, if f[i,j] is the frequency of type-ij individuals, P_A[i,I] is the probability
// of a transition from type i,j to type I,j and P_B[j,J] is the probability of a transition from
// type i,j to type i,J, this populates all r[I,J] with frequencies of type-IJ individuals after
// transitioning.)
//
// This calculation is O(AABB), but can be made O(AAB + ABB) by breaking down into two steps:
//
//   temp[I,J]  =  SUM  { f[i,J] * P_A[i,I] }
//  0≤I<A,0≤J<B   0≤i≤A
//    r [I,J]   =  SUM  { temp[I,j] * P_B[j,J] }
//  0≤I<A,0≤J<B   0≤j≤B
//
// so after the first step, temp[I,J] contains the frequencies following all i->I transitions and
// after the second step, r[I,J] contains the frequencies following all i,j->I,J transitions.
//
void Transition(int A, int B, Matrix& f, Matrix& r, Matrix& temp, Matrix& P_A, Matrix& P_B)
{
    // First do dimension-1 transition
    for (int I = 0; I < A; ++I)
    {
        for (int J = 0; J < B; ++J)
        {
            temp[I][J] = 0;
            for (int i = 0; i < A; ++i)
            {
                temp[I][J] += f[i][J] * P_A[i][I];
            }
        }
    }

    // Now do dimension-2 transition
    for (int I = 0; I < A; ++I)
    {
        for (int J = 0; J < B; ++J)
        {
            r[I][J] = 0;
            for (int j = 0; j < B; ++j)
            {
                r[I][J] += temp[I][j] * P_B[j][J];
            }
        }
    }
}

// Same as above, but for upper-left triangular matrices.
void TransitionTri(int A, int B, Matrix& f, Matrix& r, Matrix& temp, Matrix& P_A, Matrix& P_B)
{
    if (B == 1)
    {
        // First assign to temp
        for (int I = 0; I < A; ++I)
        {
            temp[I][0] = 0;
            for (int i = 0; i < A; ++i)
                temp[I][0] += f[i][0] * P_A[i][I];
        }

        // Now assign to r - this is needed because r and f might be pointing to same memory
        for (int I = 0; I < A; ++I)
            r[I][0] = temp[I][0];
    }

    else if (A == B)
    {
        // First do dimension-1 transition
        for (int I = 0; I < A; ++I)
        {
            for (int J = 0; J < B - I; ++J)
            {
                temp[I][J] = 0;
                for (int i = 0; i < A - J; ++i)
                {
                    temp[I][J] += f[i][J] * P_A[i][I];
                }
            }
        }

        // Now do dimension-2 transition
        for (int I = 0; I < A; ++I)
        {
            for (int J = 0; J < B - I; ++J)
            {
                r[I][J] = 0;
                for (int j = 0; j < B - I; ++j)
                {
                    r[I][J] += temp[I][j] * P_B[j][J];
                }
            }
        }
    }

    else if (B < A && (P.mode == "analysis" || P.mode == "analysis2")) // experimental . . .
    {
        // First do dimension-1 transition
        for (int I = 0; I < A; ++I)
        {
            for (int J = 0; J < min(B, A - I); ++J)
            {
                temp[I][J] = 0;
                for (int i = 0; i < A - J; ++i)
                {
                    temp[I][J] += f[i][J] * P_A[i][I];
                }
            }
        }

        // Now do dimension-2 transition
        for (int I = 0; I < A; ++I)
        {
            for (int J = 0; J < min(B, A - I); ++J)
            {
                r[I][J] = 0;
                for (int j = 0; j < min(B, A - I); ++j)
                {
                    r[I][J] += temp[I][j] * P_B[j][J];
                }
            }
        }
    }

    else
        throw std::runtime_error("TransitionTri can only handle flat or square triangular matrices.\n");
}

// Calculate the probability mass at k for Poisson(lambda).
double PoissonPMF(double lambda, int k)
{
    if (lambda <= 0.0)  return k == 0;
    if (k < 0)  return 0.0;
    return boost::math::pdf(boost::math::poisson(lambda), k);
}

// Calculate the probability mass at k for Binomial(n, p).
double BinomialPMF(int n, double p, int k)
{
    if (k > n || k < 0) return 0.0;
    return boost::math::pdf(boost::math::binomial(n, p), k);
}

// Return the sum of all elements in a vector.
double Sum(Vector& v)
{
    return accumulate(v.begin(), v.end(), 0.0);
}

// Return the sum of all elements in a matrix.
double Sum(Matrix& m)
{
    return accumulate(m.begin(), m.end(), 0.0, [](double r, Vector& v) { return r + Sum(v); });
}

// Rescale all values in the vector v so they sum to 1. Return the old sum.
double Normalize(Vector& v)
{
    double s = Sum(v);
    for (auto& x : v)
        x /= s;
    return s;
}

// Rescale all values in the matrix m so they sum to 1.
double Normalize(Matrix& m)
{
    double s = Sum(m);
    for (auto& v : m)
        for (auto& x : v)
            x /= s;
    return s;
}

// Determine the statistical distance between two probability vectors.
double Distance(Vector& a, Vector& b)
{
    if (a.size() != b.size())
        throw std::runtime_error("Mismatched vector size in Distance().");
    double d = 0.0;
    for (unsigned int i = 0; i < a.size(); ++i)
        d += std::fabs(a[i] - b[i]);
    return 0.5 * d;
}

// Determine the statistical distance between two probability matrices.
double Distance(Matrix& a, Matrix& b)
{
    if (a.size() != b.size())
        throw std::runtime_error("Mismatched matrix size (rows) in Distance().");
    double d = 0.0;
    for (unsigned int i = 0; i < a.size(); ++i)
    {
        if (a[i].size() != b[i].size())
        {
            cout << "Row " << i << ": " << a[i].size() << " " << b[i].size() << "\n";
            throw std::runtime_error("Mismatched matrix size (cols) in Distance().");
        }
        for (unsigned int j = 0; j < a[i].size(); ++j)
        {
            d += std::fabs(a[i][j] - b[i][j]);
        }
    }
    return 0.5 * d;
}

// Clock functions
static std::vector<double> ClockTimes;
static double C0;

double Clock()
{
    return double(clock()) / CLOCKS_PER_SEC;
}

void StartClocking()
{
    ClockTimes.clear();
    C0 = Clock();
}

void ClockCheckpoint(unsigned int cp)
{
    double C1 = Clock();
    if (cp >= ClockTimes.size()) ClockTimes.resize(cp + 1, 0.0);
    ClockTimes[cp] += C1 - C0;
    C0 = C1;
}

void ShowClockInfo()
{
    for (unsigned int i = 0; i < ClockTimes.size(); ++i)
        cout << "Checkpoint " << i << ": " << ClockTimes[i] << "\n";
}


// Show a matrix, by order of magnitude
void PrintMatrix(Matrix& m)
{
    double max_freq = 0;
    for (int i = 0; i < (int)m.size(); ++i)
        for (int j = 0; j < (int)m[i].size(); ++j)
            if (m[i][j] > max_freq)
                max_freq = m[i][j];

    if (m[0].size() == 1)
    {
        for (int a = 0; a < (int)m.size(); ++a)
        {
            if (m[a][0] == max_freq)
                cout << "*";
            else
                cout << min(9, max(0, 9 + (int)log10(m[a][0] / max_freq)));
        }
        cout << "\n";
    }
    else
    {
        cout << "  B →\n";
        for (int a = 0; a < (int)m.size(); ++a)
        {
            cout << (a == 0 ? "A" : (a == 1 ? "↓" : " ")) << ' ';
            for (int b = 0; b < (int)m[a].size(); ++b)
            {
                if (m[a][b] == max_freq)
                    cout << "*";
                else
                    cout << min(9, max(0, 9 + (int)log10(m[a][b] / max_freq)));
            }
            cout << "\n";
        }
    }
}