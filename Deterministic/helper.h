// helper.h - helper functions and classes.

#ifndef HELPER_H
#define HELPER_H

#include <vector>
#include "config.h"
using namespace std;

// Simulation parameters
extern Parameters P;

// Typedefs for convenience
typedef vector<double> Vector;
typedef vector<Vector> Matrix;

// Determine mode of fitness function.
void PrepareFitnessFunction();

// The fitness function
double Fitness(int na, double ua, int nb, double ub);

// Index: helper to iterate through indices of an IndexSet.
class IndexSet;
struct Index
{
    int i, j;
    bool good;

    // Constructor: used by IndexSet to initialise an Index.
    // s is the parent index set, i0 and j0 the initial indices, and if lim >= 0,
    // i and j do not iterate over indices such that i+j >= lim.
    Index(IndexSet& s, int i0, int j0, int ni, int lim)
     : i(i0), j(j0), good(true), is(s), n_internal(ni), limit(lim) { }

    // Increment to the next i, j, setting good to false if past-the-end.
    void operator++();

private:
    IndexSet& is;
    int n_internal;
    int limit;    
};

// IndexSet: a limited set of indices into a matrix, to restrict consideration to certain entries only.
class IndexSet
{
public:
    // Default constructor: make a null index set
    IndexSet() : use_indices(false), A(0), B(0) {}

    // Construct the index set to consider only elements of matrix m greater than or equal to threshold
    IndexSet(Matrix& m, int A_, int B_, double threshold)
    {
        Build(m, A_, B_, threshold);
    }

    // Reassign the index set to consider only elements of matrix m greater than or equal to threshold
    void Build(Matrix& m, int A_, int B_, double threshold)
    {
        A = A_;
        B = B_;
        use_indices = threshold > 0;

        if (use_indices)
        {
            // Construct set of indices into matrix for elements that exceed or equal threshold.
            indices.clear();
            for (int a = 0; a < min(A, (int)m.size()); ++a)
              for (int b = 0; b < min(B, (int)m[a].size()); ++b)
                if (m[a][b] >= threshold)
                {
                    indices.push_back(a);
                    indices.push_back(b);
                    indices.push_back(-1);
                }

            // indices[i + 2] holds k, such that indices[k], indices[k + 1] is the next pair of
            // matrix indices for which indices[k] + indices[k + 1] < indices[i] + indices[i + 1]
            for (int i = 0; i < (int)indices.size(); i += 3)
            {
                indices[i + 2] = std::numeric_limits<int>::max();
                for (int j = i + 3; j < (int)indices.size(); j += 3)
                {
                    if (indices[j] + indices[j + 1] < indices[i] + indices[i + 1])
                    {
                        indices[i + 2] = j;
                        break;
                    }
                }
            }
        }
    }

    // Construct an Index to iterate over relevant elements of the matrix
    Index begin()
    {
        if (use_indices)
            return Index(*this, indices[0], indices[1], 0, -1);
        else
            return Index(*this, 0, 0, 0, -1);
    }

    // Construct a subsidiary Index to iterate over relevant elements of the matrix,
    // such that the sum of the row + column, plus K's row & column, does not equal or exceed max
    Index begin_sub(Index K, int max)
    {
        int limit = max - K.i - K.j;
        if (use_indices)
            return Index(*this, indices[0], indices[1], 0, limit < 0 ? 0 : limit);
        else
            return Index(*this, 0, 0, 0, limit < 0 ? 0 : limit);
    }

private:
    friend struct Index;

    bool use_indices;
    int A, B;
    vector<int> indices;
};

// Return an upper-left "triangularish" matrix of dimension rows x cols, initialised to zeros.
// e.g. for row = col = 3, gives
//  0 0 0
//  0 0
//  0
// for rows = 3, cols = 2, gives
//  0 0
//  0 0
//  0
Matrix ZeroTriangularish(int rows, int cols);

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
void Transition(int A, int B, Matrix& f, Matrix& r, Matrix& temp, Matrix& P_A, Matrix& P_B);

// Same as above, but for upper-left triangular matrices.
void TransitionTri(int A, int B, Matrix& f, Matrix& r, Matrix& temp, Matrix& P_A, Matrix& P_B);

// Calculate the probability mass at k for Poisson(lambda).
double PoissonPMF(double lambda, int k);

// Calculate the probability mass at k for Binomial(n, p).
double BinomialPMF(int n, double p, int k);

// Return the sum of all elements in a vector.
double Sum(Vector& v);

// Return the sum of all elements in a matrix.
double Sum(Matrix& m);

// Rescale all values in the vector v so they sum to 1. Return the old sum.
double Normalize(Vector& v);

// Rescale all values in the matrix m so they sum to 1.
double Normalize(Matrix& m);

// Determine the statistical distance between two probability vectors.
double Distance(Vector& a, Vector& b);

// Determine the statistical distance between two probability matrices.
double Distance(Matrix& a, Matrix& b);

// Clock functions
double Clock();
void StartClocking();
void ClockCheckpoint(unsigned int cp);
void ShowClockInfo();

// cov(x, y) = E[ (x-E(x)) (y-E(y)) ]
// Online weighted covariance algorithm
struct Covariance
{
    double meanx, meany, wsum, C;

    Covariance()
     : meanx(0), meany(0), wsum(0), C(0)
    {
    }

    void operator()(double x, double y, double w)
    {
        if (w == 0)
            return;

        wsum += w;
        double dx = x - meanx;
        meanx += (w / wsum) * dx;
        meany += (w / wsum) * (y - meany);
        C += w * dx * (y - meany);
    }

    double Cov()
    {
        return C / wsum;
    }
};

// A sparse matrix
struct SparseMatrix
{
    struct Entry
    {
        int i, j;   // Row, column
        double x;   // Value
    };

    vector<Entry> d;    // Data

    SparseMatrix() {}

    SparseMatrix(Matrix& m, double threshold)
    {
        for (int a = 0; a < (int)m.size(); ++a)
            for (int b = 0; b < (int)m[a].size(); ++b)
                if (m[a][b] >= threshold)
                    d.push_back(Entry{a, b, m[a][b]});
    }
};

// Show a matrix, by order of magnitude
void PrintMatrix(Matrix& m);

#endif // HELPER_H