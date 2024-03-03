// bits.h

#ifndef BITS_H
#define BITS_H

#include "shared.h"
#include <boost/crc.hpp>
#include <numeric>
#include <string>

// Element: a transposable element.
struct Element
{
    Element(double u_, unsigned int pos_, unsigned int id_)
     : u(u_), pos(pos_), id(id_) { }

    void MutateU()
    {
        int ui = id % P.u_mut_spread.size();
        u = Clamp(P.u_min, u + R.Normal(0, P.u_mut_spread[ui] + u * P.u_mut_spread_relative[ui]), P.u_max);
    }

    double u;           // Element duplication rate
    unsigned int pos;   // Element position within chromosome
    unsigned int id;    // Element id (for internal use) - if tracking invasion, 0 for resident and 1 for invader
};

// Chromosome: a linear chromosome.
struct Chromosome
{
    // Initialise chromosome with config-specified number of elements
    Chromosome();

    // Make this chromosome the recombinant product of c1 and c2
    void Recombine(const Chromosome& c1, const Chromosome& c2);

    std::vector<Element> els;   // Elements on this chromosome
};

// Individual: an individual with transposons and a suppressor locus
class Individual
{
public:
    // Initialise individual with config-specified defaults
    Individual();

    // Display individual's chromosomes, etc.
    void Print() const;

    // Make this individual the clonal offspring of one parent, or the sexual offspring of two parents
    void Conceive(const Individual& p);
    void Conceive(const Individual& p1, const Individual& p2);

    // Keep summary stats up-to-date
    void RefreshPostBirth();
    void RefreshPostDevelopment();

    // Return a reference to the nth element in this individual
    Element& NthElement(unsigned int n);

    // Modifiers of the individual's elements: mutate, excise, or randomly insert an element
    void MutateUNthElement(unsigned int n);             // * (see below)
    void ExciseNthElement(unsigned int n);              // *
    bool InsertElement(Element& dup);                   // *

    // Modifier of the individual's suppressor: mutate the ith suppressor gene (indexed by suppressor_n x ploid)
    void MutateSuppressor(unsigned int i);              // *

    // Change a proportion of the individual's elements to mutants of a specified duplication rate and id/family.
    // If u or id are negative, do not change that property.
    void Invade(double proportion, double u, int id); // *

    // Access individual's statistics
    int NElements() const           { return n_te; }
    int NElements(int t) const      { return n_type.empty() ? n_te : n_type[t]; }
    double S() const                { return std::accumulate(sbv.begin(), sbv.end(), 0.0); }
    double S(int i) const           { return sbv[i]; }
    double Activity() const         { return activity; }
    double Activity(int f) const    { return a_fam[f]; }
    double Fitness() const          { return fitness; }

    // Fill rates with a running sum of this individual's element's
    // duplication rates, accounting for suppression, parasitism, etc; return the total.
    double DuplicationRates(std::vector<double>& rates) const;

    // Supply the checksum object with this individual's basic properties
    // such that the difference between two separate runs can be established.
    void ProcessChecksum(boost::crc_32_type& cs);

private:
    friend class Recorder;
    friend class Experiment;

    // Basic properties
    std::vector<Chromosome> chr;    // Individual's chromosomes
    std::vector<double> sg;         // Individual's suppressor locus genes

    // Statistics updated in RefreshPostBirth, and kept accurate by the starred (*) functions above
    int n_te;                       // Total number of elements (= occupied spaces in genome)
    std::vector<double> sbv;        // Suppression breeding values

    // Statistics updated in RefreshPostDevelopment, invalidated by the starred functions above
    double fitness;                 // Individual's fitness
    double activity;                // Sum of duplication rate over all elements
    std::vector<int> n_type;        // Count of elements by type: zygosity (0 hetero, 1 homo) x class of site occupied x family
    std::vector<double> a_fam;      // Element activity by family (not type!) 
};

// Helpers for tracking transposon types
int NTETypes();
int TEType(int family, int site_class, bool homozygous);
std::string TETypeName(int t);

#endif