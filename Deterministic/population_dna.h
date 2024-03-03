// population_dna.h

#ifndef POPULATION_DNA_H
#define POPULATION_DNA_H

#include "helper.h"
#include "data.h"
#include "population_rna.h"

// A population of individuals with publicly-replicating transposons
class PopulationDNA : public Population
{
public:
    PopulationDNA();

    // Populate lookup tables
    void MakeLookupTables();

    // Swap lookup tables with buffered version (for analysis2)
    void SwapLookupTables();

    // Get duplication rate of wild-type and mutant transposons
    double DuplicationRateWildType(int a, int b) const;
    double DuplicationRateMutant(int a, int b) const;

    // Update genotypes through duplication and loss of transposons.
    void Update();

    // "Special" update step for analysis mode
    void SpecialUpdate();

private:
    // Matrix of transition probability matrices
    vector<vector<SparseMatrix>> TransitionP;

    // Swap version of above
    vector<vector<SparseMatrix>> TransitionP_Swap;
};


#endif