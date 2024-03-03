// deterministic.cpp
// A deterministic population genetics model of transposon-host dynamics.
// Throughout, we denote "wild-type" transposons by "A" and "mutant" transposons by "B".

#include <iostream>
#include "config.h"
#include "data.h"
#include "helper.h"
#include "population_rna.h"
#include "population_dna.h"
using namespace std;

int main(int argc, char* argv[])
{
    double clock_start = Clock();

    // Read parameters
    P.Read(argc, argv);

    // Run through parameter sweeps
    for (int s = 0; s < P.NSweeps(); ++s)
    {
        // Determine type of fitness function
        PrepareFitnessFunction();

        // Output sweep name and parameters
        cout << P.SweepName() << ":\n";
        P.Write(cout);
        cout << "\n";

        // Run simulation
        if (P.type == "RNA")
        {
            Population pop;
            pop.Run();
        }
        else if (P.type == "DNA")
        {
            PopulationDNA pop;
            pop.Run();
        }
        else
        {
            cout << "Unrecognised type " << P.type << ". (Must be RNA or DNA.)\n";
        }

        P.NextSweep();
    }

    cout << "Simulation took " << Clock() - clock_start << " seconds.\n";
}