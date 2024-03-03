// experiment.cpp

#include "recorder.h"
#include "experiment.h"
#include "shared.h"
#include "config.h"
#include "randomizer.h"
#include "luafunc.h"
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <fstream>
#include <csignal>
using namespace std;

static volatile bool Interrupt = false;
int Signal = 0;
void InterruptHandler(int unused)
{
    (void)unused;
    Interrupt = true;
}

Experiment::Experiment()
 : offspring(P.N), running_fitness(P.N)
 {
    // Set the interrupt handler
    signal(SIGINT, InterruptHandler);
    Interrupt = false;

    // Start the recorder
    Rec.Start();

    // Loop through all reps
    for (r = 0; r < P.reps && !Interrupt; ++r)
    {
        cout << "REP " << r + 1;

        // Time this rep
        double rc0 = Clock();

        // Reset population
        population.clear();
        population.resize(P.N);

        // Record prior to first generation
        Rec.Update(population);

        StartClocking();
        bool end = false;

        for (g = 0; g < P.generations && !end && !Interrupt; ++g)
        {
            ClockCheckpoint(0);

            // 1. EXPERIMENTAL MANIPULATION
            if (!Manipulate())
                end = true;
            ClockCheckpoint(1);

            // 2. UPDATE INDIVIDUALS
            if (!UpdateIndividuals())
                end = true;
            ClockCheckpoint(2);

            // 3. RECORDING
            bool force = (end || g == P.generations - 1);
            Rec.Update(population, force);
            ClockCheckpoint(3);

            // 4. REPRODUCTION
            if (!Reproduction())
            {
                // Ensure last generation is recorded even if reproduction fails
                if (!force && (g + 1) % P.report_freq != 0)
                    Rec.Update(population, true);
                end = true;
            }
            ClockCheckpoint(4);

            // Check for user-requested quit
            if (Interrupt)
            {
                cout << "\nCaught SIGINT.\n";
                cout << "n " << Rec.Lastn() << ", u " << Rec.Lastu() << ", S " << Rec.LastS() << ", w " << Rec.Lastw() << "\n";

                while (true)
                {
                    cout << "(g) go       (s) save & continue  (q) save & quit\n";
                    cout << "(r) next rep (n) next sweep       (x) exit without saving\n";
                    char ch;
                    cin >> ch;
                    if (ch == 'g') { Interrupt = false; break; }
                    if (ch == 's') { string filename; cout << "Save to file: "; cin >> filename; Rec.Flush(filename); Interrupt = false; break; }
                    if (ch == 'q') { Signal = 1; break; }
                    if (ch == 'r') { Interrupt = false; end = true; break; }
                    if (ch == 'n') { break; }
                    if (ch == 'x') { exit(-1); }
                }
            }

            // Display progress
            if (g % 500 == 0)
            {
                if (g > 0) cout << " n " << Rec.Lastn() << ", u " << Rec.Lastu() << ", S " << Rec.LastS() << ", w " << Rec.Lastw();
                cout << "\n" << ((g / 500) * 500) / 1000.0 << "K\t";
            }
            if (g % 10 == 0)
            {
                cout << "." << flush;
            }
        }

        // Report on rep
        cout << "\n";
        ShowClockInfo();
        Rec.EndRep();
        cout << "Rep took " << Clock() - rc0 << " seconds.\n\n";
    }

    // Report on run
    Rec.EndRun();
}

bool Experiment::Manipulate()
{
    // Conduct invasions: change a proportion of elements to mutants.
    for (int inv = 0; inv < (int)P.invasions.size(); inv += 4)
    {
        int inv_g = P.invasions[inv];
        double inv_freq = P.invasions[inv + 1];
        int inv_fam = P.invasions[inv + 2];
        double inv_u = P.invasions[inv + 3];

        if (g == inv_g)
            for (auto& i : population)
                i.Invade(inv_freq, inv_u, inv_fam);
    }

    return true;
}

bool Experiment::Reproduction()
{
    // Calculate total and running fitness; check for element clearance
    double total_fitness = 0.0;
    double total_n = 0.0;
    for (unsigned int i = 0; i < population.size(); ++i)
    {
        total_fitness += population[i].Fitness();
        running_fitness[i] = total_fitness;
        total_n += population[i].NElements();
    }

    // Check for extinction or TE clearance.
    if (total_fitness / population.size() <= P.extinction_threshold)
        return false;

    if (P.halt_on_clearance && total_n == 0)
        return false;

    if (P.halt_on_n_max && P.n_max > 0 && total_n / population.size() > P.n_max)
        return false;

    // Make offspring
    for (unsigned int o = 0; o < offspring.size(); ++o)
    {
        // Sexual reproduction from two parents selected by fitness
        if (R.Bernoulli(P.sexuality))
        {
            int p1 = R.Discrete(running_fitness), p2;
            if (R.Bernoulli(P.selfing))
                p2 = p1;
            else
                p2 = R.Discrete(running_fitness);
            offspring[o].Conceive(population[p1], population[p2]);
        }

        // Asexual reproduction from one parent selected by fitness
        else
        {
            int p = R.Discrete(running_fitness);
            offspring[o].Conceive(population[p]);
        }

        // Update summary statistics for new offspring.
        offspring[o].RefreshPostBirth();
    }

    // Insitute next generation
    population.swap(offspring);

    return true;
}

bool Experiment::UpdateIndividuals()
{
    bool coevolution = (g >= P.gen_coevolution);

    for (auto& i : population)
    {
        // Calculation of duplication rates
        double total_duplication_rate = i.DuplicationRates(duplication_rates);

        // Calculate number of each event type
        int n_duplications = R.Poisson(total_duplication_rate);
        int n_u_mutations = coevolution ? R.Binomial(i.NElements(), P.u_mut_rate) : 0;
        int n_u_dup_mutations = coevolution ? R.Binomial(n_duplications, P.u_mut_rate_on_copy) : 0;
        int n_excisions = R.Binomial(i.NElements(), P.v);
        if (P.n_max > 0)
            n_duplications = min(n_duplications, P.n_max - i.NElements() + n_excisions);

        // Gather duplicating elements
        duplicating_elements.clear();
        for (int e = 0; e < n_duplications; ++e)
            duplicating_elements.push_back(i.NthElement(R.Discrete(duplication_rates)));

        // Mutate duplicating elements (to reflect inaccurate copying)
        // This should be without replacement, but not likely a problem in practice.
        for (int e = 0; e < n_u_dup_mutations; ++e)
            duplicating_elements[R.Discrete(0, duplicating_elements.size() - 1)].MutateU();

        // Element duplication rate mutation
        // This should be without replacement, but not likely a problem in practice.
        for (int e = 0; e < n_u_mutations; ++e)
            i.MutateUNthElement(R.Discrete(0, i.NElements() - 1));

        // Element excision
        for (int e = 0; e < n_excisions; ++e)
            i.ExciseNthElement(R.Discrete(0, i.NElements() - 1));

        // Insert duplicated elements
        for (auto& e : duplicating_elements)
            i.InsertElement(e);
    }

    // Suppression mutation
    if (coevolution)
    {
        for (int s = 0; s < (int)P.s0.size(); ++s)
        {
            int n_supp_mutations = R.Poisson(population.size() * P.s_mut_rate[s] * P.ploidy);
            for (int e = 0; e < n_supp_mutations; ++e)
            {
                int target = R.Discrete(0, P.ploidy * population.size() - 1);
                population[target / P.ploidy].MutateSuppressor(s * P.ploidy + target % P.ploidy);
            }
        }
    }

    // Update summary statistics post-development.
    for (auto& i : population)
        i.RefreshPostDevelopment();

    return true;
}