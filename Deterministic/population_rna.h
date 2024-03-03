// population_rna.h

#ifndef POPULATION_RNA_H
#define POPULATION_RNA_H

#include "helper.h"
#include "data.h"

// A population of individuals with privately-replicating transposons
class Population
{
public:
    Population();
    virtual ~Population();

    // Return the number of columns in row r of the frequency matrix.
    virtual int Minor(int r) const;

    // Initialize population
    virtual void Init();

    // Do lookup tables for fitness
    virtual void ConstructFitness();

    // Populate lookup tables
    virtual void MakeLookupTables();

    // Swap lookup tables with buffered version (for analysis2)
    virtual void SwapLookupTables();

    // Calculate statistics
    virtual void CalculateStatistics();

    // Get duplication rate of wild-type and mutant transposons
    virtual double DuplicationRateWildType(int a, int b) const;
    virtual double DuplicationRateMutant(int a, int b) const;

    // Calculate average fitness
    virtual double AvgW() const;

    // Calculate proportion of individuals with in-bound genotypes who are inviable
    virtual double Inviables();

    // Update genotypes through duplication and loss of transposons.
    virtual void Update();

    // Do viability selection, multiplying the frequency of each type by its relative fitness.
    // This is equivalent to killing off a fraction of each type, where the survival probability is proportional to fitness.
    virtual void Selection();

    // Reproduction
    virtual void Reproduce();

    // "Special" update step for analysis mode
    virtual void SpecialUpdate();

    // Print some basic statistics about the population
    virtual void PrintStatistics();

    // Graphically show the distribution of genotypes in the population
    virtual void ShowGenotypeDistribution();

    // Report progress through PrintStatistics(), optionally at a certain frequency
    virtual void Report(int g, string stage = "", bool always = false);

    // Record a generation
    virtual void RecordGeneration(Datafile& df, int g);

    // Record a generation, analysis2 mode
    virtual void Analysis2RecordGeneration(Datafile& df, int row, double g);

    // Initialise and run simulation with one transposon type until a certain condition is met.
    // Parameters:
    //  n_generations   maximum number of generations for which to run the simulation.
    //  epsilon         simulation terminates early if distance between two consecutive population
    //                  frequency matrices is less than epsilon. can pass negative number to override.
    //  verbose         if true, reports on progress every P.report_freq generations.
    //  df              if non-null, records statistics to df every generation.
    //  record_start    which datafile row to start recording in.
    virtual void SingleRun(int n_generations, double epsilon, bool verbose, Datafile* df = 0, int record_start = 0);

    // Introduce mutants to a wild-type-only population.
    virtual void IntroduceMutants(double mutant_freq, int n_B = -1);

    // Consolidate all transposons to wild-type.
    virtual void Consolidate();

    // Run simulation
    virtual void Run();

protected:
    int A, B;           // Number of genotypes of either transposon type
    Matrix f;           // Genotype frequencies
    Matrix w;           // Genotype fitnesses
    double N_A, N_B;    // Average copy number for either transposon type
    double V_A, V_B;    // Variance in copy number for either transposon type
    double Av_W, V_W;   // Average host fitness and variance in host fitness
    double A05, A95;    // 90% range in copy number of transposon A
    double B05, B95;    // 90% range in copy number of transposon B
    double W05, W95;    // 90% range in host fitness
    double max_freq;    // Frequency of genotype with maximum share of population
    double still, slow; // Proportion of individuals who do not survive birth / life
    double f_other;     // Proportion of individuals with genotypes beyond n_max (assumed w=0)
    double saved_sWU;   // For analysis2 mode
    double saved_b;     // For analysis2 mode

    // Lookup tables are indexed [g][e], where g is the genotype (transposon number), and e is the type
    // of event (number of duplications, number of excisions, transition genotype, or gamete copy number).
    // The value looked up is the probability of event e happening for an individual of genotype g.
    Matrix DuplicationP_A;
    Matrix DuplicationP_B;
    Matrix ExcisionP_A;
    Matrix ExcisionP_B;
    Matrix TransitionP_A;
    Matrix TransitionP_B;
    Matrix GameteP_A;
    Matrix GameteP_B;

    // Storage for swapped lookup tables
    vector<Matrix> Swap;

    // Genotypes sorted by fitness (low to high)
    struct Lookup { int a, b; double w; };
    vector<Lookup> W_lookup;
    vector<double> W_list;

    // Temporary storage
    Matrix f_temp;      // For the Transition function
    Matrix gametes;     // Storage of gametes
    Matrix offspring;   // For alternating between parental / offspring generations
    IndexSet gamete_is; // Gamete IndexSet storage
};

#endif // POPULATION_RNA_H