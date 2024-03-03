// experiment.h

#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include <vector>
#include <stdexcept>
#include "bits.h"

class Experiment
{
public:
    Experiment();

private:
    bool Manipulate();
    bool Reproduction();
    bool UpdateIndividuals();

    int g;  // current generation
    int r;  // current replication

    std::vector<Individual> population;
    std::vector<Individual> offspring;
    std::vector<double> running_fitness;
    std::vector<double> duplication_rates;

    std::vector<Element> duplicating_elements;
};

#endif