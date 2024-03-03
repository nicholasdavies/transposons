// randomizer.h
// (C) 2013 Nicholas G Davies

#ifndef RANDOMIZER_H
#define RANDOMIZER_H

#include <random>
#include <vector>

class Randomizer
{
public:
    Randomizer();

    void Reset();

    double Uniform(double min = 0.0, double max = 1.0);
    double Normal(double mean = 0.0, double sd = 1.0);
    double Normal(double mean, double sd, double clamp);
    int Discrete(int min, int max);
    int Discrete(std::vector<double>& cumulative_weights);
    void Pick(int min, int max, int n, std::vector<int>& picks);
    bool Bernoulli(double p);
    int Binomial(int n, double p);
    int Poisson(double mean);
    int Round(double x);

    void DiehardOutput(const char* filename);

private:
    std::seed_seq seed;
    std::mt19937 engine;
};

#endif
