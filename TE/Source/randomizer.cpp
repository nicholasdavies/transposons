// randomizer.cpp
// (C) 2013 Nicholas G Davies

#include "randomizer.h"
#include <cmath>
#include <fstream>

Randomizer::Randomizer()
 : seed({ 120368978, 37761590, 364135833, 399444367, 298076336 })
{
    Reset();
}

void Randomizer::Reset()
{
    engine.seed(seed);
}

double Randomizer::Uniform(double min, double max)
{
    std::uniform_real_distribution<double> dist(min, max);
    return dist(engine);
}

double Randomizer::Normal(double mean, double sd)
{
    std::normal_distribution<double> dist(mean, sd);
    return dist(engine);
}

double Randomizer::Normal(double mean, double sd, double clamp)
{
    double n;
    do n = Normal(mean, sd); while (std::fabs(n - mean) > clamp);
    return n;
}

int Randomizer::Discrete(int min, int max)
{
    std::uniform_int_distribution<int> dist(min, max);
    return dist(engine);
}

int Randomizer::Discrete(std::vector<double>& cumulative_weights)
{
    return std::distance(cumulative_weights.begin(),
                std::lower_bound(cumulative_weights.begin(),
                    cumulative_weights.end(),
                    this->Uniform(0.0, cumulative_weights.back())));
}

void Randomizer::Pick(int min, int max, int n, std::vector<int>& picks)
{
    picks.clear();

    for (int j = max - n + 1; j <= max; ++j)
    {
        int t = Discrete(min, j);
        auto ins = lower_bound(picks.begin(), picks.end(), t);
        if (ins == picks.end() || *ins != t)
            picks.insert(ins, t);
        else
            picks.push_back(j);
    }
}

bool Randomizer::Bernoulli(double p)
{
    if (p <= 0) return false;
    if (p >= 1) return true;
    std::bernoulli_distribution dist(p);
    return dist(engine);
}

int Randomizer::Binomial(int n, double p)
{
    if (p <= 0) return 0;
    std::binomial_distribution<int> dist(n, p);
    return dist(engine);
}

int Randomizer::Poisson(double mean)
{
    if (mean <= 0) return 0;
    std::poisson_distribution<int> dist(mean);
    return dist(engine);
}

int Randomizer::Round(double x)
{
    int sign = x < 0 ? -1 : 1;
    double intpart, fracpart;
    fracpart = std::modf(std::fabs(x), &intpart);
    return sign * (intpart + Bernoulli(fracpart));
}

void Randomizer::DiehardOutput(const char* filename)
{
    std::ofstream out(filename);
    for (int i = 0; i < 12000000;)
    {
        auto k = engine();
        out.write(reinterpret_cast<char*>(&k), sizeof(k));
        i += sizeof(k);
    }
}
