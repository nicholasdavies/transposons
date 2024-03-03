// recorder.h

#ifndef RECORDER_H
#define RECORDER_H

#include "data.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <boost/crc.hpp>

class Individual;

struct Covariance
{
    double ex, ey, exy;
    double sx, sy;
    int n;
    Covariance();
    void Reset();
    void Push(double x, double y);
    double Cov() const;
};

// Recorder - records details of a population.
class Recorder
{
public:
    Recorder();

    void Start();

    void Update(std::vector<Individual>& population, bool force = false);

    void Flush(std::string filename);

    void EndRep();
    void EndRun();

    double Lastn() const { return avg_n; }
    double Lastu() const { return avg_u; }
    double LastS() const { return avg_S; }
    double Lastw() const { return avg_w; }

private:
    Datafile df;

    boost::crc_32_type cs;

    int g, row;
    int rep;

    double avg_n, avg_u, avg_S, avg_w;
    std::vector<double> duplication_rates;

    void Reset();
};

#endif