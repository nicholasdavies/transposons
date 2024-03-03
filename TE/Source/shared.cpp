// shared.cpp

#include "shared.h"
#include <ctime>
#include <vector>
#include <iostream>

Parameters P;
Randomizer R;
Recorder Rec;

double Clock()
{
    return double(std::clock()) / CLOCKS_PER_SEC;
}

static std::vector<double> clock_checkpoint_times;
static double c0;
static int N_Clocks = 10;

void StartClocking()
{
    clock_checkpoint_times.clear();
    clock_checkpoint_times.resize(N_Clocks, 0.0);
    c0 = Clock();
}

void ClockCheckpoint(int cp)
{
    double c1 = Clock();
    double elapsed = c1 - c0;
    clock_checkpoint_times[cp] += elapsed;
    c0 = c1;
}

void ShowClockInfo()
{
    if (!P.check)
    {
        clock_checkpoint_times.resize(N_Clocks, 0.0);
        std::cout << "Clocking:\n";
        for (int i = 0; i < N_Clocks; ++i)
        {
            if (clock_checkpoint_times[i] > 0)
                std::cout << "Checkpoint " << i << ": " << clock_checkpoint_times[i] << "\n";
        }
    }
}


