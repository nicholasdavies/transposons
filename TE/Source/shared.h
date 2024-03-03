// shared.h

#ifndef SHARED_H
#define SHARED_H

#include "config.h"
#include "randomizer.h"
#include "recorder.h"

extern Parameters P;
extern Randomizer R;
extern Recorder Rec;
extern int Signal;

inline double Clamp(double min, double val, double max)
{
    return val < min ? min : (val > max ? max : val);
}

double Clock();

void StartClocking();
void ClockCheckpoint(int cp);
void ShowClockInfo();

#endif