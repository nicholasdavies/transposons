// Viewer - data.h

#ifndef DATA_H
#define DATA_H

#include <string>
#include <vector>
#include <map>
#include "string.hpp"

class Datafile
{
public:
    // Elements comprising a plot.
    struct PlotElements
    {
        std::vector<std::vector<std::string>> elements;
    };

    // Load the named datafile(s).
    Datafile(std::string files);

    // Create a new datafile with the given number of rows and columns and the supplied column labels.
    Datafile(int n_rep, int n_row, int n_col, std::vector<std::string> col_labels);

    // Create a new datafile with the given number of rows (per rep) and columns and the supplied column labels.
    Datafile(int n_rep, std::vector<unsigned int> n_row, int n_col, std::vector<std::string> col_labels);

    // Remove any unused rows in the datafile
    void ShrinkToFit();

    // Save the datafile to the given file
    void Save(std::string filename, std::string comments = "", bool skip_unused = true);

    // Add a plot, with compact-notation elements (colon-separated)
    template <typename... Args> void AddPlot(Args... args);

    // Set the active rep
    enum RepCodes { All = -1, Nil = -3 };
    void SelectRep(int r0, int r1 = Nil);
    void SelectRep(std::string code);
    void StartReps();
    bool DoneReps();
    void NextRep();

    // Set/get the data
    bool IsSet(int row, int col) const;
    double Get(int row, int col) const;
    void Set(int row, int col, double val);
    std::vector<double>::iterator RowBegin(int row);
    std::vector<double>::iterator RowEnd(int row);

    // Query information about the data
    unsigned int NReps() const  { return n_reps; }
    unsigned int CurrRep() const { return *rep; }
    unsigned int NRows() const  { return n_rows[*rep]; }
    unsigned int NCols() const  { return n_cols; }

    // Labels and plots
    std::vector<std::string> labels;
    std::vector<PlotElements> plots;

private:
    // Workhorse for AddPlot
    template <typename... Args>
    void ProcessPlotElements(std::string next, Args... remaining);
    void ProcessPlotElements();

    void SetBaseOffset();

    unsigned int n_reps, n_cols;            // dimensions of data
    std::vector<unsigned int> n_rows;       // rows for a given rep
    std::vector<unsigned int> used_rows;    // furthest row that have been touched for a given rep
    unsigned int base_offset;               // for selecting "active" rep
    std::vector<double> data;               // stored flat

    std::vector<int> rep_range;             // For running through reps
    std::vector<int>::iterator rep;
    std::map<std::string, int> special_reps;
};

// Template methods
template <typename... Args>
void Datafile::AddPlot(Args... args)
{
    plots.push_back(PlotElements());
    ProcessPlotElements(args...);
}

template <typename... Args>
inline void Datafile::ProcessPlotElements(std::string next, Args... remaining)
{
    plots.back().elements.push_back(unserialize_and_decode(next));
    ProcessPlotElements(remaining...);
}

inline void Datafile::ProcessPlotElements()
{ }

#endif