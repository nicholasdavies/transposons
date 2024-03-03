// Viewer - data.cpp

#include "data.h"
#include "string.hpp"
#include <numeric>
#include <fstream>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <limits>
using namespace std;

const std::string FirstLine = "# NGD2";

Datafile::Datafile(string files)
 : n_reps(0), n_cols(0), base_offset(0), rep_range(), rep(rep_range.begin())
{
    // Load all filenames specified
    vector<string> filenames = unserialize(files, ":");

    for (auto& filename : filenames)
    {
        // Open filename
        ifstream fin(filename);
        if (!fin.good())
            throw runtime_error("Could not load data file " + filename);

        // Check version
        string first_line;
        getline(fin, first_line);
        if (trim(first_line) != FirstLine)
        {
            throw runtime_error("Cannot load data files of version " + first_line);
        }

        bool datamode = false;
        unsigned int max_n_rows = 0;
        bool just_started = false;

        while (fin.good())
        {
            // get next line of data
            string line;
            getline(fin, line);
            line = trim(line);

            // line: line to skip
            if (line.empty() || line[0] == '#')
            {
                continue;
            }

            // line: data header
            else if (line == "{DATA}")
            {
                // update rep and row counts
                ++n_reps;
                n_rows.push_back(0);

                // set data mode and read labels
                datamode = true;
                getline(fin, line);
                labels = unserialize(line, "\t");
                n_cols = labels.size();

                just_started = true;
            }

            // line: plot header
            else if (line == "{PLOT}")
            {
                // set plot mode and create a new plot
                datamode = false;
                plots.push_back(PlotElements());
            }

            // line: data line
            else if (datamode)
            {
                auto row = unserialize(line, "\t");
                if (row[0] == "0" && !just_started)
                    continue;
                just_started = false;
                for (auto& d : row)
                {
                    double x;
                    istringstream iss(d);
                    iss >> x;
                    data.push_back(x);
                }
                
                ++n_rows[n_reps - 1];
                max_n_rows = max(n_rows[n_reps - 1], max_n_rows);
            }

            // line: plot line
            else
            {
                plots.back().elements.push_back(unserialize_and_decode(line));
            }
        }
    }

    used_rows = n_rows;

    SelectRep(All);
}

Datafile::Datafile(int n_rep, int n_row, int n_col, vector<string> col_labels)
 : labels(col_labels), n_reps(n_rep), n_cols(n_col), n_rows(n_reps, n_row), used_rows(n_reps, 0),
   base_offset(0), data(n_reps * n_row * n_cols, 0.0), rep_range(), rep(rep_range.begin())
{
    SelectRep(All);
}

Datafile::Datafile(int n_rep, vector<unsigned int> n_row, int n_col, vector<string> col_labels)
 : labels(col_labels), n_reps(n_rep), n_cols(n_col), n_rows(n_row), used_rows(n_reps, 0),
   base_offset(0), data(accumulate(n_row.begin(), n_row.end(), 0) * n_cols, 0.0),
   rep_range(), rep(rep_range.begin())
{
    SelectRep(All);
}

void Datafile::ShrinkToFit()
{
    unsigned int offset = 0;
    for (unsigned int r = 0; r < n_reps; ++r)
    {
        if (n_rows[r] > used_rows[r])
        {
            auto dstart = data.begin() + offset + used_rows[r] * n_cols;
            auto dend = dstart + (n_rows[r] - used_rows[r]) * n_cols;
            data.erase(dstart, dend);
            n_rows[r] = used_rows[r];
        }
        offset += n_rows[r] * n_cols;
    }
    SetBaseOffset();
}

void Datafile::Save(string filename, string comments)
{
    (void)filename;

    ofstream fout(filename);
    if (!fout.good())
        throw runtime_error("Could not open output data file " + filename);

    fout << FirstLine << "\n";

    if (!comments.empty())
    {
        fout << "\n";
        for (auto line : unserialize(comments, "\n"))
            fout << "# " << line << "\n";
        fout << "\n";
    }

    for (auto& p : plots)
    {
        fout << "{PLOT}\n";
        for (auto& line : p.elements)
        {
            fout << encode_and_serialize(line) << "\n";
        }
    }

    for (unsigned int r = 0; r < n_reps; ++r)
    {
        SelectRep(r);

        fout << "{DATA}\n";
        fout << serialize(labels, "\t") << "\n";
        for (unsigned int row = 0; row < n_rows[r]; ++row)
        {
            for (unsigned int col = 0; col < n_cols; ++col)
            {
                fout << setprecision(numeric_limits<double>::digits10 + 2) << Get(row, col);
                if (col < n_cols - 1)
                    fout << "\t";
            }
            fout << "\n";
        }
        fout << "\n";
    }
}

void Datafile::SelectRep(int r0, int r1)
{
    if (r0 == All)
    {
        r0 = 0;
        r1 = n_reps - 1;
    }

    rep_range.clear();
    if (r1 == Nil)  rep_range.push_back(r0);
    else while (r0 <= r1) rep_range.push_back(r0++);
    rep = rep_range.begin();

    SetBaseOffset();
}

void Datafile::SelectRep(string code)
{
    if (code.substr(0, 3) == "avg")
    {
        string special_range = trim(code.substr(3));
        if (special_reps.count(special_range) > 0)
        {
            // Use existing average rep
            SelectRep(special_reps[special_range]);
        }
        else
        {
            // Create a new average rep
            SelectRep(special_range);
            unsigned int new_rows = 0;
            for (StartReps(); !DoneReps(); NextRep())
            {
                new_rows = max(new_rows, NRows());
            }
            unsigned int avg_offset = data.size();
            data.insert(data.end(), new_rows * n_cols, 0.0);
            n_rows.push_back(new_rows);

            // Fill the special rep with the sum of each constituent rep
            vector<int> counts(new_rows, 0);
            for (StartReps(); !DoneReps(); NextRep())
            {
                for (unsigned int row = 0; row < NRows(); ++row)
                {
                    for (unsigned int col = 0; col < n_cols; ++col)
                        data[avg_offset + row * n_cols + col] += data[base_offset + row * n_cols + col];
                    ++counts[row];
                }
            }

            // Divide each tally by the count of entries
            for (unsigned int row = 0; row < new_rows; ++row)
                for (unsigned int col = 0; col < n_cols; ++col)
                    data[avg_offset + row * n_cols + col] /= counts[row];

            // Index the entry and select it
            special_reps[special_range] = n_rows.size() - 1;
            SelectRep(special_reps[special_range]);
        }
    }
    else
    {
        rep_range = extract_range(code, 0, n_reps - 1);
        rep = rep_range.begin();
        SetBaseOffset();
    }
}

void Datafile::StartReps()
{
    rep = rep_range.begin();
    SetBaseOffset();
}

bool Datafile::DoneReps()
{
    return rep == rep_range.end();
}

void Datafile::NextRep()
{
    ++rep;
    SetBaseOffset();
}

bool Datafile::IsSet(int row, int col) const
{
    (void) col;
    return (unsigned int)row < n_rows[*rep];
}

double Datafile::Get(int row, int col) const
{
    return data[base_offset + row * n_cols + col];
}

void Datafile::Set(int row, int col, double val)
{
    data[base_offset + row * n_cols + col] = val;
    if ((unsigned int)row + 1 > used_rows[*rep])
        used_rows[*rep] = row + 1;
}

vector<double>::iterator Datafile::RowBegin(int row)
{
    return data.begin() + base_offset + row * n_cols;
}

vector<double>::iterator Datafile::RowEnd(int row)
{
    return data.begin() + base_offset + (row + 1) * n_cols;
}

void Datafile::SetBaseOffset()
{
    if (rep != rep_range.end())
    {
        base_offset = 0;
        for (int r = 0; r < *rep; ++r)
        {
            base_offset += n_rows[r] * n_cols;
        }
    }
}
