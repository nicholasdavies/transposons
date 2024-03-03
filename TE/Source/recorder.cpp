// recorder.cpp

#include "recorder.h"
#include "bits.h"
#include <sstream>
#include <iostream>
#include <vector>
#include <array>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics/weighted_mean.hpp>
#include <boost/accumulators/statistics/weighted_variance.hpp>
#include <boost/accumulators/statistics/extended_p_square_quantile.hpp>
using namespace boost::accumulators;
using namespace std;

struct Accumulator
{
public:
    static array<double, 2> probs;

    Accumulator()
     : accum(tag::weighted_extended_p_square::probabilities = probs)
    { }

    void Reset()                    { accum = accum_t(tag::weighted_extended_p_square::probabilities = probs); }
    void Push(double x)             { accum(x, weight = 1); }
    void Push(double x, double w)   { accum(x, weight = w); }

    double Mean() const             { return weighted_mean(accum); }
    double Variance() const         { return lazy_weighted_variance(accum); }
    double Interval05() const       { return extended_p_square(accum)[0]; }
    double Interval95() const       { return extended_p_square(accum)[1]; }

private:
    typedef accumulator_set<double, stats<tag::weighted_mean,
        tag::lazy_weighted_variance, tag::weighted_extended_p_square>, double> accum_t;
    accum_t accum;
};

array<double, 2> Accumulator::probs = { 0.05, 0.95 };

Covariance::Covariance()
 : sx(0), sy(0)
{
    Reset();
}

void Covariance::Reset()
{
    ex = ey = exy = 0;
    n = 0;
}

void Covariance::Push(double x, double y)
{
    // Initialise shifts
    if (n == 0)
    {
        sx = x;
        sy = y;
    }

    ex += x - sx;
    ey += y - sy;
    exy += (x - sx) * (y - sy);
    ++n;
}

double Covariance::Cov() const
{
    return n > 0 ? (exy - ex * ey / n) / n : 0;
}

Recorder::Recorder()
 : df(0, 0, 0, {}),
   rep(0), 
   avg_n(0), avg_u(0), avg_S(0), avg_w(0)
{
}

void Recorder::Start()
{
    // Set columns in datafile
    vector<string> columns = { "g", "n", "Vn", "n05", "n95", "u", "Vu", "u05", "u95",
                                    "U", "VU", "U05", "U95", "a", "Va", "a05", "a95",
                                    "S", "VS", "S05", "S95", "w", "Vw", "w05", "w95",
                                    "x", "Vx", "x05", "x95", "d", "G", "BuW", "Bun" };

    // Copy number of each transposon type
    if (NTETypes() > 1) for (int t = 0; t < NTETypes(); ++t)
    {
        string type = TETypeName(t);
        columns.insert(columns.end(), { type, "V" + type, type + "05", type + "95" });
    }

    // Rate, suppressed rate, and activity for each transposon family
    if (P.families.size() > 1) for (int f = 0; f < (int)P.families.size(); ++f)
    {
        string fid = to_string(f);
        columns.insert(columns.end(), { "u" + fid, "Vu" + fid, "u" + fid + "05", "u" + fid + "95" });
        columns.insert(columns.end(), { "U" + fid, "VU" + fid, "U" + fid + "05", "U" + fid + "95" });
        columns.insert(columns.end(), { "a" + fid, "Va" + fid, "a" + fid + "05", "a" + fid + "95" });
    }

    // Suppression for each suppressor type
    if (P.s0.size() > 1) for (int s = 0; s < (int)P.s0.size(); ++s)
    {
        string sid = to_string(s);
        columns.insert(columns.end(), { "S" + sid, "VS" + sid, "S" + sid + "05", "S" + sid + "95" });
    }

    df = Datafile(P.reps, (P.generations / P.report_freq) + 2, columns.size(), columns);
    rep = 0;
    cs.reset();

    df.AddPlot("plot:name=Overview", "axes:0:y_min=0,x=0,y=250,w=480,h=250,y_label=Transposon abundance", "ribbon2:*:g:n05:n95:fill=$jet,x_axis=0", "line:*:g:n:stroke=#000000,x_axis=0", "axes:2:y_min=0,x=0,y=0,w=480,h=250,y_label=Duplication rate", "ribbon2:*:g:u05:u95:fill=$greens,x_axis=2", "line:*:g:u:stroke=#000000,x_axis=2", "axes:4:y_min=0,x=480,y=250,w=480,h=250,y_label=Suppression", "ribbon2:*:g:S05:S95:fill=$reds,x_axis=4", "line:*:g:S:stroke=#000000,x_axis=4", "axes:6:y_min=0,y_max=1,x=480,y=0,w=480,h=250,y_label=Fitness", "ribbon2:*:g:w05:w95:fill=$blues,x_axis=6", "line:*:g:w:stroke=#000000,x_axis=6");
    df.AddPlot("plot:name=Abundance", "axes:0:y_label=Transposon abundance,y_min=0", "ribbon2:*:g:n05:n95:fill=$jet", "line:*:g:n:stroke=#000000");
    df.AddPlot("plot:name=Duplication rate", "axes:0:y_label=Duplication rate,y_min=0", "ribbon2:*:g:u05:u95:fill=$greens", "line:*:g:u:stroke=#000000");
    df.AddPlot("plot:name=Adjusted duplication rate", "axes:0:y_label=Adjusted duplication rate,y_min=0", "ribbon2:*:g:U05:U95:fill=$greens", "line:*:g:U:stroke=#000000");
    df.AddPlot("plot:name=Suppression", "axes:0:y_label=Suppression,y_min=0", "ribbon2:*:g:S05:S95:fill=$reds", "line:*:g:S:stroke=#000000");
    df.AddPlot("plot:name=Activity", "axes:0:y_min=0,y_label=Activity", "ribbon2:*:g:a05:a95:fill=$greys", "line:*:g:a:stroke=#000000");
    df.AddPlot("plot:name=Fitness", "axes:0:y_min=0,y_max=1,x=0,y=250,w=960,h=250,y_label=Fitness / Dead", "ribbon2:*:g:w05:w95:fill=$blues,x_axis=0", "line:*:g:w:stroke=#000000,x_axis=0", "line:*:g:d:stroke=#ff0000,x_axis=0", "axes:2:y_min=0,y_max=1,x=0,y=0,w=960,h=250,y_label=Host fitness", "ribbon2:*:g:x05:x95:fill=$blues,x_axis=2", "line:*:g:x:stroke=#000000,x_axis=2");
    df.AddPlot("plot:name=Stats", "axes:0:y_min=0,y_max=1,x=0,y=250,w=480,h=250,y_label=Homozygosity", "line:*:g:G:stroke=#666666,x_axis=0", "axes:2:y_min=0,x=0,y=0,w=480,h=250,y_label=Var(n)/n", "line:*:g:Vn/n:stroke=#cc6600,x_axis=2", "axes:4:y_min=0,x=480,y=250,w=480,h=250,y_label=-BuW", "line:*:g:-BuW:stroke=#008888,x_axis=4", "axes:6:y_min=0,x=480,y=0,w=480,h=250,y_label=Bun", "line:*:g:Bun:stroke=#880000,x_axis=6");

    Reset();
}

void Recorder::Update(vector<Individual>& population, bool force)
{
    // CHECKSUM
    if (P.do_checksum)
    {
        for (auto& i : population)
            i.ProcessChecksum(cs);
    }

    if (g % P.report_freq == 0 || force)
    {
        Accumulator _nn, _uu, _UU, _aa, _SS, _w, _x, _d, _G;
        //Accumulator _ax, _SSx;
        Covariance _CovuW, _Covun;
        vector<Accumulator> _n(NTETypes());
        vector<Accumulator> _u(P.families.size());
        vector<Accumulator> _U(P.families.size());
        vector<Accumulator> _a(P.families.size());
        vector<Accumulator> _S(P.s0.size());

        for (auto& i : population)
        {
            // Per-individual
            double w = i.Fitness();
            _nn.Push(i.NElements());
            _aa.Push(i.Activity());
            _SS.Push(i.S());
            _w.Push(w);
            _x.Push(w, i.NElements());
            _d.Push(w > 0 ? 0 : 1);
            //_ax.Push(i.Activity(), i.NElements());
            //_SSx.Push(i.S(), i.NElements());

            // Per TE family
            for (int f = 0; f < (int)P.families.size(); ++f)
                _a[f].Push(i.Activity(f));

            // Per-TE type & homozygosity
            double het = 0, hom = 0;
            for (int t = 0; t < NTETypes(); ++t)
            {
                _n[t].Push(i.NElements(t));
                if (P.zygosity)
                {
                    if (t % 2 == 0)
                        het += i.NElements(t);
                    else
                        hom += i.NElements(t) / 2;
                }
            }
            if (het + hom > 0)  _G.Push(hom / (het + hom));

            // Per suppressor type
            for (int s = 0; s < (int)P.s0.size(); ++s)
                _S[s].Push(i.S(s));

            // Per-TE
            i.DuplicationRates(duplication_rates);
            auto e2 = duplication_rates.begin();
            double U = 0, prev_U_sum = 0;
            for (auto& c : i.chr)
            {
                for (auto& e : c.els)
                {
                    _uu.Push(e.u);
                    _u[e.id].Push(e.u);
                    _CovuW.Push(e.u, w);
                    _Covun.Push(e.u, i.NElements());

                    U = *e2 - prev_U_sum;
                    _UU.Push(U);
                    _U[e.id].Push(U);
                    prev_U_sum = *e2;
                    ++e2;
                }
            }   
        }

        df.Set(row,  0, g);
        df.Set(row,  1, _nn.Mean()); df.Set(row,  2, _nn.Variance()); df.Set(row,  3, _nn.Interval05()); df.Set(row,  4, _nn.Interval95());
        df.Set(row,  5, _uu.Mean()); df.Set(row,  6, _uu.Variance()); df.Set(row,  7, _uu.Interval05()); df.Set(row,  8, _uu.Interval95());
        df.Set(row,  9, _UU.Mean()); df.Set(row, 10, _UU.Variance()); df.Set(row, 11, _UU.Interval05()); df.Set(row, 12, _UU.Interval95());
        df.Set(row, 13, _aa.Mean()); df.Set(row, 14, _aa.Variance()); df.Set(row, 15, _aa.Interval05()); df.Set(row, 16, _aa.Interval95());
        df.Set(row, 17, _SS.Mean()); df.Set(row, 18, _SS.Variance()); df.Set(row, 19, _SS.Interval05()); df.Set(row, 20, _SS.Interval95());
        df.Set(row, 21,  _w.Mean()); df.Set(row, 22,  _w.Variance()); df.Set(row, 23,  _w.Interval05()); df.Set(row, 24,  _w.Interval95());
        df.Set(row, 25,  _x.Mean() / _w.Mean());
        df.Set(row, 26,  _x.Variance() / (_w.Mean() * _w.Mean()));
        df.Set(row, 27,  _x.Interval05() / _w.Mean());
        df.Set(row, 28,  _x.Interval95() / _w.Mean());
        df.Set(row, 29, _d.Mean());
        df.Set(row, 30, _G.Mean());
        df.Set(row, 31, (_CovuW.Cov() / _w.Mean()) / _uu.Variance());
        df.Set(row, 32, _Covun.Cov() / _uu.Variance());

        int c = 33;

        // Copy number for each transposon type
        if (NTETypes() > 1) for (int t = 0; t < NTETypes(); ++t, c += 4)
        {
            df.Set(row, c + 0, _n[t].Mean());
            df.Set(row, c + 1, _n[t].Variance());
            df.Set(row, c + 2, _n[t].Interval05());
            df.Set(row, c + 3, _n[t].Interval95());
        }

        // Rate, suppressed rate, and activity for each transposon family
        if (P.families.size() > 1) for (int f = 0; f < (int)P.families.size(); ++f, c += 12)
        {
            df.Set(row, c + 0, _u[f].Mean());
            df.Set(row, c + 1, _u[f].Variance());
            df.Set(row, c + 2, _u[f].Interval05());
            df.Set(row, c + 3, _u[f].Interval95());

            df.Set(row, c + 4, _U[f].Mean());
            df.Set(row, c + 5, _U[f].Variance());
            df.Set(row, c + 6, _U[f].Interval05());
            df.Set(row, c + 7, _U[f].Interval95());

            df.Set(row, c + 8, _a[f].Mean());
            df.Set(row, c + 9, _a[f].Variance());
            df.Set(row, c + 10, _a[f].Interval05());
            df.Set(row, c + 11, _a[f].Interval95());
        }

        // Suppression for each suppressor type
        if (P.s0.size() > 1) for (int s = 0; s < (int)P.s0.size(); ++s, c += 4)
        {
            df.Set(row, c + 0, _S[s].Mean());
            df.Set(row, c + 1, _S[s].Variance());
            df.Set(row, c + 2, _S[s].Interval05());
            df.Set(row, c + 3, _S[s].Interval95());
        }

        avg_n = _nn.Mean();
        avg_u = _uu.Mean();
        avg_S = _SS.Mean();
        avg_w = _w.Mean();

        ++row;
    }

    ++g;
}

void Recorder::Flush(string filename)
{
    ostringstream oss;
    oss << "incomplete run -- flushed in sweep " << P.SweepName() << ", rep " << rep << ", generation " << g << "\n\n";
    P.Write(oss);
    df.Save(filename);
}

void Recorder::EndRep()
{
    // Select next rep for datafile, if there is another left
    if (rep + 1 < P.reps)
        df.SelectRep(rep + 1);

    // Reset for next rep
    Reset();

    ++rep;
}

void Recorder::EndRun()
{
    // CHECKSUM
    if (P.do_checksum)
    {
        ostringstream ocs;
        ocs << hex << cs.checksum();

        cout << "\tfinal checksum = " << ocs.str();
        if (P.expected_checksum != ocs.str())
            cout << " (WARNING: Does not match expected checksum \"" << P.expected_checksum << "\")\n";
        else
            cout << " (MATCH)\n";
    }

    // Data file
    if (!P.fileout.empty())
    {
        ostringstream oss;
        P.Write(oss);
        df.ShrinkToFit();
        df.Save(P.fileout, oss.str());
    }
}

void Recorder::Reset()
{
    g = 0;
    row = 0;
    avg_n = avg_u = avg_S = avg_w = 0;
}