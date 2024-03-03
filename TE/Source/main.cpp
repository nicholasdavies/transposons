// A model

#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <algorithm>
#include "shared.h"
#include "experiment.h"
#include "string.hpp"
using namespace std;

void ValidateParameters();
void Test();

int main(int argc, char* argv[])
{
    // Load parameters from command line, with optional specifying of config file as first parameter
    P.Read(argc, argv);

    // Perform coherence check
    if (P.check)
    {
        cout << "Test 1\n";
        P.Read("./Check/test1.cfg");
        Experiment exp1;
        R.Reset();

        return 0;
    }

    // Perform tests
    if (P.test)
    {
        Test();
        return 0;
    }

    double c0 = Clock();

    for (int s = 0; s < P.NSweeps(); ++s)
    {
        cout << "SWEEP " << P.Sweep() << ", [" << P.SweepName() << "]\n";
        P.Write(cout);
        cout << "\n";

        ValidateParameters();

        Experiment exp;
        if (Signal) break;
        P.NextSweep();
        R.Reset();
    }

    cout << "\nSimulation took " << Clock() - c0 << " seconds.\n";
}

// Ensure Lua functions in parameters have the right arguments.
template <typename LuaFunc>
bool HasArguments(LuaFunc& func, string args)
{
    string comp = func.CompactForm();
    comp = comp.substr(0, comp.find("->"));
    return replace(comp, " ", "") == replace(args, " ", "");
}

// Ensure parameters are sensible.
void ValidateParameters()
{
    string w_arguments;
    ostringstream n_types;

    for (int t = 0; t < NTETypes(); ++t)
        n_types << TETypeName(t) << ", ";

    w_arguments = n_types.str() + "S";

    for (int f = 0; f < (int)P.families.size(); ++f)
        w_arguments += ", a" + (P.families.size() > 1 ? to_string(f) : string());

    if (P.N <= 0)
        throw runtime_error("Error: N must be greater than 0.\n");
    if (P.generations <= 0)
        throw runtime_error("Error: generations must be greater than 0.\n");
    if (P.reps <= 0)
        throw runtime_error("Error: reps must be greater than 0.\n");
    if (P.n_sites <= 0)
        throw runtime_error("Error: n_sites must be greater than 0.\n");
    if (P.ploidy != 1 && P.ploidy != 2)
        throw runtime_error("Error: ploidy must be either 1 (haploid) or 2 (diploid).\n");
    if (P.haploid_N <= 0)
        throw runtime_error("Error: haploid_N must be greater than 0.\n");
    if (P.chromosome_length < 0)
        throw runtime_error("Error: chromosome_length must be nonnegative.\n");
    if (!P.chiasma_cmf.empty() && (*min_element(P.chiasma_cmf.begin(), P.chiasma_cmf.end()) < 0 ||
                                   !is_sorted(P.chiasma_cmf.begin(), P.chiasma_cmf.end())))
        throw runtime_error("Error: chiasma_cmf must be an ordered list of nonnegative numbers.\n");
    if (P.sexuality < 0 || P.sexuality > 1)
        throw runtime_error("Error: sexuality must be between 0 and 1.\n");
    if (P.extinction_threshold < 0 || P.extinction_threshold > 1)
        throw runtime_error("Error: extinction_threshold must be between 0 and 1.\n");
    if (P.initial_elements < 0)
        throw runtime_error("Error: initial_elements must be nonnegative.\n");
    if (P.zygosity && P.ploidy != 2)
        throw runtime_error("Error: can only keep track of zygosity in diploids.\n");
    if (P.site_classes <= 0)
        throw runtime_error("Error: site_classes must be greater than 0.\n");
    if (P.site_pattern.empty() ||
         *min_element(P.site_pattern.begin(), P.site_pattern.end()) < 0 ||
         *max_element(P.site_pattern.begin(), P.site_pattern.end()) >= P.site_classes)
        throw runtime_error("Error: site_pattern must be a list of integers between 0 and site_classes - 1.\n");
    if (P.families.empty() ||
         *min_element(P.families.begin(), P.families.end()) <= -(int)P.families.size() ||
         *max_element(P.families.begin(), P.families.end()) > 2)
        throw runtime_error("Error: families must be a list of integers between -size(families) + 1 and 2.\n");
    if (P.u0 < 0)
        throw runtime_error("Error: u0 must be nonnegative.\n");
    if (P.u_mut_rate < 0)
        throw runtime_error("Error: u_mut_rate must be nonnegative.\n");
    if (P.u_mut_rate_on_copy < 0)
        throw runtime_error("Error: u_mut_rate_on_copy must be nonnegative.\n");
    if (P.u_mut_spread.size() != P.u_mut_spread_relative.size())
        throw runtime_error("Error: u_mut_spread and u_mut_spread_relative must be the same size.\n");
    for (int i = 0; i < (int)P.u_mut_spread.size(); ++i)
    {
        if (P.u_mut_spread[i] < 0)
            throw runtime_error("Error: all u_mut_spread must be nonnegative.\n");
        if (P.u_mut_spread_relative[i] < 0)
            throw runtime_error("Error: all u_mut_spread_relative must be nonnegative.\n");
    }
    if (P.u_min < 0 || P.u_max < 0 || P.u_min > P.u_max)
        throw runtime_error("Error: u_min and u_max must be nonnegative and ordered.\n");
    if (!HasArguments(P.U, "u, S"))
        throw runtime_error("Error: U must have parameters u, S.\n");
    if (P.par_hs <= 0)
        throw runtime_error("Error: par_hs must be positive.\n");
    if (P.v < 0)
        throw runtime_error("Error: v must be nonnegative.\n");
    if (P.s0.size() != P.s_mut_rate.size()
     || P.s0.size() != P.s_mut_spread.size()
     || P.s0.size() != P.s_mut_spread_relative.size()
     || P.s0.size() != P.s_min.size()
     || P.s0.size() != P.s_max.size())
        throw runtime_error("Error: s0, s_mut_rate, s_mut_spread, s_mut_spread_relative, s_min, and s_max must all be the same size.\n");
    for (int i = 0; i < (int)P.s0.size(); ++i)
    {
        if (P.s0[i] < 0)
            throw runtime_error("Error: all s0 must be nonnegative.\n");
        if (P.s_mut_rate[i] < 0)
            throw runtime_error("Error: all s_mut_rate must be nonnegative.\n");
        if (P.s_mut_spread[i] < 0)
            throw runtime_error("Error: all s_mut_spread must be nonnegative.\n");
        if (P.s_mut_spread_relative[i] < 0)
            throw runtime_error("Error: all s_mut_spread_relative must be nonnegative.\n");
        if (P.s_min[i] < 0 || P.s_max[i] < 0 || P.s_min[i] > P.s_max[i])
            throw runtime_error("Error: all s_min and s_max must be nonnegative and ordered.\n");
    }
    if (!HasArguments(P.w, w_arguments))
        throw runtime_error("Error: w must have parameters " + w_arguments + ".\n");
    if (P.invasions.size() % 4 != 0)
        throw runtime_error("Error: invasions must contain 4*n items.\n");
    if (P.report_freq <= 0)
        throw runtime_error("Error: report_freq must be greater than 0.\n");

    // Show families
    cout << "Families: ";
    for (unsigned int fam = 0; fam < P.families.size(); ++fam)
    {
        if (P.families[fam] == 1)
            cout << "Public" << fam << " ";
        else if (P.families[fam] == 2)
            cout << "Private" << fam << " ";
        else if (P.families[-P.families[fam]] == 2)
            cout << "Parasite[Private" << -P.families[fam] << "] ";
        else
            throw runtime_error("Error: invalid family specification. (parasite on non-private family, perhaps?)\n");
    }
    cout << "\n";
}

// Test various slightly more complex routines
void Test()
{
    P.te_overwrite = false;
    P.n_sites = 100;
    P.ploidy = 2;
    P.haploid_N = 2;
    P.chromosome_length = 4;
    P.u_mut_spread = vector<double>{10};

    P.u0 = 1;
    Individual ind1;
    P.u0 = 2;
    Individual ind2;

    cout << "Individual 1:\n";
    ind1.Print();

    cout << "\nIndividual 2:\n";
    ind2.Print();

    cout << "\nTest of mutation:\n";
    Individual ind3 = ind1;
    for (int i = 0; i < ind3.NElements(); ++i)
        ind3.MutateUNthElement(i);
    ind3.Print();

    cout << "\nTest of excision:\n";
    while (ind3.NElements() > 0)
        ind3.ExciseNthElement(0);
    ind3.Print();

    cout << "\nTest of insertion:\n";
    Individual ind4 = ind2;
    Element el(9, 1, 0);
    int insertions = 0;
    while (ind4.InsertElement(el)) ++insertions;
    cout << "Made " << insertions << " insertions.\n";
    ind4.Print();

    cout << "\nTest of recombination:\n";
    Individual ind5 = ind3, ind6 = ind3;
    Element el5(1, 1, 0), el6(2, 1, 0);
    for (int i = 0; i < 400; ++i)   { ind5.InsertElement(el5); ind6.InsertElement(el6); }
    for (int i = 0; i < 100; ++i)   { ind5.ExciseNthElement(100); ind6.ExciseNthElement(0); }
    for (int i = 0; i < 100; ++i)   { ind5.ExciseNthElement(200); ind6.ExciseNthElement(100); }
    ind5.Print();
    cout << "\n";
    ind6.Print();
    cout << "\n";
    ind1.Conceive(ind5, ind6);
    ind1.Print();

    P.n_sites = 40;
    P.haploid_N = 2;
    P.initial_elements = 20;
    P.zygosity = true;
    cout << "\nTest of zygosity:\n";
    Individual ind7;
    ind7.Print();
    cout << "Individual has " << ind7.NElements(0) << " heterozygote elements and " << ind7.NElements(1) << " homozygote elements.\n";

    P.site_classes = 2;
    P.site_pattern = { 0, 1 };
    P.zygosity = false;
    cout << "\nTest of site classes:\n";
    Individual ind8;
    for (int i = 0; i < 20; ++i)
        cout << "| ";
    cout << "\n";
    ind8.Print();
    cout << "Individual has " << ind8.NElements(0) << " even elements and " << ind8.NElements(1) << " odd elements.\n";

    P.zygosity = true;
    cout << "\nTest of both:\n";
    Individual ind9;
    for (int i = 0; i < 20; ++i)
        cout << "| ";
    cout << "\n";
    ind9.Print();
    cout << "Individual has " << ind9.NElements(0) << " even heterozygote, "
         << ind9.NElements(1) << " even homozygote, "
         << ind9.NElements(2) << " odd heterozygote, and "
         << ind9.NElements(3) << " odd homozygote elements.\n";
}
