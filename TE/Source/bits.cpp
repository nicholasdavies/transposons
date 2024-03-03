// bits.cpp

#include "bits.h"
#include <vector>
#include <algorithm>
#include <stdexcept>
#include <iostream>
#include <type_traits>
using namespace std;

Chromosome::Chromosome()
{
    // Choose number of initial elements
    int n = R.Poisson(P.initial_elements / (P.ploidy * P.haploid_N));
    els.reserve(n);

    // Choose element positions
    vector<int> positions;
    R.Pick(0, P.n_sites - 1, min(P.n_sites, n), positions);

    // Create elements
    for (auto& p : positions)
        els.emplace_back(P.u0, p, 0);
}

void Chromosome::Recombine(const Chromosome& c1, const Chromosome& c2)
{
    // Storage for chiasmata positions
    static vector<int> chiasmata;

    // Choose starting chromosome at random
    const Chromosome* CA = &c1, * CB = &c2;
    if (R.Discrete(0, 1)) swap(CA, CB);

    // Set up: [A_begin, A_end) and [B_begin, B_end) are spans between chiasmata on chromosomes A and B
    auto A_begin = CA->els.begin();
    auto A_end = A_begin;
    auto B_begin = CB->els.begin();
    auto B_end = B_begin;
    els.clear();

    // Pick chiasmata; we start chiasmata at 1 because a chiasma at 0 just exchanges the chromosomes
    if (P.chiasma_cmf.empty())
        R.Pick(1, P.n_sites - 1, min(P.n_sites - 1, R.Poisson(P.chromosome_length)), chiasmata);
    else
        R.Pick(1, P.n_sites - 1, min(P.n_sites - 1, R.Discrete(P.chiasma_cmf)), chiasmata);

    // Do recombination in stages, taking spans between chiasmata from alternating chromosomes
    for (auto chiasma : chiasmata)
    {
        // Update end pointers
        while (A_end != CA->els.end() && (int)A_end->pos < chiasma) ++A_end;
        while (B_end != CB->els.end() && (int)B_end->pos < chiasma) ++B_end;

        // Insert span of elements
        els.insert(els.end(), A_begin, A_end);

        // Swap chromosomes and reset
        A_begin = B_end;
        B_begin = A_end;
        swap(CA, CB);
        A_end = A_begin;
        B_end = B_begin;
    }

    // Insert remaining TEs
    els.insert(els.end(), A_begin, CA->els.end());
}

Individual::Individual()
 : chr(P.ploidy * P.haploid_N), sbv(P.s0.size(), 0)
{
    // Initialise suppressors
    if (P.ploidy == 1)
        sg = P.s0;
    else
        for (int i = 0; i < (int)P.s0.size(); ++i)
            sg.insert(sg.end(), 2, P.s0[i]);

    // Allocate space for types & families
    n_type.resize(NTETypes(), 0);
    a_fam.resize(P.families.size(), 0);

    // Initialise all summary statistics
    RefreshPostBirth();
    RefreshPostDevelopment();
}

void Individual::Print() const
{
    for (auto& c : chr)
    {
        int last = -1;
        for (auto& e : c.els)
        {
            cout << string(e.pos - last - 1, '-') << e.u;
            last = e.pos;
        }
        cout << string(P.n_sites - last - 1, '-') << "\n";
    }
}

void Individual::Conceive(const Individual& p1)
{
    *this = p1;
}

void Individual::Conceive(const Individual& p1, const Individual& p2)
{
    // For haploids, recombine parental chromosomes into child; choose suppressor genes from random parent
    if (P.ploidy == 1)
    {
        for (int c = 0; c < P.haploid_N; ++c)
            chr[c].Recombine(p1.chr[c], p2.chr[c]);
        for (int i = 0; i < (int)P.s0.size(); ++i)
            sg[i] = R.Discrete(0, 1) ? p1.sg[i] : p2.sg[i];
    }

    // For diploids, recombine parental chromosomes into child; choose one of each parent's suppressor genes
    else
    {
        for (int c = 0; c < P.haploid_N; ++c)
        {
            chr[c * 2    ].Recombine(p1.chr[c * 2], p1.chr[c * 2 + 1]);
            chr[c * 2 + 1].Recombine(p2.chr[c * 2], p2.chr[c * 2 + 1]);
        }
        for (int i = 0; i < (int)P.s0.size(); ++i)
        {
            int sel = R.Discrete(0, 3);
            sg[2 * i]     = (sel & 1) ? p1.sg[2 * i] : p1.sg[2 * i + 1];
            sg[2 * i + 1] = (sel & 2) ? p2.sg[2 * i] : p2.sg[2 * i + 1];
        }
    }
}

void Individual::RefreshPostBirth()
{
    n_te = 0;
    for (auto& c : chr)
        n_te += c.els.size();
    if (P.ploidy == 1)
        sbv = sg;
    else for (int l = 0; l < (int)sbv.size(); ++l)
        sbv[l] = (sg[2 * l] + sg[2 * l + 1]) / 2.0;
}

void Individual::RefreshPostDevelopment()
{
    activity = 0.0;
    n_type.assign(n_type.size(), 0);
    a_fam.assign(P.families.size(), 0.0);

    // Keeping track of zygosity
    if (P.zygosity)
    {
        // Iterate through elements on homologous pairs of chromosomes
        for (unsigned int c1 = 0, c2 = 1; c1 < chr.size(); c1 += 2, c2 += 2)
        {
            for (unsigned int e1 = 0, e2 = 0; e1 < chr[c1].els.size() || e2 < chr[c2].els.size(); )
            {
                if (e1 < chr[c1].els.size() && (e2 >= chr[c2].els.size() || chr[c1].els[e1].pos < chr[c2].els[e2].pos))
                {
                    // Classify e1 on c1 as a heterozygote
                    activity += chr[c1].els[e1].u;
                    a_fam[chr[c1].els[e1].id] += chr[c1].els[e1].u;
                    ++n_type[TEType(chr[c1].els[e1].id, P.site_pattern[chr[c1].els[e1].pos % P.site_pattern.size()], false)];
                    ++e1;
                }
                else if (e2 < chr[c2].els.size() && (e1 >= chr[c1].els.size() || chr[c1].els[e1].pos > chr[c2].els[e2].pos))
                {
                    // Classify e2 on c2 as a heterozygote
                    activity += chr[c2].els[e2].u;
                    a_fam[chr[c2].els[e2].id] += chr[c2].els[e2].u;
                    ++n_type[TEType(chr[c2].els[e2].id, P.site_pattern[chr[c2].els[e2].pos % P.site_pattern.size()], false)];
                    ++e2;
                }
                else
                {
                    // Classify both e1 on c1 and e2 on c2 as homozygotes
                    activity += chr[c1].els[e1].u + chr[c2].els[e2].u;
                    a_fam[chr[c1].els[e1].id] += chr[c1].els[e1].u;
                    a_fam[chr[c2].els[e2].id] += chr[c2].els[e2].u;
                    ++n_type[TEType(chr[c1].els[e1].id, P.site_pattern[chr[c1].els[e1].pos % P.site_pattern.size()], true)];
                    ++n_type[TEType(chr[c2].els[e2].id, P.site_pattern[chr[c2].els[e2].pos % P.site_pattern.size()], true)];
                    ++e1, ++e2;
                }
            }
        }
    }

    // Not keeping track of zygosity - this separate code is needed for haploids
    else
    {
        for (auto& c : chr)
        {
            for (auto& e : c.els)
            {
                activity += e.u;
                a_fam[e.id] += e.u;
                ++n_type[TEType(e.id, P.site_pattern[e.pos % P.site_pattern.size()], false)];
            }
        }
    }

    // Calculate fitness
    for (auto& n : n_type)
        P.w << n;
    P.w << S();
    for (auto& a : a_fam)
        P.w << a;
    fitness = max(0.0, P.w.Call());
}

Element& Individual::NthElement(unsigned int n)
{
    for (auto& c : chr)
    {
        if (n < c.els.size())
            return c.els[n];
        n -= c.els.size();
    }
    throw runtime_error("Attempt to access invalid TE.\n");
}

void Individual::MutateUNthElement(unsigned int n)
{
    NthElement(n).MutateU();
}

void Individual::ExciseNthElement(unsigned int n)
{
    for (auto& c : chr)
    {
        if (n < c.els.size())
        {
            --n_te;
            c.els.erase(c.els.begin() + n);
            return;
        }
        n -= c.els.size();
    }
    throw runtime_error("Attempt to remove invalid TE.\n");
}

bool Individual::InsertElement(Element& dup)
{
    unsigned int c = 0, p = 0;
    vector<Element>::iterator x;

    // If permitted to overwrite, choose a random site.
    if (P.te_overwrite)
    {
        int site = R.Discrete(0, P.ploidy * P.haploid_N * P.n_sites - 1);
        c = site / P.n_sites;
        p = site % P.n_sites;
        x = lower_bound(chr[c].els.begin(), chr[c].els.end(), p, [](auto el, auto p) { return el.pos < p; });
    }
    // If not permitted to overwrite, find an empty site, if one exists
    else
    {
        int empty_spaces = P.ploidy * P.haploid_N * P.n_sites - n_te, empty_spaces_here = 0;
        if (empty_spaces <= 0)
            return false;

        // Insert in the nth empty space
        int n = R.Discrete(0, empty_spaces - 1);

        // Locate the correct chromosome
        for (c = 0; c < chr.size(); ++c)
        {
            empty_spaces_here = P.n_sites - chr[c].els.size();
            if (n - empty_spaces_here < 0) break;
            n -= empty_spaces_here;
        }

        // Find the nth empty space, which is on chromosome c, or a random empty space on chromosome c,
        // depending on which method is likely to be faster.

        // This tuning parameter is approximately 2B/A, where A is mean time needed for one inner loop of the first method
        // and B is mean time needed for one inner loop of the second method; in total the first method takes about
        // A * (n_sites - empty_spaces_here) / 2, and the second method takes about B * n_sites / empty_spaces_here.
        // This particular value is a rough guess, made by timing a few different values.
        const double tuning = 30;

        // Find by sequential search
        if (P.n_sites * tuning > (P.n_sites - empty_spaces_here) * empty_spaces_here)
        {
            int prev = -1;
            for (x = chr[c].els.begin(); x != chr[c].els.end(); ++x)
            {
                int skipped = x->pos - prev - 1;
                n -= skipped;
                if (n < 0)
                {
                    p = (int)(x->pos) + n;
                    break;
                }
                prev = x->pos;
            }
            if (n >= 0)
            {
                if (chr[c].els.empty())
                    p = n;
                else
                    p = chr[c].els.back().pos + n + 1;
                x = chr[c].els.end();
            }
        }
        // Find by sampling and rejection
        else
        {
            while (true)
            {
                p = R.Discrete(0, P.n_sites - 1);
                x = lower_bound(chr[c].els.begin(), chr[c].els.end(), p, [](auto el, auto p) { return el.pos < p; });
                if (x == chr[c].els.end() || x->pos != p)
                    break;
            }
        }
    }

    // Overwrite existing
    if (x != chr[c].els.end() && x->pos == p)
    {
        x->u = dup.u;
        x->id = dup.id;
    }
    // Create new
    else
    {
        ++n_te;
        auto new_element = chr[c].els.insert(x, dup);
        new_element->pos = p;
    }

    return true;
}

void Individual::MutateSuppressor(unsigned int i)
{
    int n_s = i / P.ploidy;
    double new_s = Clamp(P.s_min[n_s], sg[i] + R.Normal(0, P.s_mut_spread[n_s] + sg[i] * P.s_mut_spread_relative[n_s]), P.s_max[n_s]);
    sbv[n_s] += (new_s - sg[i]) / P.ploidy;
    sg[i] = new_s;
}

void Individual::Invade(double proportion, double u, int id)
{
    for (auto& c : chr)
    {
        for (auto& e : c.els)
        {
            if (R.Bernoulli(proportion))
            {
                if (u >= 0) e.u = u;
                if (id >= 0) e.id = id;
            }
        }
    }
}

double Individual::DuplicationRates(vector<double>& rates) const
{
    // 1. Calculate suppression factors
    static vector<double> sf(P.s0.size());
    sf.assign(P.s0.size(), 0.0);
    for (int i = 0; i < (int)P.s0.size(); ++i)
        sf[i] = P.U(1.0, sbv[i]);

    // 2. Get sum of duplication rates and counts for each family
    static vector<double> a(P.families.size());
    static vector<int> n(P.families.size());
    a.assign(P.families.size(), 0.0);
    n.assign(P.families.size(), 0);
    for (auto& c : chr)
    {
        for (auto& e : c.els)
        {
            a[e.id] += e.u;
            n[e.id]++;
        }
    }

    // 3. Collect parasitism rates for each family (total suppressed activity of any parasites on a given family)
    static vector<double> par(P.families.size());
    par.assign(P.families.size(), 0.0);
    for (int f = 0; f < (int)P.families.size(); ++f)
        if (P.families[f] <= 0)
            par[-P.families[f]] += sf[f % P.s0.size()] * a[f];

    // 4. Calculate individual duplication rate for each transposon
    rates.resize(n_te, 0);
    double sum = 0;
    unsigned int i = 0;

    // Fill duplication rates vector
    for (auto& c : chr)
    {
        for (auto& e : c.els)
        {
            double u = 0.0;                         // duplication rate for this element
            double SF = sf[e.id % P.s0.size()];     // suppression factor acting on this element
            double SFa;
            int target = 0;

            switch (P.families[e.id])
            {
                case 1:     // Public
                    u = SF * a[e.id] / n[e.id];
                    break;
                case 2:     // Private
                    ///u = SF * e.u * (1.0 - P.par_pref + P.par_pref * P.par_hs / (P.par_hs + par[e.id]));
                    u = SF * e.u * (1.0 - P.par_pref * par[e.id] / (P.par_hs + par[e.id]));
                    break;
                default:    // Parasitic
                    target = -P.families[e.id];     // "host" TE family
                    SFa = sf[target % P.s0.size()];    // suppression factor acting on "host" TE
                    ///u = SFa * SF * a[target] * P.par_pref * e.u / (P.par_hs + a[e.id]);
                    u = SFa * a[target] * (P.par_pref * SF * e.u) / (P.par_hs + SF * a[e.id]);
                    break;
            }

            sum += u;
            rates[i++] = sum;
        }
    }

    return sum;
}

#define checksum_cast(x) reinterpret_cast<char*>(&x)

void Individual::ProcessChecksum(boost::crc_32_type& cs)
{
    for (auto& c : chr)
        for (auto& e : c.els)
            cs.process_bytes(checksum_cast(e), sizeof e);

    for (int i = 0; i < (int)sg.size(); ++i)
        cs.process_bytes(checksum_cast(sg[i]), sizeof sg[i]);
}

// Helpers for tracking transposon types
int NTETypes()
{
    return P.families.size() * P.site_classes * (P.zygosity ? 2 : 1);
}

int TEType(int family, int site_class, bool homozygous)
{
    return family * P.site_classes * (P.zygosity ? 2 : 1)
        + site_class * (P.zygosity ? 2 : 1)
        + (P.zygosity && homozygous ? 1 : 0);
}

string TETypeName(int t)
{
    string name = "n";

    if (P.zygosity)
        name += (t % 2 == 0 ? 'h' : 'g');

    if (P.families.size() > 1)
        name += to_string(int(t / (P.site_classes * (P.zygosity ? 2 : 1))));

    if (P.site_classes > 1)
        name += char('A' + (t % (P.site_classes * (P.zygosity ? 2 : 1))) / (P.zygosity ? 2 : 1));

    return name;
}
