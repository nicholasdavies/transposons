// config_def.h
// contains definitions of all config parameters to be loaded from the
// configuration file or the command line. Parameter names must be valid C++
// identifiers. Parameter names starting with a lowercase letter are
// guaranteed not to conflict with other identifiers.

// Experimental parameters
PARAMETER ( int, N, 100 );
PARAMETER ( int, generations, 1000 );
PARAMETER ( int, reps, 1 );
PARAMETER ( int, n_sites, 31 );
PARAMETER ( int, ploidy, 2 );
PARAMETER ( int, haploid_N, 2 );
PARAMETER ( double, chromosome_length, 0.9 );   // in morgans
PARAMETER ( vector<double>, chiasma_cmf, vector<double>{} );
PARAMETER ( double, sexuality, 1.0 );   // 1.0 = fully sexual, 0.0 = fully clonal
PARAMETER ( double, selfing, 0.0 ); // probability of mating with self, if reproducing sexually
PARAMETER ( double, extinction_threshold, 0.0 );    // extinction happens when average fitness drops to or below this value
PARAMETER ( bool, halt_on_clearance, false );   // stop simulation if TEs disappear
PARAMETER ( bool, halt_on_n_max, false );   // stop simulation if avg. TE number exceeds n_max
PARAMETER ( bool, te_overwrite, true ); // allow TEs to overwrite occupied sites. if false, a duplicating TE will "seek out" a free slot

// TE properties
PARAMETER ( double, initial_elements, 10 );
PARAMETER ( int, n_max, 0 ); // if > 0, do not allow number of TEs to go above that number
PARAMETER ( bool, zygosity, false );    // if true, keep track of zygosity of elements
PARAMETER ( int, site_classes, 1 );     // number of site classes
PARAMETER ( vector<int>, site_pattern, vector<int>{0} );   // distribution of site classes
PARAMETER ( vector<int>, families, vector<int>{2} );   // list of family attributes.
// for the above, 1 = public duplication; 2 = private duplication; -n = parasite of private family n, pointing into this vector (0-based)

// Duplication
PARAMETER ( double, u0, 0.01 );
PARAMETER ( double, u_mut_rate, 0.0 ); // per element per generation
PARAMETER ( double, u_mut_rate_on_copy, 0.0 ); // per element duplication
PARAMETER ( vector<double>, u_mut_spread, vector<double>{0.0} ); // standard deviation of mutation deltas - absolute component. indexed by family % size
PARAMETER ( vector<double>, u_mut_spread_relative, vector<double>{0.0} ); // standard deviation of mutation deltas - relative component. indexed by family % size
PARAMETER ( double, u_min, 0.0 ); // lower limit on per-element u.
PARAMETER ( double, u_max, 5.0 ); // upper limit on per-element u.
PARAMETER ( LuaFunction<double>, U, "u, S -> u" ); // Element duplication rate, host's suppression

// Cis-preference / parasitizing
PARAMETER ( double, par_hs, 1.0 ); // half-saturation constant for parasites
PARAMETER ( double, par_pref, 0.5 ); // if a private mRNA and a parasite mRNA are competing for the same ribosome/transposon, the probability that the parasite gets picked.

// Excision
PARAMETER ( double, v, 0.005 );

// Suppressor properties
PARAMETER ( vector<double>, s0, vector<double>{0.0} );  // initial suppression amount
PARAMETER ( vector<double>, s_mut_rate, vector<double>{0.0} );  // per gene per generation
PARAMETER ( vector<double>, s_mut_spread, vector<double>{0.0} );    // standard deviation of mutation deltas - absolute component
PARAMETER ( vector<double>, s_mut_spread_relative, vector<double>{0.0} );    // standard deviation of mutation deltas - relative component
PARAMETER ( vector<double>, s_min, vector<double>{0.0} );   // lower limit on s gene
PARAMETER ( vector<double>, s_max, vector<double>{1.0} );   // upper limit on s gene

// Fitness
PARAMETER ( LuaFunction<double>, w, "n, S, a -> 1 - 0.001 * n ^ 1.5" ); // copy number, suppression amount, overall activity

// Manipulation
PARAMETER ( int, gen_coevolution, 0 );  // generation at which to enable TE mutation and suppressor recruitment/loss
PARAMETER ( vector<double>, invasions, vector<double>{} );  // List of [gen, freq, fam, u, ]* : at generation gen

// Checksum and testing-related stuff
PARAMETER ( bool, do_checksum, false );
PARAMETER ( string, expected_checksum, "" );
PARAMETER ( bool, check, false );
PARAMETER ( bool, test, false );

// Recording
PARAMETER ( string, fileout, "output.ngd" );
PARAMETER ( int, report_freq, 1 );  // record every report_freq generations
