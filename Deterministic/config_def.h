// config_def.h
// contains definitions of all config parameters to be loaded from the
// configuration file or the command line. Parameter names must be valid C++
// identifiers. Parameter names starting with a lowercase letter are
// guaranteed not to conflict with other identifiers.

// Experimental parameters
PARAMETER ( int,                    generations,            5000 );
PARAMETER ( int,                    report_freq,            1000 );
PARAMETER ( int,                    n_max,                  100 );
PARAMETER ( string,                 mode,                   "static" );
PARAMETER ( string,                 type,                   "RNA" );
PARAMETER ( double,                 n_0,                    1.0 );
PARAMETER ( LuaFunction<double>,    w,                      "n -> 1 - 0.001 * n ^ 1.5" );   // can take: n; na,nb; n,u; na,ua,nb,ub
PARAMETER ( double,                 u,                      0.01 );
PARAMETER ( double,                 u_B,                    0.0101 );
PARAMETER ( double,                 v,                      0.005 );
PARAMETER ( double,                 epsilon,                1e-9 );
PARAMETER ( double,                 gamete_threshold,       1e-12 );
PARAMETER ( double,                 transition_threshold,   1e-12 );
PARAMETER ( double,                 marginal_threshold,     1e-16 );
PARAMETER ( double,                 invasion_seed,          0.01 );
PARAMETER ( double,                 invasion_threshold,     0.5 );
PARAMETER ( vector<int>,            sequential_g,           vector<int>     ({ 1000, 1000 }) );
PARAMETER ( vector<double>,         sequential_u,           vector<double>  ({ 0.01, 0.02 }) );
PARAMETER ( vector<int>,            sequential_n_max,       vector<int>     () );
PARAMETER ( int,                    analysis_gens,          100 );
PARAMETER ( int,                    analysis_gens0,         -1 );
PARAMETER ( int,                    analysis_genotype_report_freq, 1 );
PARAMETER ( double,                 analysis_seed,          1e-9 );
PARAMETER ( int,                    analysis_nb_max,        -1 );
PARAMETER ( vector<int>,            analysis_u_change_gens, vector<int> () );
PARAMETER ( double,                 analysis_u_change_amt,  0.0 );
PARAMETER ( int,                    analysis_special_gen,   -1 );
PARAMETER ( string,                 analysis_special_focus, "none" );
PARAMETER ( double,                 analysis_special_prob,  0.0 );
PARAMETER ( double,                 poisson_u_step,         0.01 );
PARAMETER ( string,                 fileout,                "./output.ngd" );
PARAMETER ( string,                 fileout_genotypes,      "" );
