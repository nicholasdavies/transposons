library(data.table)
library(ggplot2)
library(stringr)
library(Rcpp)



####################
# LOAD SIMULATIONS #
####################


source("./load_ngd.R")
source("./cost_benefit_functions.R")
sourceCpp("./cost_benefit_functions.cpp")

# Get the full list of analysis scenarios
analysis_files = list.files("./Deterministic/Runs/4-Analysis", 
    "^(pub|pri)_c([0-9]+)_u([0-9]+)_up\\.ngd$", full.names = TRUE);
analysis_matches = str_match(analysis_files, 
    "(pub|pri)_c([0-9]+)_u([0-9]+)_up\\.ngd$");
analyses = as.data.table(analysis_matches);
names(analyses) = c("file", "mode", "cost", "dupl");

# Load all analysis scenarios from the deterministic simulation
ds = list();
for (r in 1:nrow(analyses)) {
    cat(".");
    ds[[r]] = load_analysis(analyses[r, mode], analyses[r, cost], analyses[r, dupl]);
}
ds = rbindlist(ds, idcol = "sim_id");

# Reference parameters
GE = 200 # "End" generation
GC = 100 # "Change" generation for analyses with one duplication event
DD = 0.0001 # Incremental duplication rate/probability for analyses

# Transform analyses for no-selection simulations
analyses[, cost := as.numeric(paste0("0.", cost))];
analyses[, dupl := as.numeric(paste0("0.", dupl))];
analyses = merge(analyses,
    ds[stage == 3 & g == GE & tag == "new", .(mode, cost, dupl, A_B)],
    by = c("mode", "cost", "dupl"));

# Choose max. genome size for no-selection simulations, to keep 
# relative error in target_A_B below 1e-4 (see below for check).
analyses[, N_ns := 50]
analyses[cost == 5e-3, N_ns := pmax(50, round(A_B * 2.7))]
analyses[cost == 5e-4, N_ns := pmax(50, round(A_B * 2.0))]

# This file can also be used from the command line to just generate the
# no-selection simulations, if the first and last simulation ids are
# passed.
generate = FALSE;
args = commandArgs(trailingOnly = TRUE);
if (length(args) == 2) {
    r_start = as.integer(args[1]);
    r_end = as.integer(args[2]);
    generate = TRUE;
} else {
    r_start = 1;
    r_end = nrow(analyses);
}

# Load and/or conduct all no-selection simulations
ns = list();

for (r in r_start:r_end) {
    cat(paste(analyses[r, file], analyses[r, N_ns]));

    # Update those simulations that require it
    src_file = paste0("./Deterministic/Runs/4-Analysis/", analyses[r, file]);
    dst_file = paste0("./Cost-Benefit/", analyses[r, file]);
    src_time = file.mtime(src_file);
    dst_time = file.mtime(dst_file);

    N = analyses[r, N_ns];
    u = analyses[r, dupl];
    A_B = analyses[r, A_B];
    mode = analyses[r, mode];

    if (generate && (is.na(dst_time) || src_time > dst_time)) {
        sim_new = simulate(N, GE, GC, u, DD, A_B, mode, "new");
        sim_new[, tag := "new"];
        sim_old = simulate(N, GE, GC, u, DD, A_B, mode, "old");
        sim_old[, tag := "old"];
        sim = rbind(sim_new, sim_old);
        sim[, cost := analyses[r, cost]];
        fwrite(sim, dst_file);
    } else {
        cat(" -- skipping.\n");
    }
    
    # Reloading from file so there is no difference between just-generated and previously-generated simulations
    if (!is.na(dst_time)) {
        sim = fread(dst_file);
    } else {
        sim = NULL;
    }

    # Save simulation in list
    ns[[r]] = sim;
}
if (length(args) == 2) {
    quit();
}
ns = rbindlist(ns, idcol = "sim_id")
