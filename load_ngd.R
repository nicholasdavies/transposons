library(data.table)
library(ggplot2)

# data file entries:
# n -- copy number, over hosts
# u -- breeding value for activity, over transposons
# U -- actual duplication rate (accounting for suppression, parasitism, etc), over transposons
# a -- sum of non-adjusted u, over hosts
# S -- breeding value for suppression, over hosts
# w -- absolute fitness, over hosts
# x -- host relative fitness, over transposons
# d -- proportion of hosts who perish each generation (i.e. with fitness 0) -- this will be 0 if fitness function never reaches 0
# G -- homozygosity at TE insertion sites, over hosts with at least one TE
# BuW -- slope of host's relative fitness on u, over transposons
# Bun -- slope of host's copy number on u, over transposons

# load data file
load_ngd = function(path, thin = 100, gunit = 1000, tag = NULL)
{
    # Get line numbers beginning with {DATA}
    command = paste0('grep -n "^{DATA}" ', path, ' | cut -d ":" -f 1')
    skip = as.numeric(system(command, intern = TRUE))

    # Number of rows to read for each {DATA} block
    nrows = c(diff(skip) - 3, Inf)

    # Read each data block
    dt = list();
    for (i in seq_along(skip)) {
        dt[[i]] = fread(path, skip = skip[i], nrows = nrows[i], 
            blank.lines.skip = TRUE, colClasses = "numeric");
    }
    dt = rbindlist(dt, idcol = "run");

    # Thin etc
    dt = dt[g %% thin < 1];
    dt[, g := as.numeric(g) / gunit];
    
    # Set NaNs to NAs
    for (col in names(dt)) {
        dt[is.nan(get(col)), (col) := NA]
    }
    
    if (!is.null(tag)) {
        dt[, tag := ..tag]
    }
    
    return (dt[]);
}
