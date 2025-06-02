# Graphs comparing "analytical" with simulation results

# Load simulation results
source("./cost_benefit_load.R")
library(ggarchery)
library(ggpattern)
library(scales)
library(cowplot)
theme_set(theme_cowplot(font_size = 8) + theme(plot.title = element_text(size = 8, face = "plain")))



############
# PLOTTING #
############

# Material needed for the plotting is generated after the plotting section
# and saved to CSV files which are loaded for the plotting section.

# Load fitness surface
fitness = fread("./fitness.csv")

# Load evolution thresholds
private_thresh = fread("./private_thresh.csv")   # analytical solutions, though some involve numerically solving for equilibrium n and W
rna_thresh = fread("./rna_thresh.csv")           # numerical solutions
public_thresh = fread("./public_thresh.csv")     # analytical solutions, though some involve numerically solving for equilibrium n and W
dna_thresh = fread("./dna_thresh.csv")           # numerical solutions

# Function to plot the fitness surface
fitness_plot = function(fitness, boundaries, uarr1, uarr2, u_range, cost_range)
{
    f2 = rbind(
        data.table(cost = cost_range,
            threshold = 0,
            u = u_range[2]),
        fitness
    )
    # pal = c("#ffffff", "#3f3391", "#4d56b5", "#718fd0", "#91c0e7");
    # pal = c("#ffffff", "#2B1C8C", "#3F50B0", "#5483D3", "#68B7F7");
    pal = c("#ffffff", "#5141BA", "#5968CE", "#6090E3", "#68B7F7");
    # pal = c("#FFFFFF", "#EA72AD", "#EFA1A1", "#F5D095", "#FAFF89");
    # pal = c("#FFFFFF", "#F4F269", "#C1DD6B", "#8FC76E", "#5CB270");
    # pal = c("#FFFFFF", "#0077b6", "#0096c7", "#00b4d8", "#48cae4");

    pat = c("stripe", "none", "none", "none", "none");
    
    f2[, thresh_name := as.factor(threshold)]
    
    ggplot() +
        geom_ribbon_pattern(data = f2, aes(y = cost, xmin = u_range[1], xmax = u, 
            fill = thresh_name, group = thresh_name, pattern = thresh_name), 
            pattern_fill = "black", pattern_colour = NA, pattern_spacing = 0.01, 
            pattern_density = 0.2, pattern_angle = -45) +
        labs(x = "Duplication rate", y = "Cost") +
        coord_cartesian(xlim = u_range, ylim = cost_range) +
        scale_fill_manual(values = pal) +
        scale_pattern_manual(values = pat) +
        scale_x_log10(expand = c(0, 0), breaks = c(0.0001, 0.001, 0.005, 0.01, 0.05, 0.1, 0.5), labels = c("0.0001", "0.001", "0.005", "0.01", "0.05", "0.1", "0.5")) +
        scale_y_log10(expand = c(0, 0), breaks = 10^(-5:-1), labels = trans_format("log10", math_format(10^.x))) +
        geom_path(data = boundaries, aes(u, cost, colour = what, linetype = what), linewidth = 0.6) +
        annotate("rect", xmin = 9e-5, xmax = 1e-3, ymin = 4e-6, ymax = 1.05e-5, colour = "#dddd99", fill = NA, linetype = "dashed", linewidth = 0.75) +
        annotate("text", x = 3.16e-4, y = 7e-6, hjust = 0.5, vjust = 0.5, colour = "#000000", size = 6/ggplot2:::.pt, label = "TYPICAL RANGE", fontface = "italic") +
        geom_arrowsegment(data = uarr1, aes(x = u, xend = uend, y = cost, yend = cost), 
            arrows = arrow(length = unit(0.15, "cm"), type = "closed"), linewidth = 0.4, colour = "white", fill = "white") +
        geom_arrowsegment(data = uarr2, aes(x = u, xend = uend, y = cost, yend = cost), 
            arrows = arrow(length = unit(0.15, "cm"), type = "closed"), linewidth = 0.4, colour = "white", fill = "white") +
        scale_colour_manual(values = c(extinction = "black", analytical = "black", numerical = "black", simple = "black")) +
        scale_linetype_manual(values = c(extinction = "solid", analytical = "dashed", numerical = "solid", simple = "dotted")) +
        cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "none", plot.title = element_text(size = 8, face = "plain"))
}

# This is just to make boundaries of plot slightly nicer -- does not affect accuracy of the plot.
fitness2 = copy(fitness)
# Extend regions to top of plot
fitness2[cost > 9.97e-2, cost := 0.1];
# Remove one of the sections
fitness2 = fitness2[threshold != 0.01]

# Extend fitness surfaces down to 5*10^-6
thresholds = fitness2[, unique(threshold)];
for (th in thresholds) {
    c_ratio = fitness2[threshold == th, cost[1]/cost[2]];
    u_ratio = fitness2[threshold == th, u[1]/u[2]];
    cost = fitness2[threshold == th, cost[1]];
    u = fitness2[threshold == th, u[1]];
    while (cost > 5e-6) {
        cost = cost * c_ratio;
        u = u * u_ratio;
        fitness2 = rbind(fitness2, data.table(cost = cost, threshold = th, u = u))
    }
}
fitness2 = fitness2[order(-threshold, cost)]

# Assemble key boundaries on fitness plots
pri_boundaries = rbind(
    private_thresh[calc == "product" & cost > 0.001, .(what = "analytical", u, cost)],
    private_thresh[calc == "simple" & cost > 0.0016, .(what = "simple", u, cost)],
    rna_thresh[, .(what = "numerical", u, cost)],
    fitness2[threshold == 1e-6, .(what = "extinction", u, cost)]
)
# Tidy up boundaries for plotting
pri_boundaries[what == "numerical" & cost == min(cost[what == "numerical"]), cost := 0.0022]
pri_boundaries[what == "analytical" & cost == min(cost[what == "analytical"]), cost := 0.0014]
pri_boundaries[what == "simple" & cost == min(cost[what == "simple"]), cost := 0.0021]

# Assemble key boundaries on fitness plots
pub_boundaries = rbind(
    public_thresh[calc == "product", .(what = "analytical", u, cost)],
    public_thresh[calc == "simple", .(what = "simple", u, cost)],
    dna_thresh[, .(what = "numerical", u, cost)],
    public_thresh[calc == "sum" & cost < 5e-5, .(what = "numerical", u, cost)][order(-cost)],
    fitness2[threshold == 1e-6, .(what = "extinction", u, cost)]
)

# Arrows showing evolution
dodge = 0.92
uarr_costs = 6.8 * 10 ^ seq(-6, -2, by = 0.5)
pri_uarr1 = rbind(
    pri_boundaries[what == "extinction" & u < 2/3, .(cost = uarr_costs, uend = approxfun(cost, u)(uarr_costs))],
    pri_boundaries[what == "numerical" & u < 2/3, .(cost = uarr_costs, uend = approxfun(cost, u)(uarr_costs))]
)
pri_uarr1 = pri_uarr1[, .(u = 0.0001/dodge, uend = min(uend, na.rm = TRUE) * dodge), by = .(cost)]

pri_uarr2 = pri_boundaries[what == "numerical", .(cost = uarr_costs, uend = approxfun(cost, u)(uarr_costs))]
pri_uarr2 = pri_uarr2[!is.na(uend), .(u = (2/3)*dodge, uend = min(uend, na.rm = TRUE) / dodge), by = .(cost)]

pub_uarr1 = rbind(
    pub_boundaries[what == "analytical" & u < 2/3, .(cost = uarr_costs, uend = approxfun(cost, u)(uarr_costs))],
    pub_boundaries[what == "numerical" & u < 2/3, .(cost = uarr_costs, uend = approxfun(cost, u)(uarr_costs))]
)
pub_uarr1 = pub_uarr1[, .(u = 0.0001/dodge, uend = max(uend, na.rm = TRUE) * dodge), by = .(cost)]

pub_uarr2 = pub_uarr1[, .(cost, u = (2/3)*dodge, uend = uend / (dodge^2))]
pub_uarr2$u = pub_boundaries[what == "extinction", .(cost, u = pmin(2/3, u))][, approxfun(cost, u)](pub_uarr2$cost) * dodge

pri_uarr1[cost > 6.7e-6 & cost < 6.9e-6, u := 0.0008]
pub_uarr1[cost > 6.7e-6 & cost < 6.9e-6, u := 0.0008]

# Plots
p1 = fitness_plot(fitness2, pub_boundaries[what != "extinction"], pub_uarr1, pub_uarr2, u_range = c(0.0001, 2/3), cost_range = c(5e-6, 0.1)) +
    annotate("text", x = c(0.047, 0.094, 0.19), y = 0.0012, label = c("w = 0.5", "w = 0.01", "w = 0.001"),
        angle = c(56, 52, 50), hjust = 0.5, vjust = 0.5, size = 7 / ggplot2:::.pt, colour = "white", fontface = "bold") +
    labs(title = "Public replication")

p2 = fitness_plot(fitness2, pri_boundaries[what != "extinction"], pri_uarr1, pri_uarr2, u_range = c(0.0001, 2/3), cost_range = c(5e-6, 0.1)) +
    annotate("text", x = c(0.047, 0.094, 0.19), y = 0.0012, label = c("w = 0.5", "w = 0.01", "w = 0.001"),
        angle = c(56, 52, 50), hjust = 0.5, vjust = 0.5, size = 7 / ggplot2:::.pt, colour = "white", fontface = "bold") +
    labs(title = "Private replication")

p = cowplot::plot_grid(p1, p2, nrow = 1, labels = LETTERS, label_size = 8)
ggsave("./Figures/2-invasion.pdf", p, width = 20, height = 10, units = "cm", useDingbats = FALSE)







#############################
# ESTABLISH FITNESS SURFACE #
#############################

# Bounds:
# -- Cost --
# 5x10^-5 is the "believable" cost.
# Try to go down to 1x10^-5.
# Should go up to 0.1 ideally.
# -- Duplication rate --
# In other figures I have down to u = 10^-6
# Probably plot up to u = 1
# Test duplication rates up to extinction - either w = 10^-3 or w = 10^-6

################
# Fitness grid #
################

# This is just to establish an idea of the fitness surface.
# It takes about 15 min to run.
fitness_grid = list();
for (cexp in seq(-5.7, -1.7, by = 0.25)) {
    cat(cexp)
    cost = 5 * 10^cexp;
    for (uexp in seq(-4, -0.75, by = 0.25)) {
        cat(".")
        u = 5 * 10^uexp;
        sim = ext_determ(u, cost, gmax = 10000, epsilon = 1e-12, verbose = FALSE);
        fitness_grid[[length(fitness_grid) + 1]] =
            data.table(cost = cost, u = u, w = sim$w, n = sim$n, Vn = sim$Vn, n_max = sim$n_max);
        if (sim$w < 1e-6 || sim$n_max > 1800) {
            break;
        }
    }
    cat("\n")
}
fitness_grid = rbindlist(fitness_grid)
fitness_grid
# fwrite(fitness_grid, "./fitness_grid.csv")
# fitness_grid = fread("./fitness_grid.csv")

ggplot(fitness_grid) +
    geom_tile(aes(u, cost, fill = floor(log10(w)))) +
    scale_fill_viridis_c() +
    scale_x_log10() +
    scale_y_log10()

ggplot(fitness_grid) +
    geom_tile(aes(u, cost, fill = w)) +
    scale_fill_viridis_c() +
    scale_x_log10() +
    scale_y_log10()

ggplot(fitness_grid) +
    geom_line(aes(u, w, colour = log10(cost), group = cost)) +
    scale_colour_viridis_c() +
    scale_x_log10() +
    scale_y_log10()

# Observations about fitness grid
# For the very highest cost (cost ~ 0.1), things get non-monotonic; up to u = 0.5
# fitness drops to 0.965, but then for u = 0.889 fitness goes back up to 1.
# Doesn't seem to be any weird non-monotonicity for cost = 0.05 up to u = 0.889.

####################
# Fitness contours #
####################

# This section is to get a more exact fitness surface, by finding surface contours.

# Tester for uniroot
test_fitness = function(u, cost, threshold, speak = FALSE)
{
    if (speak) {
        cat(".")
    }
    ext_determ(u, cost, gmax = 10000, nmax = 3000, epsilon = 1e-12, verbose = FALSE)$w / threshold - 1
}

# Find location of fitness threshold in terms of u for a given cost
find_fitness_threshold = function(cost, fitness_grid, threshold)
{
    # Find relevant "cost" entries in fitness_grid to guide selection of u test interval
    cost_grid = fitness_grid[, sort(unique(cost))];
    lb = findInterval(cost, cost_grid);
    bounds = c(lb, lb + 1);
    bounds = bounds[bounds > 0 & bounds <= length(cost_grid)];
    costs = cost_grid[bounds];
    
    # Select u test interval
    min_u = fitness_grid[, min(u)]
    max_u = fitness_grid[, max(u)]
    u_lb = fitness_grid[cost %in% costs & (w > threshold | u == min_u), 
        .(u_lb = max(u)), by = cost][, min(u_lb)]
    u_ub = fitness_grid[cost %in% costs & (w < threshold | u == max_u), 
        .(u_ub = min(u)), by = cost][, max(u_ub)]
    
    # when there's no upper bound, u_ub is -Inf
    if (!is.finite(u_ub)) {
        warning("No upper bound in fitness grid.")
        return (NA);
    }
    
    if (!is.finite(u_lb)) {
        stop("Investigate.")
    }
    
    root = tryCatch(
        uniroot(test_fitness, c(u_lb, u_ub), cost = cost, threshold = threshold, speak = TRUE)$root, 
        error = function(e) { cat(e); cat("\n"); NA })

    cat("\n");
    return (root);
}

# This takes around three hours to run.
fitness = list();
for (cexp in seq(-1.7, -5.7, by = -0.1)) {
    cost = 5 * 10^cexp;
    print(cost);
    for (thresh in c(0.5, 0.1, 0.01, 0.001, 1e-6)) {
        print(thresh)
        fitness[[length(fitness) + 1]] = 
            data.table(cost = cost, threshold = thresh,
                u = find_fitness_threshold(cost, fitness_grid, thresh));
    }
}
fitness = rbindlist(fitness);
fitness = fitness[order(-threshold, cost)]
fitness[, u := ifelse(is.na(u), max(u, na.rm = TRUE), u), by = threshold]

# fwrite(fitness, "./fitness.csv");
# fitness = fread("./fitness.csv")

ggplot(fitness) + 
    geom_point(aes(u, cost, colour = log10(threshold))) +
    scale_x_log10() +
    scale_y_log10()



############################
# FIND INVASION THRESHOLDS #
############################




#######################
# "Analytic" approach #
#######################

# In this part we find invasion thresholds by directly calculating R, B and C
# for a "duplication event".

# Tester for uniroot
# If simulate == FALSE, get a rough equlibrium n and equilibrium dW/dn
# through analytic approximations. If simulate == TRUE, actually simulate
# a 1-type transposon model to get equilibrium n and dW/dn. Note that 
# simulate == TRUE is still part of the "analytic" approach as it then 
# calculates B and C in the normal way.
test_evolve = function(u, cost, simulate = FALSE, how = "product", mode = "private")
{
    cat(u, " ")
    results = get_RBC(u, cost, simulate, how, mode)
    results$R * results$B - results$C
}

# Establishing the thresholds here also takes a while to run.
# Thresholds, private replication
thresh = list();

for (how in c("product", "sum")) {
    for (cexp in seq(-3.8, -2.0, 0.1)) {
        cat("\n", cexp, "\n")
        cost = 5 * 10^cexp;
        thresh[[length(thresh) + 1]] = 
            data.table(cost = cost, eq_n = "simulate", calc = how,
                u = uniroot(test_evolve, c(0.2, 0.5), cost = cost, simulate = TRUE, how = how, mode = "private")$root);
    }
    
    for (cexp in c(seq(-1.9, -1.8, 0.1), -1.69897000433602)) {
        cat("\n", cexp, "\n")
        cost = 5 * 10^cexp;
        thresh[[length(thresh) + 1]] = 
            data.table(cost = cost, eq_n = "simulate", calc = how,
                u = uniroot(test_evolve, c(0.2, 0.5), cost = cost, simulate = TRUE, how = how, mode = "private")$root);
    }
}

for (cexp in c(seq(-6, -1.8, 0.1), -1.69897000433602)) {
    cat("\n", cexp, "\n")
    cost = 5 * 10^cexp;
    thresh[[length(thresh) + 1]] = 
        data.table(cost = cost, eq_n = "analytic", calc = "simple",
            u = uniroot(test_evolve, c(0.001, 0.6), cost = cost, simulate = FALSE, how = "simple", mode = "private")$root);
}
thresh = rbindlist(thresh);

# fwrite(thresh, "./private_thresh.csv")

# Thresholds, public replication
thresh = list();

for (cexp in seq(-6.0, -4.6, 0.1)) {
    cat("\n", cexp, "\n")
    cost = 5 * 10^cexp;
    interval = c(0.5 * sqrt(cost), 0.8 * sqrt(cost)); # Empirically, this works well to limit the search area
    thresh[[length(thresh) + 1]] = 
        data.table(cost = cost, eq_n = "simulate", calc = "sum",
            u = uniroot(test_evolve, interval, cost = cost, simulate = TRUE, how = "sum", mode = "public")$root);
}

for (cexp in c(seq(-4.5, -1.8, 0.1), -1.69897000433602)) {
    cat("\n", cexp, "\n")
    cost = 5 * 10^cexp;
    interval = c(0.5 * sqrt(cost), sqrt(cost)); # Empirically, this works well to limit the search area
    thresh[[length(thresh) + 1]] = 
        data.table(cost = cost, eq_n = "simulate", calc = "sum", 
            u = uniroot(test_evolve, interval, cost = cost, simulate = TRUE, how = "sum", mode = "public")$root);
}

for (cexp in seq(-6.0, -4.6, 0.1)) {
    cat("\n", cexp, "\n")
    cost = 5 * 10^cexp;
    interval = c(0.5 * sqrt(cost), 0.8 * sqrt(cost)); # Empirically, this works well to limit the search area
    thresh[[length(thresh) + 1]] = 
        data.table(cost = cost, eq_n = "simulate", calc = "product",
            u = uniroot(test_evolve, interval, cost = cost, simulate = TRUE, how = "product", mode = "public")$root);
}

for (cexp in seq(-4.5, -2.4, 0.1)) {
    cat("\n", cexp, "\n")
    cost = 5 * 10^cexp;
    interval = c(0.5 * sqrt(cost), sqrt(cost)); # Empirically, this works well to limit the search area
    thresh[[length(thresh) + 1]] = 
        data.table(cost = cost, eq_n = "simulate", calc = "product", 
            u = uniroot(test_evolve, interval, cost = cost, simulate = TRUE, how = "product", mode = "public")$root);
}

for (cexp in c(seq(-2.3, -1.8, 0.1), -1.69897000433602)) {
    cat("\n", cexp, "\n")
    cost = 5 * 10^cexp;
    interval = c(0.5 * sqrt(cost), 1.5 * sqrt(cost)); # Empirically, this works well to limit the search area
    thresh[[length(thresh) + 1]] = 
        data.table(cost = cost, eq_n = "simulate", calc = "product", 
            u = uniroot(test_evolve, interval, cost = cost, simulate = TRUE, how = "product", mode = "public")$root);
}

for (cexp in c(seq(-6, -1.8, 0.1), -1.69897000433602)) {
    cat("\n", cexp, "\n")
    cost = 5 * 10^cexp;
    interval = c(0.4 * sqrt(cost), 2 * sqrt(cost));
    thresh[[length(thresh) + 1]] = 
        data.table(cost = cost, eq_n = "analytic", calc = "simple",
            u = uniroot(test_evolve, interval, cost = cost, simulate = FALSE, how = "simple", mode = "public")$root);
}
thresh = rbindlist(thresh);

# fwrite(thresh, "./public_thresh.csv")



#######################
# Simulation approach #
#######################

# In this part we find invasion thresholds by running invasion simulations.

# Tester for uniroot
test_invade = function(u, cost, type)
{
    cat(".")
    ext_invade(u, cost, type, gens = 100, nmax = 300, nbmax = 15, epsilon = 1e-12, verbose = FALSE)
}

rna_thresh = list()
for (cexp in c(-1.69897000433602, seq(-1.8, -3.4, -0.1))) {
    cost = 5 * 10^cexp;
    print(cost);
    interval = c(0.3, 0.6);
    rna_thresh[[length(rna_thresh) + 1]] = 
        data.table(cost = cost, type = "RNA", 
            u = uniroot(test_invade, interval, cost = cost, type = "RNA")$root);
    cat("\n");
}
rna_thresh = rbindlist(rna_thresh);
# fwrite(rna_thresh, "./rna_thresh.csv")
# rna_thresh = fread("./rna_thresh.csv")

dna_thresh = list()
for (cexp in c(-1.69897000433602, seq(-1.8, -5.0, -0.1))) {
    cost = 5 * 10^cexp;
    print(cost);
    # Empirically this seems to bound the threshold
    if (cexp >= -4) {
        interval = c((cost/2)^(1/2), (cost/0.3)^(1/1.6));
    } else {
        interval = c((cost/2)^(1/2), (cost/1.7)^(1/2));
    }
    print(interval);
    dna_thresh[[length(dna_thresh) + 1]] = 
        data.table(cost = cost, type = "DNA", 
            u = uniroot(test_invade, interval, cost = cost, type = "DNA")$root);
    cat("\n");
}
dna_thresh = rbindlist(dna_thresh);

# fwrite(dna_thresh, "./dna_thresh.csv")
# dna_thresh = fread("./dna_thresh.csv")




