library(data.table)
library(rlang)
library(ggplot2)
library(scales)
library(cowplot)
library(cmocean)

theme_set(theme_cowplot(font_size = 8) + theme(plot.title = element_text(size = 8, face = "plain"), strip.background = element_blank()))

# Empirically measured copy number & transposition rates
nuzhdin_mackay = fread(
"family, n_het, n_trans
297,	28,	0
mdg3,	16,	0
Doc,	28,	6
roo,	66,	103
copia,	29,	197")

maside_et_al = fread(
"family,        X,      C2Bal,  C2MA,   C3,     C4,     n_lines, n_trans
mdg-1,          4.0,    7.0,    6.0,    9.3,    0.0,    22,     3
opus,           2.0,    4.0,    7.0,    11.9,   0.0,    21,     1
copia,          4.0,    9.0,    12.0,   14.8,   0.0,    21,     10
1731,           0.0,    0.0,    0.0,    1.0,    0.0,    20,     0
297,            3.5,    16.0,   8.0,    19.1,   1.0,    21,     4
412,            1.0,    12.0,   12.0,   8.8,    0.0,    23,     0
jockey,         4.0,    8.0,    24.0,   16.2,   1.0,    21,     1
2244,           2.8,    9.0,    5.0,    12.2,   1.0,    20,     1
pogo,           2.8,    3.0,    5.0,    6.1,    0.0,    20,     0
roo,            16.7,   19.0,   24.0,   30.2,   0.0,    20,     7
17.6,           1.6,    3.0,    1.0,    8.5,    1.0,    21,     0" )

# Below, n_trans is the number of observed transpositions; n_het is
# the number of sites in the inbred line (e.g. starting heterozygous copy number);
# 31 is the number of lines, and 163 is the number of generations.
nuzhdin_mackay[, u_central := n_trans / (n_het * 31 * 163)]
nuzhdin_mackay[, u_lo := qgamma(0.025, n_trans + 0, rate = n_het * 31 * 163)]
nuzhdin_mackay[, u_hi := qgamma(0.975, n_trans + 1, rate = n_het * 31 * 163)]
nuzhdin_mackay[, study := "Nuzhdin and Mackay 1994"]

# Below, n_het is calculated wrt a male
maside_et_al[, n_het := 0.5 * X + 0.5 * C2Bal + 0.5 * C2MA + C3 + C4]
# 116 is the number of generations and 0.373 is the proportion of euchromatin in
# chromosome 2 (as insertions in chromosome 2 were counted)
maside_et_al[, u_central := n_trans / (n_het * n_lines * 116 * 0.373)]
maside_et_al[, u_lo := qgamma(0.025, n_trans + 0, rate = n_het * n_lines * 116 * 0.373)]
maside_et_al[, u_hi := qgamma(0.975, n_trans + 1, rate = n_het * n_lines * 116 * 0.373)]
maside_et_al[, study := "Maside et al. 2000"]

empirical = rbind(
    nuzhdin_mackay[, .(family, study, u_central, u_lo, u_hi, n_het)],
    maside_et_al  [, .(family, study, u_central, u_lo, u_hi, n_het)]
)

empirical[u_central != 0, type := "full"]
empirical[u_central == 0, type := "upper"]

merge_all = function(dts) {
    Reduce(function(...) merge(..., all = TRUE, by = c("form", "u_minus_v")), dts)
}

# Log fitness functions must have n, a as first arguments, all other parameters as subsequent arguments, and
# contain exactly one expression calculating ln(w(n)), which is not wrapped in a return() statement.

lnw_ect = function(n, a, lambda, epsilon, rho, p_g) (0.5 * (n * (n - 1) - n * p_g)) * log(1 - lambda * epsilon * rho/100)
lnw_ins = function(n, a, H, s, t) n * log1p(-H * s) + 0.5 * n * (n - 1) * log1p(-H^2 * t * s) # NOTE does not account for homozygous insertions
lnw_dsy = function(n, a, Y, lambda) log(1 - Y * lambda * n)
lnw_gsy = function(n, a, k, c) log(1 - k * c * a * n)
lnw_mis = function(n, a, M, k) log(1 - M * (1.4e-7 * k * a * n)^2)
lnw_dsb = function(n, a, B, C, j) -B * (j * a * n) - C * (j * a * n)^2
lnw_com = function(n, a, lambda, epsilon, rho, p_g, H, s, t, Y, k, c, M) 0 +
    (0.5 * (n * (n - 1) - n * p_g)) * log(1 - lambda * epsilon * rho/100) +
    n * log1p(-H * s) + 0.5 * n * (n - 1) * log1p(-H^2 * t * s) +
    log(1 - Y * lambda * n) +
    log(1 - k * c * a * n) +
    log(1 - M * (1.4e-7 * k * a * n)^2)
lnw_thL = function(n, a, q) -q * (0.005*n)^2
lnw_thH = function(n, a, q) -1000 * q * (0.005*n)^2
lnw_thA = function(n, a, q) -q * (a * n)^2

# These functions give the maximum copy number n for which the fitness function has a nonnegative value
n_limit_ect = function(a, lambda, epsilon, rho, p_g) Inf
n_limit_ins = function(a, H, s, t) Inf # NOTE does not account for homozygous insertions
n_limit_dsy = function(a, Y, lambda) 1 / (Y * lambda)
n_limit_gsy = function(a, k, c) 1 / (k * c * a)
n_limit_mis = function(a, M, k) sqrt(1 / (M * (1.4e-7 * k * a)^2)) - 0.1 # the - 0.1 is for numerical reasons
n_limit_dsb = function(a, B, C, j) Inf
n_limit_com = function(a, lambda, epsilon, rho, p_g, H, s, t, Y, k, c, M) 10^6
n_limit_thL = function(a, q) Inf
n_limit_thH = function(a, q) Inf
n_limit_thA = function(a, q) Inf

# Returns a function with the same signature as lnw_func, but which returns the derivative of lnw_func with respect to its first argument n
dlnw_dn = function(lnw_func)
{
    # Extract ln fitness expression from w_func
    lnw_expression = tail(as.list(lnw_func), 1)[[1]]
    
    # Construct derivative by n function
    deriv_func = deriv(lnw_expression, "n", function.arg = names(formals(lnw_func)))
    
    # Construct and return gradient function
    gradient_args = paste0(names(formals(lnw_func)), collapse = ", ");
    func_str = paste0("function(", gradient_args, ") unname(attr(deriv_func(", gradient_args, "), 'gradient')[1, 'n'])")
    eval(str2expression(func_str))
}

# Returns a function for use in uniroot
root_func = function(lnw_func, w_factor)
{
    # Get d ln(w) / dn
    dlnw_func = dlnw_dn(lnw_func)
    
    # Construct and return root function
    gradient_args = paste0(names(formals(lnw_func)), collapse = ", ");
    func_str = paste0("function(", gradient_args, ") a + (1 + 8 * a / 2 - 3 * a^2 / 4) * ", w_factor, " * dlnw_func(", gradient_args, ")")
    eval(str2expression(func_str))
}

lnw_func = lnw_thA
form_str = "Theor A"
umv_min = 1e-6
umv_max = 1
n_pts = 101
n_min = 1e-3
n_max = 1e+13
n_limit_func = n_limit_thA
fitness_args = list(q = 0.1)

# Get equilibrium "n" curve for a given set of parameter assumptions and a range of u_minus_v
eq_n_curve = function(lnw_func, form_str, umv_min, umv_max, n_pts, n_min, n_max, n_limit_func, fitness_args)
{
    # Bounds checking
    if (umv_min <= 0 || umv_max <= 0 || umv_min >= umv_max) {
        stop("Must have 0 < umv_min < umv_max.")
    }
    
    # Get umv points
    umv_pts_master = exp(seq(log(umv_min), log(umv_max), length.out = n_pts))
    
    # Set of w factors to consider
    w_factors = c(0.1, 1, 10)
    all_results = list()
    
    for (wf in w_factors) {
        # Get root func
        f = root_func(lnw_func, wf)
        
        # Find value of u_minus_v for which equilibrium n crosses zero; this is new min/max of u_minus_v
        f2 = function(u_minus_v) {
            do.call(f, c(list(n = 0, a = u_minus_v), fitness_args))
        }
        f2_umv_min = f2(umv_min);
        f2_umv_max = f2(umv_max);
        first_n_zero = FALSE;
        last_n_zero = FALSE;
        if (f2_umv_min < 0 && f2_umv_max < 0) {
            stop("No positive equilibrium n in [umv_min, umv_max].")
        } else if (f2_umv_min < 0 && f2_umv_max > 0) {
            umv_min = uniroot(f2, c(umv_min, umv_max))$root;
            first_n_zero = TRUE;
        } else if (f2_umv_min > 0 && f2_umv_max < 0) {
            umv_max = uniroot(f2, c(umv_min, umv_max))$root;
            last_n_zero = TRUE;
        }
        
        # Get umv points
        umv_pts = umv_pts_master[umv_pts_master %between% c(umv_min, umv_max)]
        
        # Find root
        results = data.table(u_minus_v = umv_pts, n = NA_real_, w = 1)
        if (first_n_zero) { 
            results[1, n := 0];
        }
        if (last_n_zero)  { 
            results[.N, n := 0]; 
        }
        for (umv in results[is.na(n), u_minus_v]) {
            n_hard_max = do.call(n_limit_func, c(list(a = umv), fitness_args))
            eq_n = do.call(uniroot, c(list(f = f, interval = c(n_min, min(n_max, n_hard_max))), list(a = umv), unname(fitness_args)))$root
            results[u_minus_v == umv, n := eq_n]
            results[u_minus_v == umv, w := exp(do.call(lnw_func, c(list(n = eq_n, a = umv), fitness_args)))]
        }
        
        results[, form := form_str]
        names(results)[2] = paste0("n", wf)
        names(results)[3] = paste0("w", wf)
        
        all_results[[as.character(wf)]] = results
    }
    
    return (merge_all(all_results))
}

###################
# EQ. N COST PLOT #
###################

# Assemble solutions
dt = rbind(
    eq_n_curve(lnw_ect, "Ectopic recombination", 1e-6, 1, 101, 1e-3, 1e+13, n_limit_ect, 
        list(lambda = 0.005, epsilon = 0.15, rho = 1, p_g = 0.0)),
    
    eq_n_curve(lnw_ins, "Deleterious insertions", 1e-6, 1, 1001, 1e-3, 1e+13, n_limit_ins, 
        list(H = 0.1, s = 0.01, t = 0.025)),

    eq_n_curve(lnw_dsy, "DNA synthesis", 1e-6, 1, 101, 1e-3, 1e+13, n_limit_dsy, 
        list(Y = 1e-5, lambda = 0.005)),
    
    eq_n_curve(lnw_gsy, "Transcription and translation", 1e-6, 1, 101, 1e-3, 1e+13, n_limit_gsy, 
        list(k = 10, c = 6e-8)),

    eq_n_curve(lnw_mis, "Protein aggregation", 1e-6, 1, 101, 1e-3, 1e+13, n_limit_mis, 
        list(M = 7.1865, k = 10)),
    
    eq_n_curve(lnw_dsb, "Double-stranded breaks", 1e-6, 1, 101, 1e-3, 1e+13, n_limit_dsb,
        list(B = 2.32074e-3, C = 3.11734e-5, j = 10)),

    eq_n_curve(lnw_com, "Combined (RNA)", 1e-6, 1, 1001, 1e-3, 1e+13, n_limit_com,
        list(lambda = 0.005, epsilon = 0.15, rho = 1, p_g = 0.0,
            H = 0.1, s = 0.01, t = 0.025, Y = 1e-5, 
            k = 10, c = 6e-8, M = 7.1865))
    # ,
    # eq_n_curve(lnw_thH, "Theor H", 1e-3, 1, 101, 1e-3, 1e+13, n_limit_thH,
    #     list(q = 0.2))
    # ,
    # eq_n_curve(lnw_thA, "Theor A", 1e-6, 1, 101, 1e-3, 1e+13, n_limit_thA,
    #     list(q = 0.1))
)


dt[is.na(n1), n1 := 0]
dt[is.na(n10), n10 := 0]

dt[n0.1 == 0, w0.1 := 1]
dt[n1 == 0, w1 := 1]
dt[n10 == 0, w10 := 1]


dt[, form := factor(form, levels = c(
    "Transcription and translation",
    "Protein aggregation",
    "DNA synthesis",
    "Double-stranded breaks",
    "Deleterious insertions",
    "Ectopic recombination",
    "Combined (RNA)"))]

my_colours = c(
"Transcription and translation" =  "#6d2f20",
"Protein aggregation" =            "#b75347",
"DNA synthesis" =                  "#e09351",
"Double-stranded breaks" =         "#edc775",
"Deleterious insertions" =         "#94b594",
"Ectopic recombination" =          "#224b5e",
"Combined (RNA)" =                 "#000000"
)

my_linetypes = c(rep("solid", length(my_colours) - 1), rep("21", 1))
names(my_linetypes) = names(my_colours)

empirical[type == "upper", u_central := 1e-6]
empirical[type == "upper", u_lo := 1e-6]

cost_plot = ggplot(dt) +
    geom_line(aes(u_minus_v, n1, colour = form, linetype = form), size = 0.75) +
    geom_pointrange(data = empirical[type == "full"], 
        aes(xmin = u_lo, xmax = u_hi, x = u_central, y = n_het * 2), fatten = 1, size = 0.25, shape = 18) +
    # geom_pointrange(data = empirical[type == "upper"],
    #     aes(xmin = u_lo, xmax = u_hi, x = u_central, y = n_het * 2), fatten = 1, size = 0.25, shape = 18, colour = "#666666") +
    scale_x_log10(limits = c(7.8e-7, 1), breaks = 10^(-6:0), labels = trans_format("log10", math_format(10^.x)), expand = expansion(mult = c(0.01, 0.02))) +
    scale_y_log10(limits = c(6e-1, 1.7e12), breaks = 10^seq(0, 12, by = 2), labels = trans_format("log10", math_format(10^.x)), expand = expansion(mult = c(0.01, 0.05))) +
    scale_colour_manual(values = my_colours) +
    scale_linetype_manual(values = my_linetypes) +
    labs(x = "Duplication rate", y = "Equilibrium copy number", colour = "Fitness cost", linetype = "Fitness cost") +
    theme(legend.position = c(0.55, 0.87), legend.title = element_text(size = 7, margin = margin(0, 0, -3, 0)), 
        plot.margin = margin(2, 0, 0, 3),
        legend.text = element_text(size = 6), legend.key.height = unit(0.12, "cm"),
        axis.title.y = element_text(margin = margin(0, -2, 0, 0)))




##################
# DYNAMICS PLOTS #
##################

# Cost simulations
source("./load_ngd.R")
library(stringr)

load_set = function(file_pattern, cost_type, path = "./TE/Runs/2-Cost/", summarize = TRUE)
{
    fp = paste0("^", str_replace_all(file_pattern, "\\*", "([0-9]+[kmgKMG]?)"), "$")
    files = list.files(path = path, pattern = fp)
    measures = str_match(files, fp)[,2]
    measures = as.numeric(str_replace_all(measures, c("^0" = "0.", k = "000", K = "000", m = "000000", M = "000000", g = "000000000", G = "000000000")))
    files = paste0(path, files)
    
    data = list()
    for (f in seq_along(files)) {
        d = load_ngd(files[f], thin = 1, gunit = 1)
        if (summarize) {
            d = d[, .(extinction = tail(w,1) <= 1e-6, w = tail(w, 1), n = tail(n, 1), u = tail(u, 1)), by = run]
        }
        d[, cost := measures[f]]
        d[, cost_type := cost_type]
        data[[f]] = d
    }
    
    all_data = rbindlist(data)
    if (summarize) {
        all_data = all_data[order(cost, run)]
    }
    return (all_data)
}

############################
# THEORETICAL ILLUSTRATION #
############################

du = load_set("theorN_*.ngd", "Public replication", summarize = FALSE);
dl = load_set("theorN_*.ngd", "Private replication", summarize = FALSE);
dh = load_set("theorNHigh_*.ngd", "Private replication, increased cost", summarize = FALSE);
da = load_set("theorA_*.ngd", "Private replication, activity cost", summarize = FALSE);

plot_triad = function(dl, x_labels, y_labels, title, legend, nymin, nymax, nlabs)
{
    d = dl[order(cost)];
    d[, costf := factor(cost, levels = unique(cost))];
    d[, g := g / 1000];
    
    nbreaks = nlabs;
    nlabels = sapply(nlabs, function(x) sprintf("%g", x));
    oceanname = "haline";
    
    plN = ggplot(d) +
        geom_line(aes(x = g, y = n, group = interaction(run, costf), colour = costf), size = 0.25) +
        geom_vline(xintercept = 1, linetype = "dashed", size = 0.25) +
        scale_y_log10(limits = c(nymin, nymax), breaks = nbreaks, labels = nlabels) +
        scale_x_continuous(limits = c(0, 21), expand = expansion(c(0.04, 0))) +
        scale_colour_cmocean(discrete = TRUE, name = oceanname, end = 0.92) + 
        theme(legend.position = if (legend) c(0.5, 0.6) else "none", legend.title = element_text(size = 7, margin = margin(0, 0, -3, 0)),
            legend.text = element_text(size = 6), legend.key.height = unit(0.12, "cm"), plot.margin = margin(2, 0, 0, 1),
            axis.text = element_text(size = 6)) +
        labs(x = NULL, y = if (y_labels) "Copy number" else NULL, title = title, colour = expression("Cost," ~ italic(q)))

    plU = ggplot(d) +
        geom_line(aes(x = g, y = u, group = interaction(run, costf), colour = costf), size = 0.25) +
        geom_vline(xintercept = 1, linetype = "dashed", size = 0.25) +
        scale_y_log10(limits = c(0.001, 1), breaks = c(0.001, 0.01, 0.1, 1), labels = c("0.001", "0.01", "0.1", "1")) +
        scale_x_continuous(limits = c(0, 21), expand = expansion(c(0.04, 0))) +
        scale_colour_cmocean(discrete = TRUE, name = oceanname, end = 0.92) + 
        theme(legend.position = "none", plot.margin = margin(2, 0, 0, 1),
            axis.text = element_text(size = 6)) +
        labs(x = NULL, y = if (y_labels) "Duplication rate" else NULL)

    plW = ggplot(d) +
        geom_line(aes(x = g, y = w, group = interaction(run, costf), colour = costf), size = 0.25) +
        geom_vline(xintercept = 1, linetype = "dashed", size = 0.25) +
        scale_y_continuous(limits = c(0, 1)) +
        scale_x_continuous(limits = c(0, 21), expand = expansion(c(0.04, 0))) +
        scale_colour_cmocean(discrete = TRUE, name = oceanname, end = 0.92) + 
        theme(legend.position = "none", plot.margin = margin(2, 0, 0, 1),
            axis.text = element_text(size = 6)) +
        labs(x = if (x_labels) "Generations (thousands)" else "", y = if (y_labels) "Host fitness" else NULL)
    
    plot_grid(plN, plU, plW, ncol = 1, align = "v", axis = "b", rel_heights = c(1.25, 0.7, 0.9))
}

pL = plot_triad(dl, FALSE, TRUE,  "Moderate copy-\nnumber cost",  TRUE,  0.03, 10000, c(0.1, 1, 10, 100, 1000, 10000))
pH = plot_triad(dh, TRUE, FALSE,  "High copy-\nnumber cost", FALSE, 0.03, 10000, c(0.1, 1, 10, 100, 1000, 10000))
pE = plot_triad(da, FALSE, FALSE, "Moderate\nactivity cost",      FALSE, 0.03, 10000, c(0.1, 1, 10, 100, 1000, 10000))

#######################
# FIG 2 FOR MAIN TEXT #
#######################

t1 = load_ngd("./TE/Runs/2-Cost/ect_1cMMb.ngd", thin = 1)[run == 1]
t2 = load_ngd("./TE/Runs/2-Cost/ins_s1.ngd", thin = 1)[run == 1]
t3 = load_ngd("./TE/Runs/2-Cost/syn_5kb.ngd", thin = 1)[run == 1]
t4 = load_ngd("./TE/Runs/2-Cost/dsb_10.ngd", thin = 1)[run == 1]
t5 = load_ngd("./TE/Runs/2-Cost/exp_k10.ngd", thin = 1)[run == 1]
t6 = load_ngd("./TE/Runs/2-Cost/mis_k10.ngd", thin = 1)[run == 1]

fc1 = plot_grid(
    panel(t1, "n", "n05", "n95", c(0, 10120), cTE1, title = "Ectopic recombination", tx = TRUE, xlim = c(0, 3), ylab = "Copy\nnumber") + vline(1),
    panel(t1, "u", "u05", "u95", c(0, 0.031), cDup, tx = TRUE, xlim = c(0, 3), ylab = "Duplication\nrate", xlab = " ") + vline(1),
    ncol = 1, align = "v", axis = "b", rel_heights = c(1.4, 1, 1))

fc2 = plot_grid(
    panel(t2, "n", "n05", "n95", c(0, 10120), cTE1, title = "Deleterious insertions", tx = TRUE, xlim = c(0, 3)) + vline(1),
    panel(t2, "u", "u05", "u95", c(0, 0.031), cDup, tx = TRUE, xlim = c(0, 3), xlab = "Generations (thousands)") + vline(1),
    ncol = 1, align = "v", axis = "b", rel_heights = c(1.4, 1, 1))

fc3 = plot_grid(
    panel(t3, "n", "n05", "n95", c(0, 10120), cTE1, title = "DNA synthesis", tx = TRUE, xlim = c(0, 3)) + vline(1),
    panel(t3, "u", "u05", "u95", c(0, 0.031), cDup, tx = TRUE, xlim = c(0, 3), xlab = " ") + vline(1),
    ncol = 1, align = "v", axis = "b", rel_heights = c(1.4, 1, 1))

fc4 = plot_grid(
    panel(t4, "n", "n05", "n95", c(0, 10120), cTE1, title = "Double-stranded breaks", tx = TRUE, xlim = c(0, 3), ylab = "Copy\nnumber") + vline(1),
    panel(t4, "u", "u05", "u95", c(0, 0.031), cDup, tx = TRUE, xlim = c(0, 3), ylab = "Duplication\nrate", xlab = " ") + vline(1),
    ncol = 1, align = "v", axis = "b", rel_heights = c(1.4, 1, 1))

fc5 = plot_grid(
    panel(t5, "n", "n05", "n95", c(0, 10120), cTE1, title = "Transcription and translation", tx = TRUE, xlim = c(0, 3)) + vline(1),
    panel(t5, "u", "u05", "u95", c(0, 0.031), cDup, tx = TRUE, xlim = c(0, 3), xlab = "Generations (thousands)") + vline(1),
    ncol = 1, align = "v", axis = "b", rel_heights = c(1.4, 1, 1))

fc6 = plot_grid(
    panel(t6, "n", "n05", "n95", c(0, 10120), cTE1, title = "Protein aggregation", tx = TRUE, xlim = c(0, 3)) + vline(1),
    panel(t6, "u", "u05", "u95", c(0, 0.031), cDup, tx = TRUE, xlim = c(0, 3), xlab = " ") + vline(1),
    ncol = 1, align = "v", axis = "b", rel_heights = c(1.4, 1, 1))

fcc = plot_grid(fc1, fc2, fc3, fc4, fc5, fc6, nrow = 2,
    labels = "C", label_size = 8, rel_widths = c(1.2, 1, 1))

costs = rbind(
    load_set("exp_k*.ngd", "Transcription and translation"),
    load_set("mis_k*.ngd", "Protein aggregation"),
    load_set("syn_*b.ngd", "DNA synthesis"),
    load_set("dsb_*.ngd",  "Double-stranded breaks"),
    load_set("ins_s*.ngd", "Deleterious insertions"),
    load_set("ect_*cMMb.ngd", "Ectopic recombination"),
    load_set("com_*.ngd", "Combined (RNA)")
)

costs[n >= 10000, extinction := TRUE]
costs[extinction == TRUE, w := 0]

costs[cost_type == "Transcription and translation", standard_cost := 10] # k = 10
costs[cost_type == "Protein aggregation",           standard_cost := 10] # k = 10
costs[cost_type == "DNA synthesis",                 standard_cost := 5000] # 5 kb
costs[cost_type == "Double-stranded breaks",        standard_cost := 10] # j = 10
costs[cost_type == "Deleterious insertions",        standard_cost := 1] # s = 0.01 (1%)
costs[cost_type == "Ectopic recombination",         standard_cost := 1] # 1 cM/Mb
costs[cost_type == "Combined (RNA)",                standard_cost := 1] # 1 unit
costs[, rel_cost := cost / standard_cost]

csumm = costs[
    cost_type != "Double-stranded breaks",
    # cost_type != "Combined (RNA)", 
    .(w = mean(w), n = mean(n), u = mean(u)), keyby = .(cost_type, standard_cost, rel_cost)]

csumm[, cost_type := factor(cost_type, levels = c(
    "Ectopic recombination",
    "Deleterious insertions",
    "DNA synthesis",
    "Transcription and translation",
    "Protein aggregation",
    "Combined (RNA)"
    ))]

csumm = csumm[rel_cost >= 1]
lastp = csumm[rel_cost <= 10^7, .(rel_cost = max(rel_cost)), keyby = cost_type]

w_plot = ggplot(csumm[rel_cost >= 1]) +
    geom_point(aes(x = rel_cost, y = w, colour = cost_type), size = 0.6) +
    geom_line(aes(x = rel_cost, y = w, colour = cost_type, linetype = cost_type), size = 0.6) +
    geom_segment(data = lastp, aes(x = rel_cost, xend = 10^7, y = 1, yend = 1, colour = cost_type, linetype = cost_type), size = 0.6) +
    facet_wrap(~cost_type, ncol = 1) +
    scale_x_log10(limits = c(1, 1e7), breaks = 10^(0:7), labels = trans_format("log10", math_format(10^.x)), expand = expansion(c(0.01, 0.02))) +
    scale_y_continuous(expand = expansion(0.03), limits = c(0, 1), breaks = (0:4)/4, labels = c("0", "", ".5", "", "1")) +
    scale_colour_manual(values = my_colours) +
    scale_linetype_manual(values = my_linetypes) +
    labs(x = "Relative cost increase", y = "Equilibrium host fitness") +
    theme(legend.position = "none", panel.background = element_rect(fill = "#f4f4f4"), panel.grid.major = element_line(colour = "white"))

f2 = plot_grid(fcc, w_plot, nrow = 1, rel_widths = c(15, 5), 
    labels = c("C", "D"), label_size = 8)

ggsave("./Figures/2-costs-plots.pdf", f2, width = 20, height = 9, units = "cm", useDingbats = FALSE)
ggsave("./Figures/2-costs-plots.png", f2, width = 20, height = 9, units = "cm")


fSC = plot_grid(pL, pH, pE, cost_plot, nrow = 1, rel_widths = c(3.25, 2.85, 2.93, 7.5),
    labels = LETTERS[1:4], label_size = 8, label_x = c(0.15, 0.07, 0.06, 0), label_y = 1.004)

ggsave("./Figures/S1-costs.pdf", fSC, width = 16.5, height = 7.5, units = "cm", useDingbats = FALSE)
ggsave("./Figures/S1-costs.png", fSC, width = 16.5, height = 7.5, units = "cm")

################
# FOR SUPP MAT #
################

my_colours = c(
"Transcription and translation" =  "#6d2f20",
"Protein aggregation" =            "#b75347",
"DNA synthesis" =                  "#e09351",
"Double-stranded breaks" =         "#edc775",
"Deleterious insertions, 1-site" = "#94b594",
"Deleterious insertions, 2-site" = "#94b594",
"Ectopic recombination" =          "#224b5e",
"Combined (RNA), 1-site" =         "#000000",
"Combined (RNA), 2-site" =         "#000000"
)

my_linetypes = c(rep("solid", length(my_colours) - 3), rep("21", 3))
names(my_linetypes) = names(my_colours)

costs = rbind(
    load_set("exp_k*.ngd", "Transcription and translation"),
    load_set("mis_k*.ngd", "Protein aggregation"),
    load_set("syn_*b.ngd", "DNA synthesis"),
    load_set("dsb_*.ngd",  "Double-stranded breaks"),
    load_set("ins_s*.ngd", "Deleterious insertions, 1-site"),
    load_set("ins2_s*.ngd", "Deleterious insertions, 2-site"),
    load_set("ect_*cMMb.ngd", "Ectopic recombination"),
    load_set("com_*.ngd", "Combined (RNA), 1-site"),
    load_set("com2_*.ngd", "Combined (RNA), 2-site")
)

costs[n >= 10000, extinction := TRUE]
costs[extinction == TRUE, w := 0]

costs[cost_type == "Transcription and translation",  standard_cost := 10] # k = 10
costs[cost_type == "Protein aggregation",            standard_cost := 10] # k = 10
costs[cost_type == "DNA synthesis",                  standard_cost := 5000] # 5 kb
costs[cost_type == "Double-stranded breaks",         standard_cost := 10] # j = 10
costs[cost_type == "Deleterious insertions, 1-site", standard_cost := 1] # s = 0.01 (1%)
costs[cost_type == "Deleterious insertions, 2-site", standard_cost := 1] # s = 0.01 (1%)
costs[cost_type == "Ectopic recombination",          standard_cost := 1] # 1 cM/Mb
costs[cost_type == "Combined (RNA), 1-site",         standard_cost := 1] # 1 unit
costs[cost_type == "Combined (RNA), 2-site",         standard_cost := 1] # 1 unit
costs[, rel_cost := cost / standard_cost]


csumm = costs[,
    .(w = mean(w), n = mean(n), u = mean(u)), keyby = .(cost_type, standard_cost, rel_cost)]

csumm[, cost_type := factor(cost_type, levels = c(
    "Transcription and translation",
    "Protein aggregation",
    "DNA synthesis",
    "Double-stranded breaks",
    "Deleterious insertions, 1-site",
    "Deleterious insertions, 2-site",
    "Ectopic recombination",
    "Combined (RNA), 1-site",
    "Combined (RNA), 2-site"
    ))]

csumm = csumm[rel_cost >= 1]
lastp = csumm[rel_cost <= 10^7, .(rel_cost = max(rel_cost)), keyby = cost_type]

w_plot = ggplot(csumm[rel_cost >= 1]) +
    geom_point(aes(x = rel_cost, y = w, colour = cost_type), size = 0.6) +
    geom_line(aes(x = rel_cost, y = w, colour = cost_type, linetype = cost_type), size = 0.6) +
    geom_segment(data = lastp, aes(x = rel_cost, xend = 10^7, y = 1, yend = 1, colour = cost_type, linetype = cost_type), size = 0.6) +
    facet_wrap(~cost_type, ncol = 1) +
    scale_x_log10(limits = c(1, 1e7), breaks = 10^(0:7), labels = trans_format("log10", math_format(10^.x)), expand = expansion(0.01)) +
    scale_y_continuous(expand = expansion(0.03)) +
    scale_colour_manual(values = my_colours) +
    scale_linetype_manual(values = my_linetypes) +
    labs(x = "Relative cost increase", y = "Equilibrium host fitness") +
    theme(legend.position = "none", panel.background = element_rect(fill = "#f4f4f4"), panel.grid.major = element_line(colour = "white"),
        panel.spacing.y = unit(0, "cm"), axis.text = element_text(size = 6), strip.text = element_text(size = 6))

n_plot = ggplot(csumm[rel_cost >= 1 & w > 0 & w < 1]) +
    geom_point(aes(x = rel_cost, y = n, colour = cost_type), size = 0.6) +
    geom_line(aes(x = rel_cost, y = n, colour = cost_type, linetype = cost_type), size = 0.6) +
    facet_wrap(~cost_type, ncol = 1, drop = FALSE) +
    scale_x_log10(limits = c(1, 1e7), breaks = 10^(0:7), labels = trans_format("log10", math_format(10^.x)), expand = expansion(0.01)) +
    scale_y_log10() +
    scale_colour_manual(values = my_colours) +
    scale_linetype_manual(values = my_linetypes) +
    labs(x = "Relative cost increase", y = "Equilibrium copy number") +
    theme(legend.position = "none", panel.background = element_rect(fill = "#f4f4f4"), panel.grid.major = element_line(colour = "white"),
        panel.spacing.y = unit(0, "cm"), axis.text = element_text(size = 6), strip.text = element_text(size = 6))

u_plot = ggplot(csumm[rel_cost >= 1 & w > 0 & w < 1]) +
    geom_point(aes(x = rel_cost, y = u, colour = cost_type), size = 0.6) +
    geom_line(aes(x = rel_cost, y = u, colour = cost_type, linetype = cost_type), size = 0.6) +
    facet_wrap(~cost_type, ncol = 1, drop = FALSE) +
    scale_x_log10(limits = c(1, 1e7), breaks = 10^(0:7), labels = trans_format("log10", math_format(10^.x)), expand = expansion(0.01)) +
    scale_y_log10(expand = expansion(0.03)) +
    scale_colour_manual(values = my_colours) +
    scale_linetype_manual(values = my_linetypes) +
    labs(x = "Relative cost increase", y = "Equilibrium duplication rate") +
    theme(legend.position = "none", panel.background = element_rect(fill = "#f4f4f4"), panel.grid.major = element_line(colour = "white"),
        panel.spacing.y = unit(0, "cm"), axis.text = element_text(size = 6), strip.text = element_text(size = 6))

fSC = plot_grid(w_plot, n_plot, u_plot, nrow = 1, rel_widths = c(4, 4, 4), labels = LETTERS[1:3], label_size = 8)

ggsave("./Figures/supp-costs-full.pdf", fSC, width = 16, height = 14, units = "cm", useDingbats = FALSE)
ggsave("./Figures/supp-costs-full.png", fSC, width = 16, height = 14, units = "cm")
