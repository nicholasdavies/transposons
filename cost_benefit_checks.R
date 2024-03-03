# Various checks on cost-benefit calculations

library(scales)


# Load simulation results
source("./cost_benefit_load.R")

# Total variation distance
tvdist = function(p, q)
{
    0.5 * sum(abs(p-q))
}

check_dist = function(cost, u, plot = FALSE, legend = TRUE, selection = TRUE)
{
    num_cost = as.numeric(paste0("0.", cost));
    num_u    = as.numeric(paste0("0.", u));
    
    # If selection, then look at simulation-generated genotypes
    genotypes = fread(paste0("./Deterministic/Runs/6-Genotypes/pri_c", cost, "_u", u, "_old.ngd"), skip = "g\t");

    genotypes = genotypes[g == GC];
    Nw = genotypes[, max(A) + 1];
    Nm = genotypes[, max(B) + 1];
    f = matrix(0, Nw, Nm);
    f[cbind(genotypes$A + 1, genotypes$B + 1)] = genotypes$f;
    
    margin_wt = rowSums(f);
    margin_mu = colSums(f)[-1] / sum(colSums(f)[-1]);
    
    if (!selection) {
        margin_mu = eq_mut(20, 1000, num_u, 0.0);
        margin_mu = c(margin_mu, rep(0, Nm - length(margin_mu) - 1));
        cat(".")
    }

    # Check of wild-type distribution
    n = weighted.mean(0:(Nw - 1), margin_wt);
    cw = data.table(x = 0:(Nw - 1), pm = margin_wt, pt = dpois(0:(Nw - 1), n));
    
    # Check of mutant distribution
    cm = data.table(x = 1:(Nm - 1), pm = margin_mu, pt = dquasi1(1:(Nw - 1), num_u));
    
    if (plot) {
        margin_mu0 = eq_mut(20, 1000, num_u, 0.0);
        margin_mu0 = c(margin_mu0, rep(0, Nm - length(margin_mu0) - 1));
        cm[, pm0 := margin_mu0]
        
        dat = rbind(
            cw[pm/max(pt) > 5e-4 | pt/max(pt) > 5e-4][, 
            .(x = c(x, x), p = c(pm, pt), 
                a = c(rep("Exact", .N), rep("Approximation", .N)), s = "Wild-type")],
            cm[pm/max(pt) > 5e-4 | pt/max(pt) > 5e-4][, 
            .(x = c(x, x), p = c(pm, pt), 
                a = c(rep("Exact", .N), rep("Approximation", .N)), s = "Mutant, with selection")],
            cm[pm0/max(pt) > 5e-4 | pt/max(pt) > 5e-4][, 
            .(x = c(x, x), p = c(pm0, pt), 
                a = c(rep("Exact", .N), rep("Approximation", .N)), s = "Mutant, no selection")]
        )
        dat[, s := factor(s, unique(s))]
        
        blanks = dat[, .(max_p = max(p)), by = s]
        blanks[, max_p := max(max_p), by = stringr::str_sub(s, 1, 1)]
        
        ggplot(dat) +
            geom_col(aes(x, p, fill = a, group = a), 
                alpha = 0.5, position = position_dodge(width = 0.5), 
                colour = "black", linewidth = 0.2) +
            geom_blank(data = blanks, aes(0, max_p)) +
            facet_wrap(~s, nrow = 1, scales = "free") +
            scale_x_continuous(limits = c(min(dat$x) - 0.75, max(dat$x) + 0.75), 
                breaks = min(dat$x):max(dat$x), expand = c(0, 0)) +
            cowplot::theme_cowplot(font_size = 8) +
            theme(legend.position = if (legend) c(0.15, 0.8) else "none", 
                strip.background = element_blank(),
                strip.text = element_text(size = 8)) +
            labs(x = "Transposon number", y = "Frequency", fill = NULL)
    } else {
        data.table(cost = num_cost, u = num_u, w_tvd = cw[, tvdist(pm, pt)], m_tvd = cm[, tvdist(pm, pt)]);
    }
}

###########################################################
# How much do approximations differ from exact solutions  #
###########################################################

tvd = rbind(
    check_dist("0005", "001"),
    check_dist("0005", "005"),
    check_dist("0005", "010"),
    check_dist("0005", "015"),
    check_dist("0005", "020"),
    check_dist("005", "010"),
    check_dist("005", "050"),
    check_dist("005", "100"),
    check_dist("005", "150"),
    check_dist("005", "200"),
    check_dist("05", "100"),
    check_dist("05", "200"),
    check_dist("05", "300"),
    check_dist("05", "400"),
    check_dist("05", "500")
)

tvd = melt(tvd, measure.vars = c("w_tvd", "m_tvd"))
tvd[variable == "w_tvd", variable := "Wild-type"]
tvd[variable == "m_tvd", variable := "Mutant, with selection"]
tvd[, cost := as.character(cost)]
tvd[cost == "5e-04", cost := "0.0005"]

# Also do mutants without selection
tvdNS = rbind(
    check_dist("0005", "001", selection = FALSE),
    check_dist("0005", "005", selection = FALSE),
    check_dist("0005", "010", selection = FALSE),
    check_dist("0005", "015", selection = FALSE),
    check_dist("0005", "020", selection = FALSE),
    check_dist("005", "010", selection = FALSE),
    check_dist("005", "050", selection = FALSE),
    check_dist("005", "100", selection = FALSE),
    check_dist("005", "150", selection = FALSE),
    check_dist("005", "200", selection = FALSE),
    check_dist("05", "100", selection = FALSE),
    check_dist("05", "200", selection = FALSE),
    check_dist("05", "300", selection = FALSE),
    check_dist("05", "400", selection = FALSE),
    check_dist("05", "500", selection = FALSE)
)

tvdNS = melt(tvdNS, measure.vars = c("w_tvd", "m_tvd"))
tvdNS[variable == "w_tvd", variable := "Wild-type"]
tvdNS[variable == "m_tvd", variable := "Mutant, no selection"]
tvdNS[, cost := as.character(cost)]
tvdNS[cost == "5e-04", cost := "0.0005"]

tvd = rbind(tvd, tvdNS[variable == "Mutant, no selection"])

pC = ggplot(tvd) + 
    geom_line(aes(x = u, y = value, colour = cost, group = cost)) +
    geom_point(aes(x = u, y = value, colour = cost, group = cost, shape = cost)) +
    facet_wrap(~variable, nrow = 1, scales = "free_x") +
    scale_x_log10(limits = c(0.001, 0.5)) +
    cowplot::theme_cowplot(font_size = 8) +
    theme(legend.position = c(0.03, 0.7), strip.background = element_blank(),
        strip.text = element_text(size = 8)) +
    scale_colour_manual(values = c("0.0005" = "#ffa956", "0.005" = "#cc7344", "0.05" = "#993010")) +
    scale_shape_manual(values = c(2, 0, 20)) +
    labs(x = "Duplication rate", y = "Total variation distance\nfrom approximation", colour = "Cost", shape = "Cost")

check_dist_2d = function(cost, u)
{
    genotypes = fread(paste0("./Deterministic/Runs/6-Genotypes/pri_c", cost, "_u", u, "_old.ngd"), skip = "g\t");
    
    num_cost = as.numeric(paste0("0.", cost));
    num_u    = as.numeric(paste0("0.", u));

    genotypes = genotypes[g == GC];
    Nw = genotypes[, max(A) + 1];
    Nm = genotypes[, max(B) + 1];
    x = matrix(0, Nw, Nm);
    x[cbind(genotypes$A + 1, genotypes$B + 1)] = genotypes$f;
    
    # Check of joint distribution
    n = weighted.mean(0:(Nw - 1), rowSums(x));
    phi = 1 - sum(x[,1]) # fraction of hosts with any mutant transposons
    genotypes[, fa_w := dpois(A, n)];
    genotypes[, fa_m := ifelse(B == 0, 1 - phi, phi * dquasi1(B, num_u))];
    genotypes[, fa := fa_w * fa_m]
    
    C = cov.wt(cbind(genotypes$A, genotypes$B), genotypes$f, 
        method = "ML")

    data.table(cost = num_cost, u = num_u, j_tvd = genotypes[, tvdist(f, fa)],
        sAB = C$cov[1,2]/C$cov[2,2]);
}


tvd2 = rbind(
    check_dist_2d("0005", "001"),
    check_dist_2d("0005", "005"),
    check_dist_2d("0005", "010"),
    check_dist_2d("0005", "015"),
    check_dist_2d("0005", "020"),
    check_dist_2d("005", "010"),
    check_dist_2d("005", "050"),
    check_dist_2d("005", "100"),
    check_dist_2d("005", "150"),
    check_dist_2d("005", "200"),
    check_dist_2d("05", "100"),
    check_dist_2d("05", "200"),
    check_dist_2d("05", "300"),
    check_dist_2d("05", "400"),
    check_dist_2d("05", "500")
)

tvd2 = melt(tvd2, measure.vars = c("j_tvd", "sAB"))
tvd2[, cost := as.character(cost)]
tvd2[cost == "5e-04", cost := "0.0005"]

pD = ggplot(tvd2[variable == "j_tvd"]) + 
    geom_line(aes(x = u, y = value, colour = cost, group = cost)) +
    geom_point(aes(x = u, y = value, colour = cost, group = cost, shape = cost)) +
    scale_x_log10(limits = c(0.001, 0.5)) +
    cowplot::theme_cowplot(font_size = 8) +
    theme(legend.position = "none", strip.background = element_blank(),
        strip.text = element_text(size = 8), plot.title = element_text(size = 8, face = "plain", hjust = 0.5)) +
    scale_colour_manual(values = c("0.0005" = "#ffa956", "0.005" = "#cc7344", "0.05" = "#993010")) +
    scale_shape_manual(values = c(2, 0, 20)) +
    labs(x = "Duplication rate", y = "Total variation distance\nfrom marginal product", 
        title = "Joint distribution", colour = "Cost", shape = "Cost")

pE = ggplot(tvd2[variable == "sAB"]) + 
    geom_line(aes(x = u, y = value, colour = cost, group = cost)) +
    geom_point(aes(x = u, y = value, colour = cost, group = cost, shape = cost)) +
    scale_x_log10(limits = c(0.001, 0.5)) +
    cowplot::theme_cowplot(font_size = 8) +
    theme(legend.position = "none", strip.background = element_blank(),
        strip.text = element_text(size = 8), plot.title = element_text(size = 8, face = "plain", hjust = 0.5)) +
    scale_colour_manual(values = c("0.0005" = "#ffa956", "0.005" = "#cc7344", "0.05" = "#993010")) +
    scale_shape_manual(values = c(2, 0, 20)) +
    labs(x = "Duplication rate", y = "Slope of wild-type number\non mutant number", 
        title = "Joint distribution", colour = "Cost", shape = "Cost")


p = cowplot::plot_grid(
    check_dist("0005", "001", TRUE), 
    check_dist("05", "500", TRUE, legend = FALSE), 
    pC, 
    cowplot::plot_grid(pD, pE, nrow = 1, labels = c("D", "E"), label_size = 9),
    nrow = 4, labels = c("A", "B", "C", ""), label_size = 9,
    rel_heights = c(0.7, 0.7, 0.95, 1.05))

ggsave("./Figures/S-invasion_check_dist_2d.pdf", p, width = 15, height = 15, units = "cm", useDingbats = FALSE)




# How much of a difference does the marginalisation make in calculations?
# Answer below

# Turn a matrix into the product of its marginals
marginalize = function(x)
{
    rowSums(x) %o% colSums(x);
}

approximalize = function(x, u)
{
    nn = weighted.mean(0:99, rowSums(x));
    phi = 1-sum(x[,1]);
    
    mr = dpois(0:99, nn);
    mc = c(1 - phi, dquasi1(1:99, u) * phi);
    mr %o% mc;
}

approximalize_n = function(x, u)
{
    nn = weighted.mean(0:99, rowSums(x));
    
    mr = dpois(0:99, nn);
    mc = colSums(x);
    mr %o% mc;
}

approximalize_m = function(x, u)
{
    phi = 1-sum(x[,1]);
    
    mr = rowSums(x);
    mc = c(1 - phi, dquasi1(1:99, u) * phi);
    mr %o% mc;
}

check_N_calcs = function(Cost, Dupl)
{
    # Load genotype file
    str_cost = stringr::str_remove_all(sprintf("%g", Cost), "0\\.");
    str_u = stringr::str_remove_all(sprintf("%.03f", Dupl), "0\\.");
    genotypes = fread(paste0("./Deterministic/Runs/6-Genotypes/pri_c", str_cost, "_u", str_u, "_old.ngd"), skip = "g\t");
    
    # Get genotype matrix
    genotypes = genotypes[g == GC];
    Nw = genotypes[, max(A) + 1];
    Nm = genotypes[, max(B) + 1];
    f = matrix(0, Nw, Nm);
    f[cbind(genotypes$A + 1, genotypes$B + 1)] = genotypes$f;
    
    # Get calculated Ns
    N_full = Ns(f, Dupl, DD);
    N_marg = Ns(marginalize(f), Dupl, DD);
    N_appr = Ns(approximalize(f, Dupl), Dupl, DD);
    # N_appn = Ns(approximalize_n(f, Dupl), Dupl);
    # N_appm = Ns(approximalize_m(f, Dupl), Dupl);
    
    # Get Ns from simulation
    N0 = ds[cost == Cost & dupl == Dupl & mode == "pri" & stage == 2 & tag == "old" & g == GC - 1, N_B]
    DNPar1 = ds[cost == Cost & dupl == Dupl & mode == "pri" & stage == 2 & tag == "old" & g == GC, N_B] - N0
    DNOffPr1 = ds[cost == Cost & dupl == Dupl & mode == "pri" & stage == 2 & tag == "new" & g == GC, N_B] - N0
    DNOffPu1 = ds[cost == Cost & dupl == Dupl & mode == "pub" & stage == 2 & tag == "new" & g == GC, N_B] - N0
    
    cat(".")
    
    data.table(u = Dupl, cost = Cost,
        variable = factor(rep(c("N0", "DNPar1", "DNOffPr1", "DNOffPu1"), each = 4), levels = c("N0", "DNPar1", "DNOffPr1", "DNOffPu1")),
        method = factor(rep(c("Sim", "Full", "Marg", "Appx"), 4), levels = c("Sim", "Full", "Marg", "Appx")),
        value = c(N0, N_full[["N0"]], N_marg[["N0"]], N_appr[["N0"]],
            DNPar1,   N_full[["NPar1"]]   - N_full[["N0"]], N_marg[["NPar1"]]   - N_marg[["N0"]], N_appr[["NPar1"]]   - N_appr[["N0"]],
            DNOffPr1, N_full[["NOffPr1"]] - N_full[["N0"]], N_marg[["NOffPr1"]] - N_marg[["N0"]], N_appr[["NOffPr1"]] - N_appr[["N0"]],
            DNOffPu1, N_full[["NOffPu1"]] - N_full[["N0"]], N_marg[["NOffPu1"]] - N_marg[["N0"]], N_appr[["NOffPu1"]] - N_appr[["N0"]]))
}

N_data = rbind(
    check_N_calcs(.0005, 0.001),
    check_N_calcs(.0005, 0.005),
    check_N_calcs(.0005, 0.010),
    check_N_calcs(.0005, 0.015),
    check_N_calcs(.0005, 0.020),
    check_N_calcs(.005, 0.010),
    check_N_calcs(.005, 0.050),
    check_N_calcs(.005, 0.100),
    check_N_calcs(.005, 0.150),
    check_N_calcs(.005, 0.200),
    check_N_calcs(.05, 0.100),
    check_N_calcs(.05, 0.200),
    check_N_calcs(.05, 0.300),
    check_N_calcs(.05, 0.400),
    check_N_calcs(.05, 0.500)
)

# Save results as qs to retain factor ordering.
# qs::qsave(N_data, "./N_data.qs") # (with m from 1 to infinity; with u as argument to quasigeometric distribution)
# qs::qsave(N_data, "./N_data_m_1_2.qs") # with m from 1 to 2
# qs::qsave(N_data, "./N_data_Q_u(1-u).qs") # with u*(1-u) as argument to quasigeometric distribution

N_data = qs::qread("./N_data.qs")
N_data_m = qs::qread("./N_data_m_1_2.qs")
N_data_u = qs::qread("./N_data_Q_u(1-u).qs")

ggplot(N_data_u) + 
    geom_line(aes(x = u, y = value, group = method, colour = method, linetype = method)) +
    facet_grid(variable~cost, scales = "free")

# Assemble data
N_data_plot = rbind(
    N_data  [method %in% c("Sim", "Marg", "Appx") & !(variable == "DNOffPu1" & method == "Appx")],
    N_data_m[method %in% c("Appx") & variable == "DNOffPu1", .(u, cost, variable, method = "Appx", value)]
)

# Format labels
N_data_plot[, cost_label := as.character(cost)]
N_data_plot[cost == 5e-04, cost_label := "cost == 5 %*% 10^{-4}"]
N_data_plot[cost == 5e-03, cost_label := "cost == 0.005"]
N_data_plot[cost == 5e-02, cost_label := "cost == 0.05"]
N_data_plot[, cost_label := factor(cost_label, unique(cost_label))]

N_data_plot[variable == "N0",       variable_label := "italic(N)[0]"]
N_data_plot[variable == "DNPar1",   value := value / DD]
N_data_plot[variable == "DNPar1",   variable_label := "Delta*italic(N)[1]^{Par} / italic(d)"]
N_data_plot[variable == "DNOffPr1", variable_label := "Delta*italic(N)[list(Pr,1)]^{Off}"]
N_data_plot[variable == "DNOffPu1", variable_label := "Delta*italic(N)[list(Pu,1)]^{Off}"]
N_data_plot[, variable_label := factor(variable_label, c("italic(N)[0]", "Delta*italic(N)[list(Pr,1)]^{Off}", "Delta*italic(N)[list(Pu,1)]^{Off}", "Delta*italic(N)[1]^{Par} / italic(d)"))]

N_data_plot[method == "Sim", method_label := "Joint distribution"]
N_data_plot[method == "Marg", method_label := "Marginal product"]
N_data_plot[method == "Appx", method_label := "Approximation"]
N_data_plot[, method_label := factor(method_label, unique(method_label))]

ggplot(N_data_plot) + 
    geom_line(aes(x = u, y = value, group = method_label, colour = method_label, linetype = method_label, linewidth = method_label)) +
    facet_grid(variable_label ~ cost_label, scales = "free", labeller = label_parsed, switch = "y") +
    cowplot::theme_cowplot(font_size = 8) +
    theme(strip.placement = "outside", strip.background = element_blank(), legend.position = c(0.685, 0.93),
        strip.text = element_text(size = 8), panel.grid = element_line(colour = "#eeeeee", linetype = 1, linewidth = 0.2),
        legend.key.width = unit(0.5, "cm")) +
    labs(y = NULL, x = "Duplication rate,"~italic(u), colour = "Method", linetype = "Method", linewidth = "Method") +
    scale_colour_manual(values = c("#000000", "#ff0000", "#4488ff")) +
    scale_linetype_manual(values = c(1, 1, 4)) +
    scale_linewidth_manual(values = c(1.5, 0.75, 0.75))

ggsave("./Figures/S-Napprox.pdf", width = 16, height = 12, units = "cm", useDingbats = FALSE)
    
N_data


check_N_calcs(.05, 0.500)

# OK so there is a discrepancy with NOffPu1 for high u high cost.
    


# N0
N0 = ds[cost == Cost & dupl == Dupl & mode == Mode & stage == 2 & tag == "old" & g == GC - 1, N_B]

# NPar1
# Note. The small discrepancy here (highest u highest cost) does not appear to be due to large DD or large analysis seed
# as I tried making these both smaller and this had no effect on the discrepancy.
(ds[cost == Cost & dupl == Dupl & mode == Mode & stage == 2 & tag == "old" & g == GC, N_B] - N0)/DD

# NOffPr1
ds[cost == Cost & dupl == Dupl & mode == Mode & stage == 2 & tag == "new" & g == GC, N_B] - N0

# NOffPu1
ds[cost == Cost & dupl == Dupl & mode == "pub" & stage == 2 & tag == "new" & g == GC, N_B] - N0


n1f = calc36(x, u)
n1m = calc36(marginalize(x), u)
n1a = calc36(approximalize(x, u), u)
n1a_n = calc36(approximalize_n(x, u), u)
n1a_m = calc36(approximalize_m(x, u), u)

Ns(x, u)
Ns(marginalize(x), u)


# Check genotypes are indeed for stage 1
nn = ds[dupl == 0.5 & g == GC & tag == "full" & stage == 1 & mode == "pri", a]
ds[dupl == 0.5 & g == GC & tag == "full" & stage == 1 & mode == "pri", a]
weighted.mean(0:99, rowSums(x))

ds[dupl == 0.5 & g == GC & tag == "full" & stage == 1 & mode == "pri", b]
weighted.mean(0:99, colSums(x))


# Do N0 calc
calcN0 = function(x, u)
{
    N = 0;
    D = 0;
    for (n in 0:99) {
      for (m in 1:(99-n)) {
        for (i in 0:25) {
          for (j in 0:25) {
            N = N + x[n+1, m+1] * dpois(i, n*u) * dpois(j, m*u) * (m + j) * (n + m + i + j)
            D = D + x[n+1, m+1] * dpois(i, n*u) * dpois(j, m*u) * (m + j)
          }
        }
      }
      cat(".");
    }
    cat("\n");
    N/D
}

# Do calc [36]
# NB this is N_1^Par, the initial difference to the parent transposon(s).
calc36 = function(x, u)
{
    N = 0;
    D = 0;
    u = num_u;
    for (n in 0:99) {
      for (m in 1:(99-n)) {
        for (i in 0:25) {
          for (j in 0:25) {
            N = N + x[n+1, m+1] * m * dpois(i, n*u) * dpois(j, m*u) * (m + j) * (n + m + i + j + 1)
            D = D + x[n+1, m+1] * m * dpois(i, n*u) * dpois(j, m*u) * (m + j)
          }
        }
      }
      cat(".");
    }
    cat("\n");
    N/D
}

n1f = calc36(x, u)
n1m = calc36(marginalize(x), u)
n1a = calc36(approximalize(x, u), u)
n1a_n = calc36(approximalize_n(x, u), u)
n1a_m = calc36(approximalize_m(x, u), u)

n0f = calcN0(x, u)
n0m = calcN0(marginalize(x), u)
n0a = calcN0(approximalize(x, u), u)
n0a_n = calcN0(approximalize_n(x, u), u)
n0a_m = calcN0(approximalize_m(x, u), u)
(nn * (1 - u^2) + (1 + u)^2 + u * (1 - u) / (1 + u)) / (1-u)

n1f - n0f
n1m - n0m
n1a - n0a
n1a_n - n0a_n
n1a_m - n0a_m

# So, to answer the question,
# How much of a difference does the marginalisation make in calculations?
# Marginalizing (removing correlation) doesn't have a huge effect; the biggest 
# effect is approximating the mutant distribution. Approximating the wild-type 
# distribution has in fact no impact (only depends on E(n)) for this \Delta N.



# What is m[1]
# Reason for looking into this was just the idea that instead of a quasigeometric
# distribution with parameter u, maybe something like u(1-u) would be a better
# approximation. What the below shows is that it would be a better approximation,
# but for u(1-u) it goes from approx[1] > exact[1] to approx[1] < exact[1],
# which is not as clean as for u where at least we always have approx[1] > exact[1].
# I think rather than make too much of a hash of things, better to stick with
# the relatively cleaner quasigeometric with parameter u. 

m_dist = function(cost, u)
{
    genotypes = fread(paste0("./Deterministic/Runs/6-Genotypes/pri_c", cost, "_u", u, "_old.ngd"), skip = "g\t");
    genotypes = genotypes[g == GC];
    Nw = genotypes[, max(A) + 1];
    Nm = genotypes[, max(B) + 1];
    x = matrix(0, Nw, Nm);
    x[cbind(genotypes$A + 1, genotypes$B + 1)] = genotypes$f;

    # Check of mutant distribution
    mx = colSums(x)[-1];
    mx = mx / sum(mx);
    cm = data.table(x = 1:(Nm - 1), pm = mx);
    
    return (cm$pm)
}

c(1 - m_dist("0005", "001")[1], 0.001 - 0.001^2)
c(1 - m_dist("0005", "005")[1], 0.005 - 0.005^2)
c(1 - m_dist("0005", "010")[1], 0.010 - 0.010^2)
c(1 - m_dist("0005", "015")[1], 0.015 - 0.015^2)
c(1 - m_dist("0005", "020")[1], 0.020 - 0.020^2)
c(1 - m_dist("005", "010")[1], 0.010 - 0.010^2)
c(1 - m_dist("005", "050")[1], 0.050 - 0.050^2)
c(1 - m_dist("005", "100")[1], 0.100 - 0.100^2)
c(1 - m_dist("005", "150")[1], 0.150 - 0.150^2)
c(1 - m_dist("005", "200")[1], 0.200 - 0.200^2)
c(1 - m_dist("05", "100")[1], 0.100 - 0.100^2)
c(1 - m_dist("05", "200")[1], 0.200 - 0.200^2)
c(1 - m_dist("05", "300")[1], 0.300 - 0.300^2)
c(1 - m_dist("05", "400")[1], 0.400 - 0.400^2)
c(1 - m_dist("05", "500")[1], 0.500 - 0.500^2)
