# Analysis of benefits and costs of transposon duplication
# Follows notebook entry of 17 June 2023

# Load simulation results
source("./cost_benefit_load.R")


# Check error in target_A_B for simulations
# (Checks whether no-selection simulations [ns] are a good approximation to
# with-selection simulations [ds] in terms of the target A_B)
ns[, target_A_B := 2 * omega / (1 - u)]
ggplot(ns[stage == 1 & g == GE]) +
    geom_point(aes(x = u, y = abs(A_B - target_A_B) / target_A_B, colour = cost)) +
    geom_hline(aes(yintercept = 1e-04)) +
    facet_wrap(~mode) + scale_x_log10() + scale_y_log10()
# Acceptable



###################################
# GRAPHICAL COST-BENEFIT ANALYSIS #
###################################


# Assumes ns, ds, GC, GE, and DD are available in global environment
cost_benefit = function(SimID, Dupl, Mode, Cost, UseVF = FALSE, UseFF = FALSE, ApproxQuantities = TRUE, PubCode = 1, plot = TRUE)
{
    # Read parameters
    if (missing(SimID) && !missing(Dupl) && !missing(Mode) && !missing(Cost)) {
        SimID = ds[uA == Dupl & mode == Mode & cost == Cost, unique(sim_id)]
        if (length(SimID) != 1) {
            stop("Could not find unique SimID with requested properties ", values(Dupl, Mode, Cost));
        }
    } else if (!missing(SimID) && missing(Dupl) && missing(Mode) && missing(Cost)) {
        Dupl = ds[sim_id == SimID, uA[1]];
        Mode = ds[sim_id == SimID, mode[1]];
        Cost = ds[sim_id == SimID, cost[1]];
        if (is.na(Dupl)) {
            stop("Could not find ", values(SimID));
        }
    }

    # Measured cost and benefit of a single additional duplication event among mutants with probability DD
    measured_cost = ds[sim_id == SimID & stage == 1 & tag == "old", -(b[GE] - b[GC])]
    measured_bene = ds[sim_id == SimID & stage == 1 & tag == "new", b[GE]]

    # Cost/benefit per duplication event
    # Note we use the TE number in stage 1, prior to duplication, as those are the 
    # TEs to which the additional duplication occurs
    measured_cost_per_DE = measured_cost / (ds[sim_id == SimID & stage == 1 & tag == "old" & g == GC - 1, b] * DD)
    measured_bene_per_DE = measured_bene / (ds[sim_id == SimID & stage == 1 & tag == "old" & g == GC - 1, b] * DD)

    # Determine/calculate a few key quantities. Mean host copy number at 
    # equilibrium at stage 1 and stage 2:
    eq_n1 = ds[sim_id == SimID & stage == 1 & g == GC - 1 & tag == "new", n];
    # Slope of relative fitness on copy number at stage 2:
    swn2 = ds[sim_id == SimID & stage == 2 & g == GC - 1 & tag == "full", swn];
    if (ApproxQuantities) {
        VF2 = 1 + 2.6 * Dupl + 1.4 * Cost - 1.9 * Dupl * Cost;
        eq_n2 = eq_n1 * (1 + Dupl);
        eq_N2 = eq_n2 + VF2;
        W = 1 / (1 + Dupl);
    } else {
        VF2 = ds[sim_id == SimID & stage == 2 & g == GC - 1 & tag == "new", Vn/n];
        eq_n2 = ds[sim_id == SimID & stage == 2 & g == GC - 1 & tag == "new", n];
        eq_N2 = ds[sim_id == SimID & stage == 2 & g == GC - 1 & tag == "new", N];
        W = ds[sim_id == 62 & stage == 2 & g == GC - 1 & tag == "new", W];
    }
    
    # Calculation of "fudge factors"
    # FF1 is used to lower the Delta N associated with the cost, which here 
    # is calculated too high.
    if (UseFF) {
        FF1 = 0.444 * VF2^2 - 0.877 * VF2 + 1.433;
    } else {
        FF1 = 1;
    }

    # Cost: dN, difference in experienced copy number for "parent" mutant 
    # versus wild-type TEs following a duplication event. 
    # Initial dN is DD * (1 + u) / (1 - u).
    # Then each subsequent step goes by e^s (1 + u) / 2.
    # This is experienced by all original TEs.
    cost_trace = data.table(g = GC:GE);
    cost_trace[, dN_th := DD * (1 + Dupl) / (1 - Dupl) * 1/FF1];

    # Adjusting for recombination, selection, and proliferation. 
    cost_trace[, dN_th := dN_th * ((1 + Dupl) * exp(swn2) / 2) ^ (0:(.N-1))];

    cost_trace = merge(cost_trace,
        ns[sim_id == SimID & stage == 2 & tag == "old", .(g, dN_ns = N_B - N_B[g == GC - 1])],
        by = "g");
    cost_trace = merge(cost_trace,
        ds[sim_id == SimID & stage == 2 & tag == "old", .(g, dN_ds = N_B - N_B[g == GC - 1])],
        by = "g");
    
    # dW, difference in experienced relative fitness for mutant versus
    # wild-type TEs following a duplication event (focusing on new TEs).
    # Notes: for dW_ds, similarly to dN_ds could do W_B - W_A but less accurate.
    cost_trace[, dW_th := dN_th * swn2];
    
    cost_trace = merge(cost_trace,
        ds[sim_id == SimID & stage == 2 & tag == "old", .(g, dW_ds = W_B - W_B[g == GC - 1])],
        by = "g");
    # In order to see how the difference in dW for mutant TEs changes their
    # number over time, we need to make this relative to W for wild-type TEs.
    # Let's notate this relative value of W as Z.
    cost_trace[, dZ_th := dW_th / W];
    cost_trace = merge(cost_trace,
        ds[sim_id == SimID & stage == 2 & tag == "old", 
            .(g, dZ_ds = (W_B - W_B[g == GC - 1]) / W_A[g == GC - 1])],
        by = "g");
    # dZ is now a small negative number, with e.g. -0.01 indicating that mutant
    # TEs have a transposon relative fitness of 0.99; this accounts for both
    # counterselection and proliferation, so e.g. if the time series of dZ was
    # 0.98, 0.99, 1, then at the end of the process the number of mutants would
    # be 0.98*0.99*1 times the initial number of mutants.
    
    # Finally, calculate cost per duplication event.
    # This is (1 - prod(1 + dZ)) / DD, but
    # for accuracy we calculate it as expm1(-sum(log1p(dZ))) / DD
    # since dZ may be very small. 
    cost_trace[, cost_th := -expm1(cumsum(log1p(dZ_th))) / DD];
    cost_trace[, cost_ds := -expm1(cumsum(log1p(dZ_ds))) / DD];

    # Benefit: dN, difference in experienced copy number for "offspring" mutant 
    # versus wild-type TEs following a duplication event. 
    # Initial dN is 1/(1 + u) for private replication. 
    # For public replication, there are a number of options to choose from.
    # Then each subsequent step goes by e^s (1 + u) / 2.
    # This is experienced by a new fraction DD of original TEs.
    # Notes: for dN_ds, could do N_B - N_A instead but because N_A changes 
    # slightly due to nonzero mutants, and eventually equalizes with N_B, the 
    # approach below is more accurate.
    bene_trace = data.table(g = GC:GE);
    if (Mode == "pri") {
        # See WORK - APPROXIMATIONS, "dN for benefit - private replication" below
        bene_trace[, dN_th := 1 / (1 + Dupl)];
    } else {
        # See WORK - APPROXIMATIONS, "dN for benefit - public replication" below.
        dN_calc = with(list(nn = eq_n1, u = Dupl), {
            dN_11 = nn * (1 + u) / (exp(nn) - 1) - u * (2 + u) / (1 + u)
            dN_22 = (exp(nn) * (nn * (1 + 4*u) * (1 + nn + nn*u) - 4*u) - nn + 4*u) / (exp(nn) * (nn + 4 * (nn - 1) * u) - nn + 4*u) - nn * (1 + u) - (1 + u * (1 + u) * (7 + 4*u)) / ((1 + u) * (1 + 2 * u))
            dN_33 = (exp(nn) * (18*u^2*(-nn+1) + nn^2 + 4*u*(nn^2-nn) + 9*u^2*nn^2 + nn^3 + 5*u*nn^3 + 13*u^2*nn^3 + 9*u^3*nn^3) - 18*u^2 + 4*u*nn - nn^2) /
                (exp(nn) * (18*u^2*(-nn+1) + nn^2 + 4*u*(nn^2-nn) + 9*u^2*nn^2) - 18*u^2 + 4*u*nn - nn^2) - (nn + 10/3 + (3 + nn) * u - 1/(1 + u) - (4 + 8*u) / (3 + 6*u + 9*u^2))
            return (list(dN_11 = dN_11, dN_22 = dN_22, dN_33 = dN_33))
        });
        
        if (PubCode == 1) {
            bene_trace[, dN_th := eq_n1 / (exp(eq_n1) - 1) - Dupl];
        } else if (PubCode == 2) {
            bene_trace[, dN_th := eq_n1 / (exp(eq_n1) - 1) - 2 * Dupl];
        } else if (PubCode == 3) {
            bene_trace[, dN_th := eq_n1 / (exp(eq_n1) - 1) - Dupl / (1 + Dupl)];
        } else if (PubCode == 11) {
            bene_trace$dN_th = dN_calc$dN_11;
        } else if (PubCode == 22) {
            bene_trace$dN_th = dN_calc$dN_22;
        } else if (PubCode == 33) {
            bene_trace$dN_th = dN_calc$dN_33;
        }
    }
    
    # Adjusting for recombination. 
    bene_trace[, dN_th := dN_th * ((1 + Dupl) * exp(swn2) / 2) ^ (0:(.N-1))];

    bene_trace = merge(bene_trace,
        ns[sim_id == SimID & stage == 2 & tag == "new", .(g, dN_ns = N_B - N_B[g == GC - 1])],
        by = "g");
    bene_trace = merge(bene_trace,
        ds[sim_id == SimID & stage == 2 & tag == "new", .(g, dN_ds = N_B - N_B[g == GC - 1])],
        by = "g");

    # Add relatedness...
    if (Mode == "pri") {
        bene_trace[, R_th := 1];
        bene_trace[, R_ds := 1];
    } else {
        # Calculation of approximate value for relatedness -- see below
        bene_trace[, R_th := ((exp(-eq_n1) + eq_n1 - 1)*4*Dupl*eq_n1 + 
                (1 - exp(-eq_n1))*eq_n1^2 + (2*(1 - exp(-eq_n1)) + eq_n1^2 - 2*eq_n1)*9*Dupl^2) / 
                ((1 + 2*Dupl + 3*Dupl^2)*eq_n1^3)];
        bene_trace[, R_ds := ds[sim_id == SimID & tag == "new" & stage == 1 & g == GC, r]];
    }

    # dW, difference in experienced relative fitness for mutant versus
    # wild-type TEs following a duplication event (focusing on new TEs).
    # Notes: for dW_ds, similarly to dN_ds could do W_B - W_A but less accurate.
    bene_trace[, dW_th := dN_th * swn2];
    bene_trace = merge(bene_trace,
        ds[sim_id == SimID & stage == 2 & tag == "new", .(g, dW_ds = W_B - W_B[g == GC - 1])],
        by = "g");

    # In order to see how the difference in dW for mutant TEs changes their
    # number over time, we need to make this relative to W for wild-type TEs.
    # Let's notate this relative value of W as Z.
    bene_trace[, dZ_th := dW_th / W];
    bene_trace = merge(bene_trace,
        ds[sim_id == SimID & stage == 2 & tag == "new", 
            .(g, dZ_ds = (W_B - W_B[g == GC - 1]) / W_A[g == GC - 1])],
        by = "g");
    # Finally, calculate benefit. This is R * prod(1 + dZ), but
    # for accuracy we calculate it as R * exp(sum(log1p(dZ))) since dZ may be
    # very small. The factor of W comes in because the average transposon after
    # proliferation has less worth than the average transposon at birth. The 
    # duplication event creates an extra transposon that comes into being
    # after proliferation. But all TEs that exist after proliferation are
    # subject to selection. If the duplication event created a new TE prior
    # to the proliferation phase, then that TE would be worth more as it would
    # have the chance to duplicate itself during proliferation, counteracting
    # the effect of the selection phase.
    bene_trace[, bene_th := R_th * exp(cumsum(log1p(dZ_th))) * W];
    bene_trace[, bene_ds := R_ds * exp(cumsum(log1p(dZ_ds))) * W];

    # Generate plot
    val_data = rbind(
        cost_trace[, .(g, variable = "cost", kind = "Th", value = cost_th)],
        cost_trace[, .(g, variable = "cost", kind = "DS", value = cost_ds)],
        bene_trace[, .(g, variable = "bene", kind = "Th", value = bene_th)],
        bene_trace[, .(g, variable = "bene", kind = "DS", value = bene_ds)]
    );
    
    n_data = rbind(
        cost_trace[, .(g, variable = "cost", kind = "Th", value = dN_th)],
        cost_trace[, .(g, variable = "cost", kind = "DS", value = dN_ds)],
        bene_trace[, .(g, variable = "bene", kind = "Th", value = dN_th)],
        bene_trace[, .(g, variable = "bene", kind = "DS", value = dN_ds)]
    )

    w_data = rbind(
        cost_trace[, .(g, variable = "cost", kind = "Th", value = dW_th)],
        cost_trace[, .(g, variable = "cost", kind = "DS", value = dW_ds)],
        bene_trace[, .(g, variable = "bene", kind = "Th", value = dW_th)],
        bene_trace[, .(g, variable = "bene", kind = "DS", value = dW_ds)]
    )
    
    R = bene_trace[1, R_th];
    Wpar = W + cost_trace$dW_ds;
    Woff = W + bene_trace$dW_ds;
    
    cost_direct = (1 / DD) * (1 - prod(Wpar / W));
    bene_direct = R * W * prod(Woff / W);
    
    cost_approx1 = -2 * swn2 * (1+Dupl)^2 / (1-Dupl);
    cost_approx2 = -2 * swn2 * (1+Dupl)^2 / ((1-Dupl) * (2 - exp(swn2) * (1 + Dupl)));
    bene_approx1 = R * W * (1 + 2*swn2);
    
    # Calculate B
    if (Mode == "pri") {
        bene_approx2 = (1 / (1 + Dupl)) * (1 + 2*swn2 / (2 - exp(swn2) * (1 + Dupl)));
    } else {
        bene_approx2 = with(list(nn = eq_n2, u = Dupl, beta = swn2), {
            # get around issues with exp(nn) overflowing
            if (nn < 100) { 
                dN1 = (exp(nn) * (nn * (1 + 4*u) * (1 + nn + nn*u) - 4*u) - nn + 4*u) / 
                    (exp(nn) * (nn + 4 * (nn - 1) * u) - nn + 4*u);
            } else {
                dN1 = (nn * (1 + 4*u) * (1 + nn + nn*u) - 4*u) / 
                    (nn + 4 * (nn - 1) * u);
            }

            dN = dN1 - 
                nn * (1 + u) - 
                (1 + u * (1 + u) * (7 + 4*u)) / 
                ((1 + u) * (1 + 2 * u));
            R = ((exp(-nn) + nn - 1)*4*u*nn + 
                    (1 - exp(-nn))*nn^2 + (2*(1 - exp(-nn)) + nn^2 - 2*nn)*9*u^2) / 
                    ((1 + 2*u + 3*u^2)*nn^3)
            B = R * (1/(1+u) + (2*beta / (2 - exp(beta) * (1 + u))) * dN);
            return (B)
        });
    }
    
    if (!plot) {
        return (data.table(sim_id = SimID, mode = Mode, cost = Cost, dupl = Dupl,
            effect = c("cost", "benefit", "cost", "benefit", "cost", "benefit"), 
            method = c("measured", "measured", "approx1", "approx1", "approx2", "approx2"),
            value = c(measured_cost_per_DE, measured_bene_per_DE, cost_approx1, bene_approx1, cost_approx2, bene_approx2)))
    }
    
    p1 = ggplot(val_data) +
        geom_line(aes(g, value, colour = variable, linetype = kind)) +
        annotate("point", x = 147, y = measured_cost_per_DE, colour = "#55cccc", shape = 4, size = 4) +
        annotate("point", x = 153, y = measured_bene_per_DE, colour = "#ff6666", shape = 4, size = 4) +
        annotate("point", x = 172, y = cost_approx1, colour = "#55cccc", shape = 0, size = 4) +
        annotate("point", x = 178, y = bene_approx1, colour = "#ff6666", shape = 0, size = 4) +
        annotate("point", x = 172, y = cost_approx2, colour = "#55cccc", shape = 1, size = 4) +
        annotate("point", x = 178, y = bene_approx2, colour = "#ff6666", shape = 1, size = 4) +
        annotate("point", x = 147, y = cost_direct, colour = "#55cccc", shape = 1, size = 3) +
        annotate("point", x = 153, y = bene_direct, colour = "#ff6666", shape = 1, size = 3) +
        labs(x = "Generation", y = "Value", title = "Cost and benefit", subtitle = values(SimID, Mode, Dupl, Cost)) +
        cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "bottom") +
        ylim(0, NA)
    
    n_annotations = n_data[g %in% c(100, 105),
        .(diff = round(value[kind == "Th"] / value[kind == "DS"], 4),
            y = value[kind == "Th"] * 0.5), by = .(g, variable)];
    
    p2 = ggplot(n_data[g <= 105]) +
        geom_line(aes(g, value, colour = variable, linetype = kind)) +
        geom_label(data = n_annotations, aes(g, y, label = diff), size = 2, label.padding = unit(0.1, "lines")) +
        labs(x = "Generation", y = "Value", title = "N_B excess") +
        facet_wrap(~variable, scales = "free_y", nrow = 2) +
        cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "bottom") 

    w_annotations = w_data[g %in% c(100, 105),
        .(diff = round(value[kind == "Th"] / value[kind == "DS"], 4),
            y = value[kind == "Th"]), by = .(g, variable)];

    p3 = ggplot(w_data[g <= 105]) +
        geom_line(aes(g, -value, colour = variable, linetype = kind)) +
        geom_label(data = w_annotations, aes(g, -y*0.5, label = diff), size = 2, label.padding = unit(0.1, "lines")) +
        labs(x = "Generation", y = "Value", title = "W_B deficit") +
        facet_wrap(~variable, scales = "free_y", nrow = 2) +
        cowplot::theme_cowplot(font_size = 8) +
        theme(legend.position = "bottom")

    cowplot::plot_grid(p1, p2, p3, nrow = 1)
}

# Generate full set of cost-benefit figures in one PDF
pdf(paste0("./cost_benefit.pdf"), width = 12, height = 4, useDingbats = FALSE)
for (i in 1:124) {
    print(cost_benefit(i, UseVF = FALSE, ApproxQuantities = TRUE, PubCode = 22, UseFF = FALSE))
}
dev.off()

results = list()
for (i in 1:124) {
    results[[i]] = cost_benefit(i, UseVF = FALSE, ApproxQuantities = TRUE, PubCode = 22, UseFF = FALSE, plot = FALSE)
}
results = rbindlist(results)

ggplot(results[mode == "pri"]) +
    geom_line(aes(x = dupl, y = value, linetype = effect, colour = method)) +
    facet_wrap(~cost, scales = "free_x", ncol = 1)


# PLOTS FOR DIDACTIC FIGURE ON HOW THE INVASION ANALYSIS WORKS

SimID = 24
b0 = ds[sim_id == SimID & tag == "new" & g == GC & stage == 1, b]
dn = ds[sim_id == SimID & tag == "new" & g %between% c(GC - 1, GC + 10) & stage == 2, 
    .(g = g - GC, N_A, N_B = c(NA, tail(N_B, -1)), W_A, W_B = c(NA, tail(W_B, -1)), db = c(NA, tail(b, -1)))]
dn[, V_B := db / (DD * b0)]

do = ds[sim_id == SimID & tag == "old" & g %between% c(GC - 1, GC + 10) & stage == 2, 
    .(g = g - GC, N_A, N_B = c(NA, tail(N_B, -1)), W_A, W_B = c(NA, tail(W_B, -1)), db = c(NA, tail(b - b[1], -1)))]
do[, N_B := (N_B - N_A) / DD + N_A]
do[, W_B := (W_B - W_A) / DD + W_A]
do[, V_B := -db / (DD * b0)]

theme_set(cowplot::theme_cowplot(font_size = 8))

p1 = ggplot(do) +
    geom_line(aes(x = g, y = N_A, colour = "Wild-type", linetype = "Wild-type", linewidth = "Wild-type")) + 
    geom_line(aes(x = g, y = N_B, colour = "Mutant", linetype = "Mutant", linewidth = "Mutant")) +
    scale_colour_manual(values = c("Wild-type" = "black", "Mutant" = "#bb5555")) +
    scale_linetype_manual(values = c("Wild-type" = "42", "Mutant" = "solid")) +
    scale_linewidth_manual(values = c("Wild-type" = 0.5, "Mutant" = 0.5)) +
    theme(legend.position = c(1, 1), legend.justification = c(1, 1), legend.key.width = unit(0.45, "cm"), plot.title = element_text(face = "plain", size = 8)) +
    scale_x_continuous(limits = c(-1, 10), breaks = c(0, 2, 4, 6, 8, 10)) +
    scale_y_continuous(limits = c(5.9, 7.2), breaks = c(6, 6.5, 7)) +
    labs(x = "Generations after\nduplication event", y = "Copy number", colour = NULL, linetype = NULL, linewidth = NULL, title = ~"Existing mutant: " * italic(N))

p2 = ggplot(do) +
    geom_line(aes(x = g, y = W_A, colour = "Wild-type", linetype = "Wild-type", linewidth = "Wild-type")) + 
    geom_line(aes(x = g, y = W_B, colour = "Mutant", linetype = "Mutant", linewidth = "Mutant")) +
    scale_colour_manual(values = c("Wild-type" = "black", "Mutant" = "#bb5555")) +
    scale_linetype_manual(values = c("Wild-type" = "42", "Mutant" = "solid")) +
    scale_linewidth_manual(values = c("Wild-type" = 0.5, "Mutant" = 0.5)) +
    theme(legend.position = "none", plot.title = element_text(face = "plain", size = 8)) +
    scale_x_continuous(limits = c(-1, 10), breaks = c(0, 2, 4, 6, 8, 10)) +
    scale_y_continuous(limits = c(0.937, 0.973)) +
    labs(x = "Generations after\nduplication event", y = "Relative fitness", colour = NULL, linetype = NULL, linewidth = NULL, title = ~"Existing mutant: " * italic(W))

p3 = ggplot(do) +
    geom_line(aes(x = g, y = V_B), colour = "#bb5555") +
    labs(x = "Generations after\nduplication event", y = "Reproductive\nvalue", title = "Cost") +
    theme(plot.title = element_text(face = "plain", size = 8)) +
    scale_x_continuous(limits = c(-1, 10), breaks = c(0, 2, 4, 6, 8, 10)) +
    ylim(-0.01, 1.01)

p4 = ggplot(dn) +
    geom_line(aes(x = g, y = N_A, colour = "Wild-type", linetype = "Wild-type", linewidth = "Wild-type")) + 
    geom_line(aes(x = g, y = N_B, colour = "Mutant", linetype = "Mutant", linewidth = "Mutant")) +
    scale_colour_manual(values = c("Wild-type" = "black", "Mutant" = "#c9b8ec")) +
    scale_linetype_manual(values = c("Wild-type" = "42", "Mutant" = "solid")) +
    scale_linewidth_manual(values = c("Wild-type" = 0.5, "Mutant" = 0.5)) +
    theme(legend.position = c(1, 1), legend.justification = c(1, 1), legend.key.width = unit(0.45, "cm"), plot.title = element_text(face = "plain", size = 8)) +
    scale_x_continuous(limits = c(-1, 10), breaks = c(0, 2, 4, 6, 8, 10)) +
    scale_y_continuous(limits = c(5.9, 7.2), breaks = c(6, 6.5, 7)) +
    labs(x = "Generations after\nduplication event", y = "Copy number", colour = NULL, linetype = NULL, linewidth = NULL, title = ~"New mutant: " * italic(N))

p5 = ggplot(dn) +
    geom_line(aes(x = g, y = W_A, colour = "Wild-type", linetype = "Wild-type", linewidth = "Wild-type")) + 
    geom_line(aes(x = g, y = W_B, colour = "Mutant", linetype = "Mutant", linewidth = "Mutant")) +
    scale_colour_manual(values = c("Wild-type" = "black", "Mutant" = "#c9b8ec")) +
    scale_linetype_manual(values = c("Wild-type" = "42", "Mutant" = "solid")) +
    scale_linewidth_manual(values = c("Wild-type" = 0.5, "Mutant" = 0.5)) +
    theme(legend.position = "none", plot.title = element_text(face = "plain", size = 8)) +
    scale_x_continuous(limits = c(-1, 10), breaks = c(0, 2, 4, 6, 8, 10)) +
    scale_y_continuous(limits = c(0.937, 0.973)) +
    labs(x = "Generations after\nduplication event", y = "Relative fitness", colour = NULL, linetype = NULL, linewidth = NULL, title = ~"New mutant: " * italic(W))

p6 = ggplot(dn) +
    geom_line(aes(x = g, y = V_B), colour = "#c9b8ec") +
    labs(x = "Generations after\nduplication event", y = "Reproductive\nvalue", title = "Benefit") +
    theme(plot.title = element_text(face = "plain", size = 8)) +
    scale_x_continuous(limits = c(-1, 10), breaks = c(0, 2, 4, 6, 8, 10)) +
    ylim(-0.01, 1.01)

pp = cowplot::plot_grid(p1, p2, p3, p4, p5, p6, nrow = 2, labels = c("B", "D", "F", "C", "E", "G"), label_size = 9, rel_widths = c(3.75, 3.75, 4))
ggsave("./Figures/S-invasion-analysis-illust-inset.pdf", pp, width = 11.5, height = 7, units = "cm", useDingbats = FALSE)

# Heat map
genotypes = fread(paste0("./Deterministic/Runs/6-Genotypes/pri_c005_u030_old.ngd"), skip = "g\t");
genotypes = genotypes[g == GC];

p_heat = ggplot(genotypes[A <= 12 & B <= 12]) +
    geom_tile(aes(x = A, y = B, fill = log10(f))) +
    scale_fill_viridis_c() +
    scale_x_continuous(expand = c(0, 0), breaks = seq(0, 12, 2)) +
    scale_y_continuous(expand = c(0, 0), breaks = seq(0, 12, 2)) +
    labs(x = "Wild-type", y = "Mutant", fill = ~log[10] ~ "freq.") +
    theme(legend.position = c(0.65, 0.75), legend.box.background = element_rect(fill = "white", colour = "black"), legend.box.margin = margin(3, 3, 3, 3))

p_wt = ggplot(genotypes[A <= 12][, .(sum(f)), by = A]) +
    geom_col(aes(x = A, y = V1)) +
    scale_x_continuous(expand = c(0, 0), breaks = seq(0, 12, 2)) +
    scale_y_continuous(expand = c(0, 0)) +
    labs(x = NULL, y = "Freq.")

p_mut = ggplot(genotypes[B <= 12][, .(sum(f)), by = B]) +
    geom_col(aes(x = B, y = V1)) +
    annotate("rect", xmin = -0.5, ymin = 1e-5, xmax = 0.5, ymax = 0.05, fill = "white") +
    scale_x_continuous(expand = c(0, 0), breaks = seq(0, 12, 2)) +
    scale_y_continuous(trans = scales::trans_new(
            name = "custom",
            transform = function(x) ifelse(x < 1e-5, x * 8e4, 0.8 + x/5),
            inverse = function(x) ifelse(x < 0.8, x / 8e4, (x - 0.8) * 5)
        ), breaks = c(seq(0, 1e-5, 5e-6), 0.5, 1),
        labels = expression(0, ~5 %*% 10^-6, ~1 %*% 10^-5, 0.5, 1),
        expand = c(0, 0),
        limits = c(0, 1)) +
    labs(x = NULL, y = "Freq.") +
    coord_flip()

p_heat_multi = cowplot::plot_grid(p_wt, ggdraw(), p_heat, p_mut, nrow = 2,
    rel_widths = c(3, 1.222222), rel_heights = c(1, 3), align = "hv")

ggsave("./Figures/S-invasion-analysis-illust-heat.pdf", p_heat_multi, width = 9.5, height = 9, units = "cm", useDingbats = FALSE)
ggsave("./Figures/S-invasion-analysis-illust-Nplot.pdf", p4, width = 6, height = 5, units = "cm", useDingbats = FALSE)


Nw = genotypes[, max(A) + 1];
Nm = genotypes[, max(B) + 1];
x = matrix(0, Nw, Nm);
x[cbind(genotypes$A + 1, genotypes$B + 1)] = genotypes$f;

# Check of wild-type distribution
n = weighted.mean(0:(Nw - 1), rowSums(x)*(0:(Nw-1)));
cw = data.table(x = 0:(Nw - 1), pm = rowSums(x), pt = dpois(0:(Nw - 1), n));

# Check of mutant distribution
mx = colSums(x)[-1];
mx = mx / sum(mx);
cm = data.table(x = 1:(Nm - 1), pm = mx, pt = dquasi1(1:(Nw - 1), num_u));


######################################
# DERIVING APPROXIMATIONS USED ABOVE #
######################################


##########################
# Approximation for Vn/n #
##########################

# Get data table with deterministic simulation measurements
x = ds[tag == "old" & stage == 2 & g == GC - 1, .(sim_id, u = uA, cost, Vn, n)];

# Derive linear approximation -- we want this to pass through 1 at u = 0 & cost = 0
glm(I(Vn/n - 1) ~ 0 + u * cost, data = x)

# Coefficients:
#      u    cost  u:cost  
#  2.597   1.413  -1.904

# Coefficients can be rounded:
ggplot(x) + 
    geom_point(aes(1 + 2.6 * u + 1.4 * cost - 1.9 * u * cost, Vn/n, colour = log10(cost))) + 
    geom_abline(aes(slope = 1, intercept = 0))

# Final result: Vn/n ~= 1 + 2.6u + 1.4c - 1.9uc

# A sort of "fine" approximation for Vn/n for stage == 1 is 1 + u/(1+u):
x = ds[tag == "old" & stage == 1 & g == GC - 1, .(sim_id, u = uA, cost, Vn, n)];
ggplot(x) + 
    geom_point(aes(1 + u / (1 + u), Vn/n, colour = log10(cost))) + 
    geom_abline(aes(slope = 1, intercept = 0))

# Or, can we derive VF1 from VF2?
# Yes. VF1 is roughly VF2 minus 2*u.
x = ds[tag == "old" & stage %in% c(1, 2) & g == GC - 1, .(u = uA[1], cost = cost[1], vf1 = Vn[1]/n[1], vf2 = Vn[2]/n[2]), by = sim_id];
ggplot(x) + 
    geom_point(aes(vf2 - 2*u, vf1, colour = log10(cost))) + 
    geom_abline(aes(slope = 1, intercept = 0))

# Then, can bring that into the approximation for VF2 derived above.
ggplot(x) + 
    geom_point(aes(1 + 0.6 * u + 1.4 * cost - 1.9 * u * cost, vf1, colour = log10(cost))) + 
    geom_abline(aes(slope = 1, intercept = 0))

ggplot(x) + 
    geom_point(aes(1 + 0.63 * u + 1.3 * cost - 2.2 * u * cost, vf1, colour = log10(cost))) + 
    geom_abline(aes(slope = 1, intercept = 0))

ggplot(x) + 
    geom_point(aes(vf1, (1 + 0.6 * u + 1.4 * cost - 1.9 * u * cost - vf1)/vf1, colour = log10(cost)))

ggplot(x) + 
    geom_point(aes(vf1, (1 + 0.6 * u + 1.4 * cost - 1.9 * u * cost - vf1)/(vf1 - 1), colour = log10(cost)))


#####################################
# Approximation for dW/dn and dW/dN #
#####################################

# Note: also check for GC, GC+1 (perhaps)
x = ds[mode == "pri" & g == GC - 1 & stage == 1 & tag == "old", .(sim_id, u = uA, cost, swn, sWN, sWN_B, n, N, a, VF = Vn/n)]

x[, an1 := dWdn(n, cost)]
x[, an2 := dWdn_poisson(n, cost)]
x[, an3 := dWdn_poisson(n, cost, VF)]
x[, aN1 := dWdn(N, cost)]
x[, aN2 := dWdn_poisson(N, cost)]
x[, aN3 := dWdn_poisson(N, cost, VF)]
x[, aN4 := dWdN_mutant(a, u, cost, VF)]

ggplot(melt(x, measure.vars = names(x)[names(x) %like% "^an"])) +
    geom_abline(aes(intercept = 0, slope = 1)) +
    geom_path(aes(x = swn, y = value, colour = variable, group = cost, linetype = as.factor(cost))) +
    facet_wrap(~variable, nrow = 2)

ggplot(melt(x, measure.vars = names(x)[names(x) %like% "^aN"])) +
    geom_abline(aes(intercept = 0, slope = 1)) +
    geom_path(aes(x = sWN, y = value, colour = variable, group = cost, linetype = as.factor(cost))) +
    facet_wrap(~variable)


#################################
# Approximation for relatedness #
#################################

# One approximation for relatedness is just B_B / N_B, but that is E(B)/E(N) 
# rather than E(B/N) so it doesn't quite work. Here I derive an approximation
# for relatedness analytically.

# Calculate relatedness numerically based on a Poisson distribution of wild-type 
# TEs (with mean wt_n) and a quasigeometric distribution of mutant TEs (with 
# duplication rate parameter u).
relatedness = function(wt_n, u)
{
    # Choose an upper bound for wild-type TEs
    max_n = max(20, 3 * wt_n)
    
    # Take a weighted mean of j/(i+j), i.e. relatedness, where j is mutant TE
    # number and i is wild-type TE number, from the mutant TE perspective.
    s = 0;
    w = 0;
    for (i in 0:max_n) {
        for (j in 1:50) {
            # Below, j/(i+j) is relatedness, dpois is wild-type density,
            # dquasi1 is mutant density, and j is number of mutant observers
            s = s + j/(i+j) * dpois(i, wt_n) * dquasi1(j, u) * j;
            # Sum up weights of each observation that goes into s
            w = w + dpois(i, wt_n) * dquasi1(j, u) * j;
        }
    }
    
    # Return weighted mean
    s/w
}

# Get deterministic simulation data, including "true" measured relatedness r
# and "naive" approximation B_B / N_B
x = ds[tag == "new" & stage == 1 & g == GC & mode == "pub", .(sim_id, cost, u = uA, n, r, r_naive = B_B / N_B, vf = Vn/n)]
x[, k := n/(vf - 1)]

# Approximation based on numerical calculation
x[, r_num := mapply(relatedness, n, u)]

# Obtained by assuming i follows a Poisson distribution, and j can only be 1, 1:2, or 1:3
x[, r_approx_pois1 := (1 - exp(-n))/n]
x[, r_approx_pois2 := (n + 4*u*(n - 1) - exp(-n)*(n - 4*u))/((1 + 2*u)*n^2)]
x[, r_approx_pois3 := ((exp(-n) + n - 1)*4*u*n + (1 - exp(-n))*n^2 + (2*(1 - exp(-n)) + n^2 - 2*n)*9*u^2) / ((1 + 2*u + 3*u^2)*n^3)]

# Negative-binomial versions of the above
x[, r_approx_nb1 := (k - (k/(k + n))^k * (k + n))/((-1 + k)*n)]
x[, r_approx_nb2 := ((k/(k + n))^k * (k + n)*(4*k*u + 2*n - k*n + 4*u*n) + k*(-4*k*u + (-2 + k)*(1 + 4*u)*n)) /
    ((-2 + k)*(-1 + k)*(1 + 2*u)*n^2)]
x[, r_approx_nb3 := (-(k/(k + n))^k * (k + n)*(18*k^2*u^2 - 4*k*(-3 + k - 9*u)*u*n + (6 + k^2 + 6*u*(2 + 3*u) - k*(5 + 4*u))*n^2) + 
    k*(18*k^2*u^2 - 2*(-3 + k)*k*u*(2 + 9*u)*n + (-3 + k)*(-2 + k)*(1 + u*(4 + 9*u))*n^2))/
    ((-3 + k)*(-2 + k)*(-1 + k)*(1 + u*(2 + 3*u))*n^3)]

# Check approximations visually
x = melt(x, measure.vars = names(x)[names(x) %like% "^r_"]);

ggplot(x) +
    geom_abline(aes(intercept = 0, slope = 1)) +
    geom_line(aes(x = r, y = value, colour = variable, group = cost, linetype = as.factor(cost))) +
    facet_wrap(~variable)

ggplot(x[variable %in% c("r_approx_pois3", "r_approx_nb3")]) +
    geom_abline(aes(intercept = 0, slope = 1)) +
    geom_line(aes(x = r, y = value, colour = variable, group = paste(cost, variable), linetype = as.factor(cost)))



########################################
# dN for benefit - private replication #
########################################

# After duplication event
nn = 5
u = 0.001
d = 1e-4

N = 0; D = 0;
for (n in 0:50) {
  for (m in 1:10) {
    for (i in 0:10) {
      for (j in 0:4) {
        for (s in 0:m) {
          N = N + dpois(n, nn) * dquasi1(m, u) * dpois(i, n*u) * dpois(j, m*u) * dbinom(s, m, d) * s * (n + m + i + j + s)
          D = D + dpois(n, nn) * dquasi1(m, u) * dpois(i, n*u) * dpois(j, m*u) * dbinom(s, m, d) * s
        }
      }
    }
  }
}
N/D

N = 0; D = 0;
for (n in 0:50) {
  for (m in 1:10) {
    for (i in 0:10) {
      for (j in 0:4) {
        N = N + dpois(n, nn) * dquasi1(m, u) * dpois(i, n*u) * dpois(j, m*u) * (d * m * (n + m + i + j) + d^2 * m^2 + d*(1 - d) * m)
        D = D + dpois(n, nn) * dquasi1(m, u) * dpois(i, n*u) * dpois(j, m*u) * d * m
      }
    }
  }
}
N/D

N = 0; D = 0;
for (n in 0:50) {
  for (m in 1:10) {
    for (i in 0:10) {
      N = N + dpois(n, nn) * dquasi1(m, u) * dpois(i, n*u) * (d * m * (n + m + i + m*u) + d^2 * m^2 + d*(1 - d) * m)
      D = D + dpois(n, nn) * dquasi1(m, u) * dpois(i, n*u) * d * m
    }
  }
}
N/D

N = 0; D = 0;
for (n in 0:50) {
  for (m in 1:10) {
    N = N + dpois(n, nn) * dquasi1(m, u) * (d * m * (n + m + n*u + m*u) + d^2 * m^2 + d*(1 - d) * m)
    D = D + dpois(n, nn) * dquasi1(m, u) * d * m
  }
}
N/D

N = 0; D = 0;
for (n in 0:50) {
  for (m in 1:10) {
    N = N + dpois(n, nn) * dquasi1(m, u) * (m * d * (n + n * u + 1 - d) + m^2 * (d^2 + d * (1 + u)))
    D = D + dpois(n, nn) * dquasi1(m, u) * m * d
  }
}
N/D

N = 0; D = 0;
for (n in 0:50) {
    N = N + dpois(n, nn) * (1/(1-u) * d * (n + n * u + 1 - d) + (1 + u) / (1 - u)^2 * (d^2 + d * (1 + u)))
    D = D + dpois(n, nn) * d * 1 / (1 - u)
}
N/D

N = (1/(1-u) * d * (nn + nn * u + 1 - d) + (1 + u) / (1 - u)^2 * (d^2 + d * (1 + u)))
D = d / (1 - u)
N/D

# Fully simplified result:
Post = (2 + nn + u + 2 * d * u + u^2 - nn * u^2)/(1 - u)

# Before duplication event
N = 0; D = 0;
for (n in 0:50) {
  for (m in 1:10) {
    for (i in 0:10) {
      for (j in 0:5) {
        N = N + dpois(n, nn) * dquasi1(m, u) * dpois(i, n*u) * dpois(j, m*u) * (m + j) * (n + m + i + j)
        D = D + dpois(n, nn) * dquasi1(m, u) * dpois(i, n*u) * dpois(j, m*u) * (m + j)
      }
    }
  }
}
N/D

N = 0; D = 0;
for (n in 0:50) {
  for (m in 1:10) {
    for (i in 0:10) {
      N = N + dpois(n, nn) * dquasi1(m, u) * dpois(i, n*u) * (m * (n + m + i) + (m*u) * (n + 2*m + i) + (m*u * (m*u + 1)))
      D = D + dpois(n, nn) * dquasi1(m, u) * dpois(i, n*u) * (m + m*u)
    }
  }
}
N/D

N = 0; D = 0;
for (n in 0:50) {
  for (m in 1:10) {
    N = N + dpois(n, nn) * dquasi1(m, u) * (m^2 * (1 + 2*u + u^2) + m * (n * (1 + 2*u + u^2) + u))
    D = D + dpois(n, nn) * dquasi1(m, u) * m * (1 + u)
  }
}
N/D

N = 0; D = 0;
for (n in 0:50) {
    N = N + dpois(n, nn) * ((1+u)/(1-u)^2 * (1 + 2*u + u^2) + 1/(1-u) * (n * (1 + 2*u + u^2) + u))
    D = D + dpois(n, nn) * 1/(1-u) * (1 + u)
}
N/D

N = (1+u)/(1-u)^2 * (1 + 2*u + u^2) + 1/(1-u) * (nn * (1 + 2*u + u^2) + u)
D = 1/(1-u) * (1 + u)
N/D

# Fully simplified result:
Pre = nn * (1 + u) + (1 + u * (4 + u * (2 + u)))/(1 - u^2)

# Then dN = Post - Pre, which simplifies to:
1/(1+u) + 2*d*u/(1-u)
# Setting d to 0 yields:
1/(1+u)

((1+u)/(1-u)^2*(1+u)^2+1/(1-u)*(nn*(1+u)^2+u))/((1+u)/(1-u))

((1+u)^3 + (1-u)*(nn*(1+u)^2+u)) / ((1-u)*(1+u))

((1+u)^2 + (1-u)*(nn*(1+u)+u/(1+u))) / (1-u)

(1 + 2*u + u^2 + nn*(1-u^2) + u*(1-u)/(1+u)) / (1-u)


#######################################
# dN for benefit - public replication #
#######################################

### After duplication event ###

nn = 2.5
u = 0.5
d = 1e-8

N = 0; D = 0;
for (n in 0:50) {
  for (m in 1:30) {
    for (i in 0:50) {
      for (j in 0:30) {
          s = 0:m
          N = N + dpois(n, nn) * dquasi1(m, u) * dpois(i, n*u) * dpois(j, m*u) * sum(dbinom(s, m, d * m / (n + m)) * s * (n + m + i + j + s))
          D = D + dpois(n, nn) * dquasi1(m, u) * dpois(i, n*u) * dpois(j, m*u) * sum(dbinom(s, m, d * m / (n + m)) * s)
      }
    }
  }
}
N/D

N = 0; D = 0;
for (n in 0:50) {
  for (m in 1:30) {
    for (i in 0:50) {
      for (j in 0:30) {
          N = N + dpois(n, nn) * dquasi1(m, u) * m * dpois(i, n*u) * dpois(j, m*u) * m/(n+m) * (n + m + i + j + 1)
          D = D + dpois(n, nn) * dquasi1(m, u) * m * dpois(i, n*u) * dpois(j, m*u) * m/(n+m)
      }
    }
  }
}
N/D

N = 0; D = 0;
for (n in 0:50) {
  for (m in 1:30) {
    for (i in 0:50) {
      for (j in 0:30) {
        N = N + dpois(n, nn) * dquasi1(m, u) * dpois(i, n*u) * dpois(j, m*u) * (d^2 * m^4 / (n + m)^2 + (d * m / (n + m)) * (1 - d * m / (n + m)) * m + d * m^2 / (n + m) * (n + m + i + j))
        D = D + dpois(n, nn) * dquasi1(m, u) * dpois(i, n*u) * dpois(j, m*u) * d * m^2 / (n + m)
      }
    }
  }
}
N/D

N = 0; D = 0;
for (n in 0:50) {
  for (m in 1:30) {
    for (i in 0:50) {
      N = N + dpois(n, nn) * dquasi1(m, u) * dpois(i, n*u) * (d^2 * m^4 / (n + m)^2 + (d * m / (n + m)) * (1 - d * m / (n + m)) * m + d * m^2 / (n + m) * (n + m + i + m * u))
      D = D + dpois(n, nn) * dquasi1(m, u) * dpois(i, n*u) * d * m^2 / (n + m)
    }
  }
}
N/D

N = 0; D = 0;
for (n in 0:100) {
  for (m in 1:60) {
    N = N + dpois(n, nn) * dquasi1(m, u) * (d^2 * m^4 / (n + m)^2 + (d * m / (n + m)) * (1 - d * m / (n + m)) * m + d * m^2 / (n + m) * (n + m + n * u + m * u))
    D = D + dpois(n, nn) * dquasi1(m, u) * d * m^2 / (n + m)
  }
}
N/D

# That might be a bit insoluble. Let's divide through by d and eliminate any terms in d still remaining.
N = 0; D = 0;
for (n in 0:80) {
  for (m in 1:60) {
    N = N + dpois(n, nn) * dquasi1(m, u) * (m^2 / (n + m) + m^2 / (n + m) * (n + m + n * u + m * u))
    D = D + dpois(n, nn) * dquasi1(m, u) * m^2 / (n + m)
  }
}
N/D
PostNum = N/D

# OK, but I don't think I can solve that. 
# Let's take just m = 1:3.
N = 0; D = 0;
for (n in 0:50) {
    N = N + dpois(n, nn) * (1 + (1 - u) / (n + 1) + (4 + 4 * (1 - u) / (n + 2)) * u + (8 + 9 * (1 - u) / (n + 3)) * u^2 - 4 * u^3 - 9 * u^4)
    D = D + dpois(n, nn) * (1 - u) * (1/(n + 1) + 4 * u/(n + 2) + 9 * u^2 / (n + 3))
}
N/D

N = -18*u^2 + 18*exp(nn)*u^2*(-nn+1) + 4*u*nn - nn^2 + exp(nn)*nn^2 + 4*exp(nn)*u*(nn^2-nn) + 9*exp(nn)*u^2*nn^2 + exp(nn)*nn^3 + 5*exp(nn)*u*nn^3 + 13*exp(nn)*u^2*nn^3 + 9*exp(nn)*u^3*nn^3
D = -18*u^2 + 18*exp(nn)*u^2*(-nn+1) + 4*u*nn - nn^2 + exp(nn)*nn^2 + 4*exp(nn)*u*(nn^2-nn) + 9*exp(nn)*u^2*nn^2
N/D

# Fully simplified result (m = 1:3)
N = exp(nn) * (18*u^2*(-nn+1) + nn^2 + 4*u*(nn^2-nn) + 9*u^2*nn^2 + nn^3 + 5*u*nn^3 + 13*u^2*nn^3 + 9*u^3*nn^3) - 18*u^2 + 4*u*nn - nn^2
D = exp(nn) * (18*u^2*(-nn+1) + nn^2 + 4*u*(nn^2-nn) + 9*u^2*nn^2) - 18*u^2 + 4*u*nn - nn^2
Post = N/D

# Can also try m = 1:2
N = 0; D = 0;
for (n in 0:50) {
    N = N + dpois(n, nn) * ((2 + n)/(1 + n) + ((10 + 15*n + 4*n^2)*u)/((1 + n)*(2 + n)) - ((6 + n)*u^2)/(2 + n) - 4*u^3)
    D = D + dpois(n, nn) * ((1 - u)*(1/(1 + n) + (4*u)/(2 + n)))
}
N/D

# Intermediate steps from SI
N = 0; D = 0;
for (n in 0:50) {
    N = N + dpois(n, nn) * ((1-u)*(n*(1+u)+2+u)/(n+1) + 4*(1-u)*u*(n*(1+u)+3+2*u)/(n+2))
    D = D + dpois(n, nn) * (4*u*(1-u)/(n+2) + (1-u)/(n+1))
}
N/D

N = 0; D = 0;
for (n in 0:50) {
    N = N + dpois(n, nn) * (4 + 10*u - 6*u^2 - 8*u^3 + n * (4 + 15*u - 7*u^2 - 12*u^3) + n^2 * (1 + 4*u - u^2 - 4*u^3))/((n+1)*(n+2))
    D = D + dpois(n, nn) * ((1-u)*(2+4*u) + n * (1-u)*(1+4*u))/((n+1)*(n+2))
}
N/D

N = 0; D = 0;
for (n in 0:50) {
    N = N + nn^n * exp(-nn) / factorial(n+2) * (4 + 10*u - 6*u^2 - 8*u^3 + n * (4 + 15*u - 7*u^2 - 12*u^3) + n^2 * (1 + 4*u - u^2 - 4*u^3))
    D = D + nn^n * exp(-nn) / factorial(n+2) * ((1-u)*(2+4*u) + n * (1-u)*(1+4*u))
}
N/D

N = 0; D = 0;
for (n in 0:50) {
    N = N + nn^-2 * nn^(n+2) * exp(-nn) / factorial(n+2) * (4 + 10*u - 6*u^2 - 8*u^3 + n * (4 + 15*u - 7*u^2 - 12*u^3) + n^2 * (1 + 4*u - u^2 - 4*u^3))
    D = D + nn^-2 * nn^(n+2) * exp(-nn) / factorial(n+2) * ((1-u)*(2+4*u) + n * (1-u)*(1+4*u))
}
N/D

N = 0; D = 0;
for (n in 0:50) {
    N = N + nn^(n+2) * exp(-nn) / factorial(n+2) * (4 + 10*u - 6*u^2 - 8*u^3 + n * (4 + 15*u - 7*u^2 - 12*u^3) + n^2 * (1 + 4*u - u^2 - 4*u^3))
    D = D + nn^(n+2) * exp(-nn) / factorial(n+2) * ((1-u)*(2+4*u) + n * (1-u)*(1+4*u))
}
N/D


N = (exp(-nn)*(1 - u)*(4*u - 4*exp(nn)*u - nn + exp(nn)*nn + 4*exp(nn)*u*nn + exp(nn)*nn^2 + 5*exp(nn)*u*nn^2 + 4*exp(nn)*u^2*nn^2))/nn^2
D = (exp(-nn)*(1 - u)*(4*u*(1 + exp(nn)*(nn - 1)) + (exp(nn) - 1)*nn))/nn^2
N/D

# Fully simplified result (m = 1:2)
N = 4 * u - 4 * exp(nn) * u - nn + exp(nn) * nn + 4 * exp(nn) *u * nn + exp(nn) * nn^2 + 5 * exp(nn) * u * nn^2 + 4 * exp(nn) * u^2 * nn^2
D = 4 * u * (1 + exp(nn) * (nn - 1)) + (exp(nn) - 1) * nn
N/D
Post = (exp(nn)*(nn*(1 + 4*u)*(1 + nn + nn*u) - 4*u) - nn + 4*u) / (exp(nn)*(nn + 4*(nn - 1)*u) - nn + 4*u)



(nn - 4*u - exp(nn)*(nn*(nn*(1 + 5*u) + 1) + 4*u*(nn*(1 + nn*u) - 1))) / (nn - 4*u - exp(nn)*(nn - 4*u*(nn + 1)))

# Or m = 1:1 only
N = 0;
D = 0;
for (n in 0:50) {
    N = N + dpois(n, nn) * ((1 - u)*(2 + n + u + n*u))/(1 + n)
    D = D + dpois(n, nn) * (1 - u)/(1 + n)
}
N/D

N = (exp(-nn) * (1 - u) * (-1 + exp(nn) * (1 + nn + u*nn)))/nn
D = (exp(-nn) * (1 - u) * (-1 + exp(nn)))/nn
N/D

# Fully simplified result (m = 1:1)
N = exp(nn) * (nn + u*nn + 1) - 1
D = exp(nn) - 1
N/D
Post = (exp(nn) * (nn + u*nn + 1) - 1) / (exp(nn) - 1)


# Or, assuming negative binomial distribution
vf = 1.2;
k = nn/(vf - 1);

N = 0; D = 0;
for (n in 0:80) {
  for (m in 1:60) {
    N = N + dnbinom1(n, nn, vf) * dquasi1(m, u) * (m^2 / (n + m) + m^2 / (n + m) * (n + m + n * u + m * u))
    D = D + dnbinom1(n, nn, vf) * dquasi1(m, u) * m^2 / (n + m)
  }
}
N/D

# Negative binomial, m = 1:3
N = 0; D = 0;
for (n in 0:50) {
    N = N + dnbinom1(n, nn, vf) * (1 + (1 - u) / (n + 1) + (4 + 4 * (1 - u) / (n + 2)) * u + (8 + 9 * (1 - u) / (n + 3)) * u^2 - 4 * u^3 - 9 * u^4)
    D = D + dnbinom1(n, nn, vf) * (1 - u) * (1/(n + 1) + 4 * u/(n + 2) + 9 * u^2 / (n + 3))
}
N/D

N = ((-1 + u)*(6*nn^3*(1 + 5*u + 13*u^2 + 9*u^3 + 
          (1 + u*(2 + 3*u))*(k/(k + nn))^k) + 
       k*nn^2*((6 + 54*u^2 - 4*u*(-6 + nn) - 5*nn)*(k/(k + nn))^k - 
          (1 + u*(4 + 9*u))*(6 + 11*(1 + u)*nn)) + 
       k^3*(-18*u^2 + 2*u*(2 + 9*u)*nn - (1 + u*(4 + 9*u))*nn^2 - 
          (1 + u)*(1 + u*(4 + 9*u))*nn^3 + 
          (k/(k + nn))^k*(18*u^2 - 4*u*nn + nn^2)) + 
       k^2*nn*(-12*u - 54*u^2 + 5*nn + 20*u*nn + 45*u^2*nn + 
          6*nn^2 + 30*u*nn^2 + 78*u^2*nn^2 + 
          54*u^3*nn^2 + 
          (k/(k + nn))^k*(6*u*(2 + 9*u) - (5 + 8*u)*nn + nn^2))))/
   ((-3 + k)*(-2 + k)*(-1 + k)*nn^3)
D = ((-1 + u)*((k/(k + nn))^k*(k + nn)*
        (18*k^2*u^2 - 4*k*(-3 + k - 9*u)*u*nn + 
          (6 + k^2 + 6*u*(2 + 3*u) - k*(5 + 4*u))*nn^2) + 
       k*(-18*k^2*u^2 + 2*(-3 + k)*k*u*(2 + 9*u)*nn - 
          (-3 + k)*(-2 + k)*(1 + u*(4 + 9*u))*nn^2)))/
   ((-3 + k)*(-2 + k)*(-1 + k)*nn^3)
N/D

# Negative binomial, m = 1:3, fully simplified result
Post = (6*nn^3*(1 + 5*u + 13*u^2 + 9*u^3 + 
        (1 + u*(2 + 3*u))*(k/(k + nn))^k) + 
     k*nn^2*((6 + 54*u^2 - 4*u*(-6 + nn) - 5*nn)*(k/(k + nn))^k - 
        (1 + u*(4 + 9*u))*(6 + 11*(1 + u)*nn)) + 
     k^3*(-18*u^2 + 2*u*(2 + 9*u)*nn - (1 + u*(4 + 9*u))*nn^2 - 
        (1 + u)*(1 + u*(4 + 9*u))*nn^3 + 
        (k/(k + nn))^k*(18*u^2 - 4*u*nn + nn^2)) + 
     k^2*nn*(-12*u - 54*u^2 + 5*nn + 20*u*nn + 45*u^2*nn + 
        6*nn^2 + 30*u*nn^2 + 78*u^2*nn^2 + 
        54*u^3*nn^2 + (k/(k + nn))^k*
         (6*u*(2 + 9*u) - (5 + 8*u)*nn + nn^2)))/
   ((k/(k + nn))^k*(k + nn)*(18*k^2*u^2 - 4*k*(-3 + k - 9*u)*u*nn + 
        (6 + k^2 + 6*u*(2 + 3*u) - k*(5 + 4*u))*nn^2) + 
     k*(-18*k^2*u^2 + 2*(-3 + k)*k*u*(2 + 9*u)*nn - 
        (-3 + k)*(-2 + k)*(1 + u*(4 + 9*u))*nn^2))

# Negative binomial, m = 1:2
N = 0; D = 0;
for (n in 0:50) {
    N = N + dnbinom1(n, nn, vf) * ((2 + n)/(1 + n) + ((10 + 15*n + 4*n^2)*u)/((1 + n)*(2 + n)) - ((6 + n)*u^2)/(2 + n) - 4*u^3)
    D = D + dnbinom1(n, nn, vf) * ((1 - u)*(1/(1 + n) + (4*u)/(2 + n)))
}
N/D

N = ((-1 + u)*(2*nn^2*(-1 - 5*u - 4*u^2 - (1 + 2*u)*(k/(k + nn))^k) + 
       k^2*(4*u - nn - 4*u*nn - (1 + u)*(1 + 4*u)*nn^2 + 
          (k/(k + nn))^k*(-4*u + nn)) + 
       k*nn*((k/(k + nn))^k*(-2 - 8*u + nn) + (1 + 4*u)*(2 + 3*(1 + u)*nn))))/
   ((-2 + k)*(-1 + k)*nn^2)
D = ((-1 + u)*((k/(k + nn))^k*(k + nn)*(-4*k*u + (-2 + k - 4*u)*nn) + 
       k*(4*k*u - (-2 + k)*(1 + 4*u)*nn)))/((-2 + k)*(-1 + k)*nn^2)
N/D

# Negative binomial, m = 1:2, fully simplified result
Post = (2*nn^2*(-1 - 5*u - 4*u^2 - (1 + 2*u)*(k/(k + nn))^k) + 
     k^2*(4*u - nn - 4*u*nn - (1 + u)*(1 + 4*u)*nn^2 + 
        (k/(k + nn))^k*(-4*u + nn)) + 
     k*nn*((k/(k + nn))^k*(-2 - 8*u + nn) + (1 + 4*u)*(2 + 3*(1 + u)*nn)))/
   ((k/(k + nn))^k*(k + nn)*(-4*k*u + (-2 + k - 4*u)*nn) + 
     k*(4*k*u - (-2 + k)*(1 + 4*u)*nn))

# Negative binomial, m = 1:1
N = 0; D = 0;
for (n in 0:50) {
    N = N + dnbinom1(n, nn, vf) * ((1 - u)*(2 + n + u + n*u))/(1 + n)
    D = D + dnbinom1(n, nn, vf) * (1 - u)/(1 + n)
}
N/D

N = ((-1 + u)*(-k + nn - k*nn + u*nn - k*u*nn + (k/(k + nn))^k*(k + nn)))/((-1 + k)*nn)
D = ((-1 + u)*(-k + (k/(k + nn))^k*(k + nn)))/((-1 + k)*nn)
N/D

# Negative binomial, m = 1:1, fully simplified result
Post = 1 - ((-1 + k)*(1 + u)*nn)/(-k + (k/(k + nn))^k*(k + nn))



### Before duplication event ###

# Note this is written differently than pre duplication event for private, but
# it is the same conceptually and has the same answer (though they were 
# simplified differently).
N = 0; D = 0;
for (n in 0:50) {
  for (m in 1:30) {
    for (k in 0:35) {
      for (l in 0:k) {
        N = N + dpois(n, nn) * dquasi1(m, u) * dpois(k, (n+m)*u) * dbinom(l, k, m / (n + m)) * (m + l) * (n + m + k)
        D = D + dpois(n, nn) * dquasi1(m, u) * dpois(k, (n+m)*u) * dbinom(l, k, m / (n + m)) * (m + l)
      }
    }
  }
}
N/D

N = 0; D = 0;
for (n in 0:50) {
  for (m in 1:30) {
    for (k in 0:35) {
      N = N + dpois(n, nn) * dquasi1(m, u) * dpois(k, (n+m)*u) * (m + k * m / (n + m)) * (n + m + k)
      D = D + dpois(n, nn) * dquasi1(m, u) * dpois(k, (n+m)*u) * (m + k * m / (n + m))
    }
  }
}
N/D

N = 0; D = 0;
for (n in 0:50) {
  for (m in 1:30) {
    N = N + dpois(n, nn) * dquasi1(m, u) * (m^2 * (u^2 + 2*u + 1) + m * (n * (u^2 + 2*u + 1) + u))
    D = D + dpois(n, nn) * dquasi1(m, u) * (m + u * m)
  }
}
N/D

N = 0; D = 0;
for (n in 0:50) {
    N = N + dpois(n, nn) * ((1 + u) * (u^2 + 2*u + 1) / (1 - u)^2 + (n * (u^2 + 2*u + 1) + u) / (1 - u))
    D = D + dpois(n, nn) * (1 + u) / (1 - u)
}
N/D

# Fully simplified result (m = 1:infinity)
Pre = (1 + u)^2 / (1 - u) + nn * (1 + u) + u / (1 + u)

# So, above is solvable over m = 1:infinity.
# However, this may not work well with above approximations of m = 1:1, 1:2, 1:3
# So will make those approximations here as well.
# Note that all the expressions for "Pre" work equally well if the wild-type TEs
# are differently distributed, since the expression only relies on the first moment
# of the wild-type distribution.
# Taking m = 1:3
N = 0; D = 0;
for (n in 0:50) {
  for (m in 1:3) {
    N = N + dpois(n, nn) * dquasi1(m, u) * (m^2 * (u^2 + 2*u + 1) + m * (n * (u^2 + 2*u + 1) + u))
    D = D + dpois(n, nn) * dquasi1(m, u) * (m + u * m)
  }
}
N/D

N = 0; D = 0;
for (n in 0:50) {
    N = N + dpois(n, nn) * ( 1 + n + 3*(2 + n)*u + (13 + 4*n)*u^2 + 5*u^3 - (16 + 5*n)*u^4 - 3*(3 + n)*u^5 )
    D = D + dpois(n, nn) * ( 1 + 2*u + 2*u^2 - 2*u^3 - 3*u^4 )
}
N/D

N = ( 1 + nn + 3*(2 + nn)*u + (13 + 4*nn)*u^2 + 5*u^3 - (16 + 5*nn)*u^4 - 3*(3 + nn)*u^5 )
D = ( 1 + 2*u + 2*u^2 - 2*u^3 - 3*u^4 )
Pre = nn + 10/3 + (3 + nn) * u - 1/(1 + u) - (4 + 8*u) / (3 + 6*u + 9*u^2)

# Taking m = 1:2
N = 0; D = 0;
for (n in 0:50) {
  for (m in 1:2) {
    N = N + dpois(n, nn) * dquasi1(m, u) * (m^2 * (u^2 + 2*u + 1) + m * (n * (u^2 + 2*u + 1) + u))
    D = D + dpois(n, nn) * dquasi1(m, u) * (m + u * m)
  }
}
N/D

N = 0; D = 0;
for (n in 0:50) {
    N = N + dpois(n, nn) * (1 - u) * (1 + 7*u + 11*u^2 + 4*u^3 + n * (1 + u)^2 * (1 + 2*u))
    D = D + dpois(n, nn) * (1 + 2*u - u^2 - 2*u^3)
}
N/D

Pre = nn * (1 + u) + (1 + u * (1 + u) * (7 + 4*u)) / ((1 + u) * (1 + 2*u))

# Taking m = 1:1
N = 0; D = 0;
for (n in 0:50) {
  for (m in 1:1) {
    N = N + dpois(n, nn) * dquasi1(m, u) * (m^2 * (u^2 + 2*u + 1) + m * (n * (u^2 + 2*u + 1) + u))
    D = D + dpois(n, nn) * dquasi1(m, u) * (m + u * m)
  }
}
N/D

N = 0; D = 0;
for (n in 0:50) {
    N = N + dpois(n, nn) * (1 - u) * (1 + 3*u + u^2 + n * (1 + u)^2)
    D = D + dpois(n, nn) * (1 - u^2)
}
N/D

Pre = 2 + nn + u + nn * u - 1/(1 + u)

### Summary ###
# Post3, Post2, Post1 are approximations to PostNum with 3, 2, 1 terms in
# quasigeometric distribution.
Post3 = (exp(nn) * (18*u^2*(-nn+1) + nn^2 + 4*u*(nn^2-nn) + 9*u^2*nn^2 + nn^3 + 5*u*nn^3 + 13*u^2*nn^3 + 9*u^3*nn^3) - 18*u^2 + 4*u*nn - nn^2) /
    (exp(nn) * (18*u^2*(-nn+1) + nn^2 + 4*u*(nn^2-nn) + 9*u^2*nn^2) - 18*u^2 + 4*u*nn - nn^2)
Post2 = (exp(nn)*(nn*(1 + 4*u)*(1 + nn + nn*u) - 4*u) - nn + 4*u) / (exp(nn)*(nn + 4*(nn - 1)*u) - nn + 4*u)
Post1 = (exp(nn) * (nn + u*nn + 1) - 1) / (exp(nn) - 1)

# Pre3, Pre2, Pre1 approximations to Pre0 with 3, 2, 1 terms in qg distribution.
# Pre0 is full-precision, but finite-term versions of Pre may be better when 
# paired with finite-term versions of Post.
Pre0 = (1 + u)^2 / (1 - u) + nn * (1 + u) + u / (1 + u)
Pre3 = nn + 10/3 + (3 + nn) * u - 1/(1 + u) - (4 + 8*u) / (3 + 6*u + 9*u^2)
Pre2 = nn * (1 + u) + (1 + u * (1 + u) * (7 + 4*u)) / ((1 + u) * (1 + 2*u))
Pre1 = 2 + nn + u + nn * u - 1/(1 + u)

# Comparing results
PostNum - Pre0 # "true" result
Post1 - Pre1
Post2 - Pre2
Post3 - Pre3

Post1 - Pre0 # Probably not worth considering
Post2 - Pre0 # Probably not worth considering
Post3 - Pre0 # Considered, but doesn't work well for high-u

# Three versions worth considering from the above derivation.
dN_11 = nn * (1 + u) / (exp(nn) - 1) - u * (2 + u) / (1 + u)
dN_22 = (exp(nn) * (nn * (1 + 4*u) * (1 + nn + nn*u) - 4*u) - nn + 4*u) / (exp(nn) * (nn + 4 * (nn - 1) * u) - nn + 4*u) - nn * (1 + u) - (1 + u * (1 + u) * (7 + 4*u)) / ((1 + u) * (1 + 2 * u))
dN_33 = (exp(nn) * (18*u^2*(-nn+1) + nn^2 + 4*u*(nn^2-nn) + 9*u^2*nn^2 + nn^3 + 5*u*nn^3 + 13*u^2*nn^3 + 9*u^3*nn^3) - 18*u^2 + 4*u*nn - nn^2) /
    (exp(nn) * (18*u^2*(-nn+1) + nn^2 + 4*u*(nn^2-nn) + 9*u^2*nn^2) - 18*u^2 + 4*u*nn - nn^2) - (nn + 10/3 + (3 + nn) * u - 1/(1 + u) - (4 + 8*u) / (3 + 6*u + 9*u^2))

# Alternative derivation yielding a simpler expression
# From notebook p.66, 6 Sep 2023 
# N_B excess for public - from notebook
ds[mode == "pub" & tag == "new" & stage %in% 1:2 & g %in% 99:100, 
    .(dN_ds = N_B[4] - N_B[2], guess = (1 + n[1] - exp(-n[1])) / (1 - exp(-n[1])) - N_B[1]), 
    by = .(sim_id, u = uA)]

# Can substitute n + 1 + u/(1+u) for N_B (not proven) - this is a guess.
# Note also that a better guess for N_B is n + Vn/n as shown above -- this must somehow be provable
ds[mode == "pub" & tag == "new" & stage %in% 1:2 & g %in% 99:100, 
    .(dN_ds = N_B[4] - N_B[2], N_B_pre = N_B[1], 
        N_B_pre_guess = n[1] + 1 + uA[1] / (1 + uA[1]),
        N_B_pre_guess2 = n[1] + Vn[1]/n[1],
        guess = (1 + n[1] - exp(-n[1])) / (1 - exp(-n[1])) - (n[1] + 1 + uA[1] / (1 + uA[1]))), 
    by = .(sim_id, u = uA)]

# With the substitution n + 1 + 2u, the expression simplifies
# (1 + n - exp(-n)) / (1 - exp(-n)) - (n + 1 + u/(1+u))
# (1 + n - exp(-n) - (n + 1)*(1 - exp(-n))) / (1 - exp(-n)) - u/(1+u)
# (1 + n - exp(-n) - (n + 1 - nexp(-n) - exp(-n))) / (1 - exp(-n)) - u/(1+u)
# (1 + n - exp(-n) - n - 1 + nexp(-n) + exp(-n))) / (1 - exp(-n)) - u/(1+u)
# nexp(-n) / (1 - exp(-n)) - u/(1+u)
# n / (exp(n) - 1) - u/(1+u)
ds[mode == "pub" & tag == "new" & stage %in% 1:2 & g %in% 99:100, 
    .(dN_ds = N_B[4] - N_B[2], guess = n[1] / (exp(n[1]) - 1) - uA[1] / (1 + uA[1])), 
    by = .(sim_id, u = uA)]

# But following from above let's also try substituting n + Vn/n
ds[mode == "pub" & tag == "new" & stage %in% 1:2 & g %in% 99:100, 
    .(dN_ds = N_B[4] - N_B[2], guess = (1 + n[1] - exp(-n[1])) / (1 - exp(-n[1])) - (n[1] + Vn[1]/n[1])), 
    by = .(sim_id, u = uA)]

# Simplifying...
# (1 + n - exp(-n)) / (1 - exp(-n)) - n - Vn/n
# (1 + n - exp(-n) - n(1 - exp(-n))) / (1 - exp(-n)) - Vn/n
# (1 + n - exp(-n) - n + nexp(-n)) / (1 - exp(-n)) - Vn/n
# (1 - exp(-n) + nexp(-n)) / (1 - exp(-n)) - Vn/n
# nexp(-n) / (1 - exp(-n)) + 1 - Vn/n
# n / (exp(n) - 1) + 1 - Vn/n
ds[mode == "pub" & tag == "new" & stage %in% 1:2 & g %in% 99:100, 
    .(dN_ds = N_B[4] - N_B[2], guess = n[1] / (exp(n[1]) - 1) + 1 - Vn[1]/n[1]), 
    by = .(sim_id, u = uA)]

# Choosing from among these

# x = ns[mode == "pub" & g %in% c(GC - 1, GC) & tag == "new",
#     .(dN = N_B[4] - N_B[2]), by = .(sim_id, u, cost)]
# x = merge(x, ds[mode == "pub" & g == GC - 1 & stage == 1 & tag == "new", .(sim_id, nn = n, vf = Vn/n)], by = "sim_id")
x = ds[mode == "pub" & g %in% c(GC - 1, GC) & tag == "new",
    .(nn = n[1], dN = N_B[5] - N_B[2], vf = Vn[1]/n[1]), by = .(sim_id, u = uA, cost)]
x[, k := nn/(vf - 1)]

x[, dN_11 := nn * (1 + u) / (exp(nn) - 1) - u * (2 + u) / (1 + u)];
x[, dN_22 := (exp(nn) * (nn * (1 + 4*u) * (1 + nn + nn*u) - 4*u) - nn + 4*u) / (exp(nn) * (nn + 4 * (nn - 1) * u) - nn + 4*u) - nn * (1 + u) - (1 + u * (1 + u) * (7 + 4*u)) / ((1 + u) * (1 + 2 * u))];
x[, dN_33 := (exp(nn) * (18*u^2*(-nn+1) + nn^2 + 4*u*(nn^2-nn) + 9*u^2*nn^2 + nn^3 + 5*u*nn^3 + 13*u^2*nn^3 + 9*u^3*nn^3) - 18*u^2 + 4*u*nn - nn^2) /
    (exp(nn) * (18*u^2*(-nn+1) + nn^2 + 4*u*(nn^2-nn) + 9*u^2*nn^2) - 18*u^2 + 4*u*nn - nn^2) - (nn + 10/3 + (3 + nn) * u - 1/(1 + u) - (4 + 8*u) / (3 + 6*u + 9*u^2))];

x[, dN_11_nb := (1 - ((-1 + k)*(1 + u)*nn)/(-k + (k/(k + nn))^k*(k + nn))) - (2 + nn + u + nn * u - 1/(1 + u))];
x[, dN_22_nb := ((2*nn^2*(-1 - 5*u - 4*u^2 - (1 + 2*u)*(k/(k + nn))^k) + 
     k^2*(4*u - nn - 4*u*nn - (1 + u)*(1 + 4*u)*nn^2 + 
        (k/(k + nn))^k*(-4*u + nn)) + 
     k*nn*((k/(k + nn))^k*(-2 - 8*u + nn) + (1 + 4*u)*(2 + 3*(1 + u)*nn)))/
   ((k/(k + nn))^k*(k + nn)*(-4*k*u + (-2 + k - 4*u)*nn) + 
     k*(4*k*u - (-2 + k)*(1 + 4*u)*nn))) - (nn * (1 + u) + (1 + u * (1 + u) * (7 + 4*u)) / ((1 + u) * (1 + 2*u)))];
x[, dN_33_nb := ((6*nn^3*(1 + 5*u + 13*u^2 + 9*u^3 + 
        (1 + u*(2 + 3*u))*(k/(k + nn))^k) + 
     k*nn^2*((6 + 54*u^2 - 4*u*(-6 + nn) - 5*nn)*(k/(k + nn))^k - 
        (1 + u*(4 + 9*u))*(6 + 11*(1 + u)*nn)) + 
     k^3*(-18*u^2 + 2*u*(2 + 9*u)*nn - (1 + u*(4 + 9*u))*nn^2 - 
        (1 + u)*(1 + u*(4 + 9*u))*nn^3 + 
        (k/(k + nn))^k*(18*u^2 - 4*u*nn + nn^2)) + 
     k^2*nn*(-12*u - 54*u^2 + 5*nn + 20*u*nn + 45*u^2*nn + 
        6*nn^2 + 30*u*nn^2 + 78*u^2*nn^2 + 
        54*u^3*nn^2 + (k/(k + nn))^k*
         (6*u*(2 + 9*u) - (5 + 8*u)*nn + nn^2)))/
   ((k/(k + nn))^k*(k + nn)*(18*k^2*u^2 - 4*k*(-3 + k - 9*u)*u*nn + 
        (6 + k^2 + 6*u*(2 + 3*u) - k*(5 + 4*u))*nn^2) + 
     k*(-18*k^2*u^2 + 2*(-3 + k)*k*u*(2 + 9*u)*nn - 
        (-3 + k)*(-2 + k)*(1 + u*(4 + 9*u))*nn^2))) - (nn + 10/3 + (3 + nn) * u - 1/(1 + u) - (4 + 8*u) / (3 + 6*u + 9*u^2))];

x[, dN_3 := nn / (exp(nn) - 1) - u / (1 + u)]; # Works quite well for ds
x[, dN_guess_ns := nn / (exp(nn) - 1) - u * (1 - u) / (2 * (1 + u))]; # Works very well for ns

x[, dN_ord2 := -((1 + nn + nn*u*(3 + 4*u) + u*(6 + 13*u))/(1 + 2*u*(1 + u))) + (18*(-1 + exp(nn))*u^2 + exp(nn)*nn^3*(1 + 4*u + 8*u^2) - 2*nn*u*(2*(-1 + u) + exp(nn)*(2 + 7*u)) + nn^2*(-1 + u + exp(nn)*(1 + u*(3 + 5*u))))/(18*(-1 + exp(nn))*u^2 - 2*nn*u*(2*(-1 + u) + exp(nn)*(2 + 7*u)) + nn^2*(-1 + u + exp(nn)*(1 + u*(3 + 5*u))))]
x[, dN_ord2inf := 2 + 4/(-1 + u) + u + 1/(1 + u) - nn*(1 + u) + (18*(-1 + exp(nn))*u^2 + exp(nn)*nn^3*(1 + 4*u + 8*u^2) - 
   2*nn*u*(2*(-1 + u) + exp(nn)*(2 + 7*u)) + nn^2*(-1 + u + exp(nn)*(1 + u*(3 + 5*u))))/
 (18*(-1 + exp(nn))*u^2 - 2*nn*u*(2*(-1 + u) + exp(nn)*(2 + 7*u)) + nn^2*(-1 + u + exp(nn)*(1 + u*(3 + 5*u))))]
x[, dN_ord1 := (16*(-1 + exp(nn))*u^2 - 16*nn*u^2 + nn^2*(1 + 2*u - (3 + exp(nn))*u^2))/((1 + 2*u)*(nn*(-1 + u) + 4*u + exp(nn)*(nn - 4*u + 3*nn*u)))]

x[, dN_yetanother := nn^2/(nn*(-1 + u) + 4*u + exp(nn)*(nn - 4*u + 3*nn*u))]

# Check approximations visually
x = melt(x, measure.vars = names(x)[names(x) %like% "^dN_"]);
ggplot(x) +
    geom_abline(aes(intercept = 0, slope = 1)) +
    geom_line(aes(x = dN, y = value, colour = variable, group = cost, linetype = as.factor(cost))) +
    facet_wrap(~variable)

ggplot(x[variable %in% c("dN_22", "dN_33_nb", "dN_3", "dN_new")]) +
    geom_abline(aes(intercept = 0, slope = 1)) +
    geom_line(aes(x = dN, y = value, colour = variable, group = paste(variable, cost), linetype = as.factor(cost)))

ggplot(x[variable %in% c("dN_22", "dN_33_nb", "dN_3", "dN_new") & cost == 5e-04]) +
    geom_abline(aes(intercept = 0, slope = 1)) +
    geom_line(aes(x = dN, y = value, colour = variable, group = paste(variable, cost), linetype = as.factor(cost)))

ggplot(x[variable %in% c("dN_22", "dN_33_nb", "dN_3", "dN_new") & cost == 5e-03]) +
    geom_abline(aes(intercept = 0, slope = 1)) +
    geom_line(aes(x = dN, y = value, colour = variable, group = paste(variable, cost), linetype = as.factor(cost)))

ggplot(x[variable %in% c("dN_22", "dN_33_nb", "dN_3", "dN_new") & cost == 5e-02]) +
    geom_abline(aes(intercept = 0, slope = 1)) +
    geom_line(aes(x = dN, y = value, colour = variable, group = paste(variable, cost), linetype = as.factor(cost)))

################################################
# dN for cost - public and private replication #
################################################

# After duplication event
nn = 1.7
u = 0.01

N = 0; D = 0;
for (n in 0:50) {
  for (m in 1:10) {
    for (i in 0:25) {
      for (j in 0:5) {
        N = N + dpois(n, nn) * dquasi1(m, u) * m * dpois(i, n*u) * dpois(j, m*u) * (m + j) * (n + m + i + j + 1)
        D = D + dpois(n, nn) * dquasi1(m, u) * m * dpois(i, n*u) * dpois(j, m*u) * (m + j)
      }
    }
  }
}
paste(N/D)

N = 0; D = 0;
for (n in 0:50) {
  for (m in 1:10) {
    for (i in 0:25) {
      for (j in 0:5) {
        N = N + dpois(n, nn) * dquasi1(m, u) * m * dpois(i, n*u) * dpois(j, m*u) * (m * (n + m + i + 1) + j * (n + 2*m + i + 1) + j^2) 
        D = D + dpois(n, nn) * dquasi1(m, u) * m * dpois(i, n*u) * dpois(j, m*u) * (m + j)
      }
    }
  }
}
paste(N/D)

N = 0; D = 0;
for (n in 0:50) {
  for (m in 1:10) {
    for (i in 0:25) {
      N = N + dpois(n, nn) * dquasi1(m, u) * m * dpois(i, n*u) * (m * (n + m + i + 1) + m*u * (n + 2*m + i + 1) + m*u*(m*u+1)) 
      D = D + dpois(n, nn) * dquasi1(m, u) * m * dpois(i, n*u) * (m + m*u)
    }
  }
}
paste(N/D)

N = 0; D = 0;
for (n in 0:50) {
  for (m in 1:10) {
    N = N + dpois(n, nn) * dquasi1(m, u) * m * (m * (n + m + n*u + 1) + m*u * (n + 2*m + n*u + 1) + m*u*(m*u+1)) 
    D = D + dpois(n, nn) * dquasi1(m, u) * m * (m + m*u)
  }
}
paste(N/D)

N = 0; D = 0;
for (n in 0:50) {
  for (m in 1:10) {
    N = N + dpois(n, nn) * dquasi1(m, u) * (m^2 * (n * (1 + u)^2 + 1 + 2*u) + m^3*(1 + u)^2) 
    D = D + dpois(n, nn) * dquasi1(m, u) * m^2 * (1 + u)
  }
}
paste(N/D)

N = 0; D = 0;
for (n in 0:50) {
    N = N + dpois(n, nn) * ((1 + u)/(1 - u)^2 * (n * (1 + u)^2 + 1 + 2*u) + (1 + 4*u + u^2)/(1 - u)^3 * (1 + u)^2) 
    D = D + dpois(n, nn) * (1 + u)/(1 - u)^2 * (1 + u)
}
paste(N/D)

N = ((1 + u)/(1 - u)^2 * (nn * (1 + u)^2 + 1 + 2*u) + (1 + 4*u + u^2)/(1 - u)^3 * (1 + u)^2) 
D = (1 + u)/(1 - u)^2 * (1 + u)
paste(N/D)

# Simplified result:
Post = (nn * (1 - u^2) + (2 + 6*u + 3*u^2 + u^3) / (1 + u))/(1 - u)
paste(Post)

# Before duplication event
# Same as above, so reuse:
Pre = (1 + u)^2 / (1 - u) + nn * (1 + u) + u / (1 + u)

# Then dN = Post - Pre, which simplifies to:
(1 + u) / (1 - u)


# Check
x = ns[g %in% c(GC - 1, GC) & tag == "old" & stage == 2,
    .(dN.ns = N_B[2] - N_B[1]), by = .(sim_id, u, cost, mode)]
x = merge(x, 
    ds[g %in% c(GC - 1, GC) & tag == "old" & stage == 2,
        .(dN.ds = N_B[2] - N_B[1], vf = Vn[1]/n[1]), by = .(sim_id)],
    by = "sim_id")

# Note -- the dNs with "_ff" suffixes include "fudge factors" based on the 
# variance factor Vn/n to fit better. These fudge factors were fitted below.
x[, dN_qg := DD * (1 + u) / (1 - u)];
x[, dN_qg_ff := DD * (1 + u) / (1 - u) / (0.444 * vf^2 - 0.877 * vf + 1.433)];
x[, dN_emp := DD * (1 + 2*u) / (1 - u^2)];
x[, dN_emp_ff := DD * (1 + 2*u) / (1 - u^2) / (0.325 * vf^2 - 0.628 * vf + 1.303)];

# Check approximations visually
x = melt(x, measure.vars = names(x)[names(x) %like% "^dN_"]);

ggplot(x) +
    geom_abline(aes(intercept = 0, slope = 1)) +
    geom_line(aes(x = dN.ds - DD, y = value - DD, colour = variable, group = cost, linetype = as.factor(cost))) +
    facet_wrap(~variable) +
    scale_x_log10() + scale_y_log10()

ggplot(x) +
    geom_abline(aes(intercept = 0, slope = 1)) +
    geom_line(aes(x = dN.ns - DD, y = value - DD, colour = variable, group = cost, linetype = as.factor(cost))) +
    facet_wrap(~variable) +
    scale_x_log10() + scale_y_log10()

### Fitting of fudge factors (not used in paper)

# Quasigeometric case
# Get data
x = ds[tag == "old" & stage == 2 & g %in% c(GC - 1, GC), .(sim_id, u = uA, cost, N_B)];
x = x[, .(dN_meas = N_B[2] - N_B[1], dN_theor = DD * (1 + u[1]) / (1 - u[1])), by = .(sim_id, u, cost)]
x = merge(x,
    ds[tag == "old" & stage == 2 & g == GC - 1, .(sim_id, n, b, swn, sWN, sWN_B, Vn, Va, Vb, Vab)],
    by = "sim_id");

# Fit model
glm(I(dN_theor/dN_meas) ~ I(Vn/n) + I((Vn/n)^2), data = x)

# Plot to check
xx = seq(1, 2.4, 0.1)
ggplot(x) +
    geom_point(aes(Vn/n, dN_theor/dN_meas, colour = log10(cost))) +
    annotate("line", x = xx, y = 0.444 * xx^2 - 0.877 * xx + 1.433)

# "Empirical" (guessed) case
# Get data
x = ds[tag == "old" & stage == 2 & g %in% c(GC - 1, GC), .(sim_id, u = uA, cost, N_B)];
x = x[, .(dN_meas = N_B[2] - N_B[1], dN_theor = DD * (1 + 2 * u[1]) / (1 - u[1]^2)), by = .(sim_id, u, cost)]
x = merge(x,
    ds[tag == "old" & stage == 2 & g == GC - 1, .(sim_id, n, b, swn, sWN, sWN_B, Vn, Va, Vb, Vab)],
    by = "sim_id");

# Fit model
glm(I(dN_theor/dN_meas) ~ I(Vn/n) + I((Vn/n)^2), data = x)

# Plot to check
xx = seq(1, 2.5, 0.1)
ggplot(x) +
    geom_point(aes(Vn/n, dN_theor/dN_meas, colour = log10(cost))) +
    annotate("line", x = xx, y = 0.325 * xx^2 - 0.628 * xx + 1.303)



#############################################
# dN for cost more explicity based on delta #
#############################################

nn = 2.5
u = 0.5
d = 1e-4

N = 0; D = 0;
for (n in 0:50) {
  for (m in 1:30) {
    for (i in 0:50) {
      for (j in 0:30) {
          k = 0:m
          N = N + dpois(n, nn) * dquasi1(m, u) * dpois(i, n*u) * dpois(j, m*u) * sum(dbinom(k, m, d) * (m + j) * (n + m + i + j + k))
          D = D + dpois(n, nn) * dquasi1(m, u) * dpois(i, n*u) * dpois(j, m*u) * sum(dbinom(k, m, d) * (m + j))
      }
    }
  }
}
paste(N/D)

N = 0; D = 0;
for (n in 0:50) {
  for (m in 1:30) {
    for (i in 0:50) {
      for (j in 0:40) {
          N = N + dpois(n, nn) * dquasi1(m, u) * dpois(i, n*u) * dpois(j, m*u) * (m + j) * (n + m + i + j + m*d)
          D = D + dpois(n, nn) * dquasi1(m, u) * dpois(i, n*u) * dpois(j, m*u) * (m + j)
      }
    }
  }
}
paste(N/D)

N = 0; D = 0;
for (n in 0:50) {
    for (m in 1:30) {
        for (i in 0:50) {
            N = N + dpois(n, nn) * dquasi1(m, u) * dpois(i, n*u) * (m*u * (m*u + 1) + m*u * (n + 2*m + i + m*d) + m * (n + m + i + m*d))
            D = D + dpois(n, nn) * dquasi1(m, u) * dpois(i, n*u) * (m + m*u)
        }
    }
}
paste(N/D)

N = 0; D = 0;
for (n in 0:100) {
    for (m in 1:100) {
        N = N + dpois(n, nn) * dquasi1(m, u) * (m^2 * ((1 + u)^2 + d * (1 + u)) + m * (n * (1 + u)^2 + u))
        D = D + dpois(n, nn) * dquasi1(m, u) * (m * (1 + u))
    }
}
paste(N/D)

N = 0; D = 0;
for (n in 0:50) {
    N = N + dpois(n, nn) * ((1+u)/((1-u)^2) * ((1 + u)^2 + d * (1 + u)) + 1/(1-u) * (n * (1 + u)^2 + u))
    D = D + dpois(n, nn) * (1/(1-u) * (1 + u))
}
paste(N/D)


N = ((1+u)/((1-u)^2) * ((1 + u)^2 + d * (1 + u)) + 1/(1-u) * (nn * (1 + u)^2 + u))
D = (1/(1-u) * (1 + u))
paste(N/D)

N = ((1+u)/(1-u) * ((1 + u)^2 + d * (1 + u)) + (nn * (1 + u)^2 + u))
D = (1 + u)
paste(N/D)

N = nn * (1 - u^2) + (1 + u)^2 + d * (1 + u) + u * (1 - u) / (1 + u)
D = 1 - u
paste(N/D)

# Results in d(1+u)/(1-u)



################################
# dN in subsequent generations #
################################

# Generation t
u = 0.1
d = 1e-2
beta = -0.0123

# This is a relatively simple case
N = 0; D = 0;
for (s in 0:1) {
    N = N + dbinom(s, 1, d) * (1 + s) * (1 + s)
    D = D + dbinom(s, 1, d) * (1 + s)
}
paste(N/D)

# Final answer
paste((1 + 3*d) / (1 + d))

# Generation t + 1
N = 0; D = 0;
for (s in 0:1) {
    for (k in 0:(s+1)) {
        for (l in 0:100) {
            N = N + dbinom(s, 1, d) * exp(beta*s) * dbinom(k, s+1, 0.5) * dpois(l, k*u) * (k + l) * (k + l)
            D = D + dbinom(s, 1, d) * exp(beta*s) * dbinom(k, s+1, 0.5) * dpois(l, k*u) * (k + l)
        }
    }
}
paste(N/D)

N = 0; D = 0;
for (s in 0:1) {
    for (k in 0:(s+1)) {
        for (l in 0:100) {
            N = N + dbinom(s, 1, d) * exp(beta*s) * dbinom(k, s+1, 0.5) * dpois(l, k*u) * (k^2 + 2*k*l + l^2)
            D = D + dbinom(s, 1, d) * exp(beta*s) * dbinom(k, s+1, 0.5) * dpois(l, k*u) * (k + l)
        }
    }
}
paste(N/D)

N = 0; D = 0;
for (s in 0:1) {
    for (k in 0:(s+1)) {
        N = N + dbinom(s, 1, d) * exp(beta*s) * dbinom(k, s+1, 0.5) * (k^2 + 2*k^2*u + k*u*(k*u+1))
        D = D + dbinom(s, 1, d) * exp(beta*s) * dbinom(k, s+1, 0.5) * (k * (1 + u))
    }
}
paste(N/D)

N = 0; D = 0;
for (s in 0:1) {
    for (k in 0:(s+1)) {
        N = N + dbinom(s, 1, d) * exp(beta*s) * dbinom(k, s+1, 0.5) * (k^2 * (1 + 2*u + u^2) + k*u)
        D = D + dbinom(s, 1, d) * exp(beta*s) * dbinom(k, s+1, 0.5) * (k * (1 + u))
    }
}
paste(N/D)

N = 0; D = 0;
for (s in 0:1) {
    N = N + dbinom(s, 1, d) * exp(beta*s) * (((s+1)/4 + (s+1)^2/4) * (1 + u)^2 + (s+1)/2 * u)
    D = D + dbinom(s, 1, d) * exp(beta*s) * (s+1)/2 * (1 + u)
}
paste(N/D)

N = 0; D = 0;
for (s in 0:1) {
    N = N + dbinom(s, 1, d) * exp(beta*s) * ((s^2 + 3*s + 2) * (1 + u)^2 + 2*(s + 1) * u)
    D = D + dbinom(s, 1, d) * exp(beta*s) * 2*(s + 1) * (1 + u)
}
paste(N/D)

N = (1 - d) * (1 + 3*u + u^2) + d * exp(beta) * (3 + 8*u + 3*u^2)
D = (1 - d) * (1 + u) + d * exp(beta) * 2 * (1 + u)

