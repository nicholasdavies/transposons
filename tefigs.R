library(ggplot2)
library(data.table)
library(cowplot)
library(scales)

source("./load_ngd.R")
theme_set(theme_cowplot(font_size = 8) + theme(plot.title = element_text(size = 8, face = "plain")))

# colours
cDup = "#7caf5c"
cFit = "#59386c"
cTE1 = "#802417"
cTE2 = "#c06636"
cTE3 = "#ce9344"
cTE4 = "#e8b960"
cTE5 = "#508ea2"
cTE6 = "#17486f"
cAct = "#447861"
cPaf = "#b1a1cc";
cSuf = "#738e8e";


# Figure panel for TE dynamics
panel = function(dt, y, ymin, ymax, ylim, colour, title = NULL, xlab = NULL, ylab = NULL, xlim = c(NA, NA), tx = TRUE)
{
    if (!is.vector(dt)) {
        dt = list(dt)
    }
    
    plot = ggplot();
    
    for (i in seq_along(y)) {
        if (!is.na(ymin[i])) {
            plot = plot + 
                geom_ribbon(data = dt[[i]], aes_string(x = "g", ymin = ymin[i], ymax = ymax[i]), fill = colour[i], alpha = 0.4)
        }
        plot = plot +
            geom_line(data = dt[[i]], aes_string(x = "g", y = y[i]), colour = colour[i])
    }
    
    if (tx == FALSE) {
        plot = plot + theme(axis.text.y = element_blank())
    }
    
    plot +
        scale_y_continuous(limits = ylim) +
        scale_x_continuous(limits = xlim, expand = expansion(0.02)) +
        labs(x = xlab, y = ylab, title = title)
}

vline = function(x)
{
     geom_vline(aes(xintercept = x), linetype = "dashed", size = 0.25)
}

############
# FIGURE 1 #
############

ta = load_ngd("./Deterministic/Runs/2-Static-vs-Sequential-DNA/1-static.ngd")
tb = load_ngd("./Deterministic/Runs/2-Static-vs-Sequential-DNA/2-coexistence.ngd")
tc = load_ngd("./TE/Runs/1-Tragedy/1_mut_pub_u01.ngd")
td = load_ngd("./Deterministic/Runs/1-Static-vs-Sequential-RNA/1-static.ngd")
te = load_ngd("./Deterministic/Runs/1-Static-vs-Sequential-RNA/2-sequential.ngd")
tf = load_ngd("./TE/Runs/1-Tragedy/1_mut_pri_u01.ngd")

tb[, uu := (u * n + U * m) / (n + m)]
te[, uu := (u * n + U * m) / (n + m)]

# Swap around n and m in te for plotting
te[, N := n]
te[, N05 := n05]
te[, N95 := n95]

te[g %between% c(8.01, 12), c("n", "n05", "n95") := .(m, m05, m95)]
te[g %between% c(8.01, 12), c("m", "m05", "m95") := .(N, N05, N95)]

te[g %between% c(16.01, 20), c("n", "n05", "n95") := .(m, m05, m95)]
te[g %between% c(16.01, 20), c("m", "m05", "m95") := .(N, N05, N95)]

f2a = plot_grid(
    panel(ta, "n", "n05", "n95", c(0, 150), cTE1, ylab = "Copy number", title = "One type") +
        draw_label(quote(italic(a) == ".005"), x = 20, y = 125, hjust = 1, vjust = 1, size = 8, color = cTE1),
    panel(ta, "u", NA, NA, c(0, 0.021), cDup, ylab = "Duplication rate") +
        scale_y_continuous(limits = c(0, 0.021), breaks = c(0, 0.01, 0.02)),
    panel(ta, "w", "w05", "w95", c(0, 1), cFit, ylab = "Host fitness", xlab = " "),
    ncol = 1, align = "v", axis = "b", rel_heights = c(1.4, 1, 1))

f2b = plot_grid(
    panel(list(tb, tb), c("n", "m"), c("n05", "m05"), c("n95", "m95"), c(0, 150), c(cTE1, cTE2), title = "   Multiple types", tx = FALSE) +
        draw_label(quote(italic(a) == ".005"), x = 20, y = 125, hjust = 1, vjust = 1, size = 8, color = cTE1) +
        draw_label(quote(italic(a) == ".010"), x = 20, y =  40, hjust = 1, vjust = 1, size = 8, color = cTE2),
    panel(tb, "uu", NA, NA, c(0, 0.021), cDup, tx = FALSE) +
        scale_y_continuous(limits = c(0, 0.021), breaks = c(0, 0.01, 0.02)),
    panel(tb, "w", "w05", "w95", c(0, 1), cFit, xlab = "Generations (thousands)", tx = FALSE),
    ncol = 1, align = "v", axis = "b", rel_heights = c(1.4, 1, 1))

f2c = plot_grid(
    panel(tc, "n", "n05", "n95", c(0, 150), cTE1, title = "   Free mutation", tx = FALSE) + vline(4),
    panel(tc, "u", "u05", "u95", c(0, 0.021), cDup, tx = FALSE) + vline(4) +
        scale_y_continuous(limits = c(0, 0.021), breaks = c(0, 0.01, 0.02)),
    panel(tc, "w", "w05", "w95", c(0, 1), cFit, xlab = " ", tx = FALSE) + vline(4),
    ncol = 1, align = "v", axis = "b", rel_heights = c(1.4, 1, 1))

f2abc = cowplot::plot_grid(f2a, NULL, f2b, NULL, f2c, NULL, nrow = 1, labels = c("E", "", "F", "", "G", ""), label_size = 8, rel_widths = c(1.3, 0.05, 1, 0.05, 1, 0.15))


f2d = plot_grid(
    panel(td, "n", "n05", "n95", c(0, 600), cTE1, ylab = "Copy number", title = "One type") +
        draw_label(quote(italic(a) == ".005"), x = 20, y = 200, hjust = 1, vjust = 1, size = 8, color = cTE1),
    panel(td, "u", NA, NA, c(0, 0.05), cDup, ylab = "Duplication rate"),
    panel(td, "w", "w05", "w95", c(0, 1), cFit, ylab = "Host fitness", xlab = " "),
    ncol = 1, align = "v", axis = "b", rel_heights = c(1.4, 1, 1))

f2e = plot_grid(
    panel(list(te[g %between% c(0, 8)], te[g %between% c(4, 12)], te[g %between% c(8, 16)], te[g %between% c(12, 20)], te[g %between% c(16, 20)]), 
        c("n", "m", "n", "m", "n"), 
        c("n05", "m05", "n05", "m05", "n05"), 
        c("n95", "m95", "n95", "m95", "n95"), 
        c(0, 600), c(cTE1, cTE2, cTE3, cTE4, cTE5), title = "   Multiple types", tx = FALSE) +
        draw_label(quote(italic(a) == ".010"), x =  9, y = 300, hjust = 1, vjust = 1, size = 8, color = cTE2) +
        draw_label(quote(italic(a) == ".015"), x = 13, y = 400, hjust = 1, vjust = 1, size = 8, color = cTE3) +
        draw_label(quote(italic(a) == ".020"), x = 17, y = 500, hjust = 1, vjust = 1, size = 8, color = cTE4) +
        draw_label(quote(italic(a) == ".025"), x = 20, y = 600, hjust = 1, vjust = 1, size = 8, color = cTE5),
    panel(te, "uu", NA, NA, c(0, 0.05), cDup, tx = FALSE),
    panel(te, "w", "w05", "w95", c(0, 1), cFit, xlab = "Generations (thousands)", tx = FALSE),
    ncol = 1, align = "v", axis = "b", rel_heights = c(1.4, 1, 1))

f2f = plot_grid(
    panel(tf, "n", "n05", "n95", c(0, 600), cTE1, xlim = c(0, 20), title = "   Free mutation", tx = FALSE) + vline(4),
    panel(tf, "u", "u05", "u95", c(0, 0.05), cDup, xlim = c(0, 20), tx = FALSE) + vline(4),
    panel(tf, "w", "w05", "w95", c(0, 1), cFit, xlab = " ", xlim = c(0, 20), tx = FALSE) + vline(4),
    ncol = 1, align = "v", axis = "b", rel_heights = c(1.4, 1, 1))

f2def = cowplot::plot_grid(NULL, f2d, NULL, f2e, NULL, f2f, nrow = 1, labels = c("", "H", "", "I", "", "J"), label_size = 8, rel_widths = c(0.15, 1.3, 0.05, 1, 0.05, 1))

title1 = ggdraw() + draw_label("Public replication", size = 8, fontface = "bold")
title2 = ggdraw() + draw_label("Private replication", size = 8, fontface = "bold")


f2 = cowplot::plot_grid(title1, title2, f2abc, f2def, nrow = 2, ncol = 2, rel_heights = c(1, 20))

ggsave("./Figures/1-tragedy-plots.pdf", f2, width = 20, height = 8, units = "cm", useDingbats = FALSE)
ggsave("./Figures/1-tragedy-plots.png", f2, width = 20, height = 8, units = "cm")


############
# FIGURE 3 #
############

# Load files
s1 = load_ngd("./TE/Runs/5-Parasite/parasite_com_phi090_len05_u01.ngd");
s2 = load_ngd("./TE/Runs/5-Parasite/parasite_com_phi100_len05_u01.ngd");
s3 = load_ngd("./TE/Runs/3-Suppression/csupp_GN_co7_k100.ngd");
s4 = load_ngd("./TE/Runs/3-Suppression/csupp_GP_co7_k100.ngd");

# Process files
s1[, n0 := nh0 + ng0]
s1[, n005 := nh005 + ng005]
s1[, n095 := nh095 + ng095]
s1[, n1 := nh1 + ng1]
s1[, n105 := nh105 + ng105]
s1[, n195 := nh195 + ng195]
s2[, n0 := nh0 + ng0]
s2[, n005 := nh005 + ng005]
s2[, n095 := nh095 + ng095]
s2[, n1 := nh1 + ng1]
s2[, n105 := nh105 + ng105]
s2[, n195 := nh195 + ng195]

# Parasitised fraction
s1[, p :=   ifelse(n0 == 0, NA_real_, 0.9 * a1 / (100 + a1))]
s1[, p05 := ifelse(n0 == 0, NA_real_, 0.9 * a105 / (100 + a105))]
s1[, p95 := ifelse(n0 == 0, NA_real_, 0.9 * a195 / (100 + a195))]

s2[, p :=   ifelse(n0 == 0, NA_real_, 1.0 * a1 / (100 + a1))]
s2[, p05 := ifelse(n0 == 0, NA_real_, 1.0 * a105 / (100 + a105))]
s2[, p95 := ifelse(n0 == 0, NA_real_, 1.0 * a195 / (100 + a195))]

# Suppressed fraction
s3[, s := S / (10 + S)]
s3[, s05 := S05 / (10 + S05)]
s3[, s95 := S95 / (10 + S95)]

s4[, s := S / (10 + S)]
s4[, s05 := S05 / (10 + S05)]
s4[, s95 := S95 / (10 + S95)]

# For plotting
s3[, U05 := pmax(U05, 1e-3)]
s4[, U05 := pmax(U05, 1e-3)]



f4a = plot_grid(
    panel(list(s1, s1), c("n0", "n1"), c("n005", "n105"), c("n095", "n195"), 
        c(0, 1500), c(cTE1, cTE5), ylab = "Copy\nnumber", title = quote(paste("Parasitism, ", italic(p), " = 0.9")), xlim = c(0, 40)) +
        draw_label("Replication-\ncompetent\nelement", x = 19, y = 50, hjust = 0, vjust = 0, size = 7, color = cTE1) +
        draw_label("Parasitic\nelement", x = 19, y = 1000, hjust = 0, vjust = 0, size = 7, color = cTE5),
    panel(s1, "u0", "u005", "u095", c(0, 3), cDup, ylab = "Duplication\nrate", xlim = c(0, 40)),
    panel(s1, "p", "p05", "p95", c(0, 1), cPaf, ylab = "Parasitized\nfraction", xlim = c(0, 40)),
    panel(s1, "w", "w05", "w95", c(0, 1), cFit, ylab = "Host\nfitness", xlab = " ", xlim = c(0, 40)),
    ncol = 1, align = "v", axis = "b", rel_heights = c(1.5, 1, 1, 1.2))

f4b = plot_grid(
    panel(list(s2, s2), c("n0", "n1"), c("n005", "n105"), c("n095", "n195"), 
        c(0, 525), c(cTE1, cTE5), ylab = NULL, title = quote(paste("Parasitism, ", italic(p), " = 1.0")), xlim = c(0, 40)),
    panel(s2, "u0", "u005", "u095", c(0, 3), cDup, ylab = NULL, xlim = c(0, 40)),
    panel(s2, "p", "p05", "p95", c(0, 1), cPaf, ylab = NULL, xlim = c(0, 40)),
    panel(s2, "w", "w05", "w95", c(0, 1), cFit, ylab = NULL, xlab = " ", xlim = c(0, 40)),
    ncol = 1, align = "v", axis = "b", rel_heights = c(1.5, 1, 1, 1.2))

f4c = plot_grid(
    panel(s3, "n", "n05", "n95",
        c(0, 2000), cTE1, ylab = "Copy\nnumber", title = "Transcriptional", xlim = c(0, 40)),
    panel(list(s3, s3), c("U", "u"), c("U05", "u05"), c("U95", "u95"), c(0, 0.265), c(cAct, cDup), ylab = "Duplication\nrate", xlim = c(0, 40)) + 
        scale_y_log10(limits = c(1e-3, 1e4), breaks = c(1e-2, 1, 1e2, 1e4), labels = trans_format("log10", math_format(10^.x))),
    panel(s3, "s", "s05", "s95", c(0, 1), cSuf, ylab = "Suppressed\nfraction", xlim = c(0, 40)),
    panel(s3, "w", "w05", "w95", c(0, 1), cFit, ylab = "Host\nfitness", xlab = " ", xlim = c(0, 40)),
    ncol = 1, align = "v", axis = "b", rel_heights = c(1.5, 1, 1, 1.2)) +
    draw_label("suppression", x = 0.58, y = 0.938, hjust = 0.5, vjust = 0.5, size = 8, color = "#000000")

f4d = plot_grid(
    panel(s4, "n", "n05", "n95",
        c(0, 1000), cTE1, ylab = NULL, title = "Post-transcriptional", xlim = c(0, 40)),
    panel(list(s4, s4), c("U", "u"), c("U05", "u05"), c("U95", "u95"), c(0, 0.265), c(cAct, cDup), ylab = NULL, xlim = c(0, 40)) + 
        scale_y_log10(limits = c(1e-3, 1e4), breaks = c(1e-2, 1, 1e2, 1e4), labels = trans_format("log10", math_format(10^.x))),
    panel(s4, "s", "s05", "s95", c(0, 1), cSuf, ylab = NULL, xlim = c(0, 40)),
    panel(s4, "w", "w05", "w95", c(0, 1), cFit, ylab = NULL, xlab = " ", xlim = c(0, 40)),
    ncol = 1, align = "v", axis = "b", rel_heights = c(1.5, 1, 1, 1.2)) +
    draw_label("suppression", x = 0.58, y = 0.938, hjust = 0.5, vjust = 0.5, size = 8, color = "#000000")

f4 = plot_grid(f4a, NULL, f4b, NULL, f4c, NULL, f4d, 
    nrow = 1, labels = c("B", "", "C", "", "E", "", "F"), label_size = 8, rel_widths = c(1.18, 0.073, 1, 0.291, 1.18, 0.073, 1)) + 
    draw_label("Generations (thousands)", x = 0.293, y = 0.03, hjust = 0.5, vjust = 0.5, size = 8, color = "#000000") +
    draw_label("Generations (thousands)", x = 0.793, y = 0.03, hjust = 0.5, vjust = 0.5, size = 8, color = "#000000")

ggsave("./Figures/3-parasupp-plots.pdf", f4, width = 16.5, height = 8.5, units = "cm", useDingbats = FALSE)
ggsave("./Figures/3-parasupp-plots.png", f4, width = 16.5, height = 8.5, units = "cm")





######################################
# SI - PARASITE SENSITIVITY ANALYSIS #
######################################

load_parasites = function(len, u)
{
    phis = c("050", "090", "091", "092", "093", "094", "095", "096", "097", "098", "099", "100");
    filenames = paste0("./TE/Runs/5-Parasite/parasite_com_phi", phis, "_len", len, "_u", u, ".ngd");
    s = list();
    for (i in seq_along(phis)) {
        s[[i]] = load_ngd(filenames[[i]], tag = as.numeric(phis[[i]])/100);
    }
    s = rbindlist(s);
    
    # Process
    s[, n0 := nh0 + ng0]
    s[, n005 := nh005 + ng005]
    s[, n095 := nh095 + ng095]
    s[, n1 := nh1 + ng1]
    s[, n105 := nh105 + ng105]
    s[, n195 := nh195 + ng195]
    
    # Parasitised fraction
    s[, p :=   ifelse(n0 == 0, NA_real_, tag * a1 / (100 + a1))]
    s[, p05 := ifelse(n0 == 0, NA_real_, tag * a105 / (100 + a105))]
    s[, p95 := ifelse(n0 == 0, NA_real_, tag * a195 / (100 + a195))]
    
    return (s[])
}

paraplot = function(s, title, xmax = NA)
{
    xadj = function()
        list(
            scale_x_continuous(limits = c(0, xmax), expand = expansion(0, 0)),
            theme(panel.spacing = unit(0.7, "lines"), plot.margin = margin(4, 6, 4, 3, unit = "pt"))
        )
    
    plot_grid(
        panel(list(s, s), c("n0", "n1"), c("n005", "n105"), c("n095", "n195"), 
            c(0, 2000), c(cTE1, cTE5), ylab = "Copy\nnumber", title = title) + facet_wrap(~paste0("italic(p)==", tag), nrow = 1, labeller = label_parsed) + xadj(),
        panel(s, "u0", "u005", "u095", c(0, 0.3), cDup, ylab = "Duplication\nrate") + facet_wrap(~paste0("p=", tag), nrow = 1) +
            theme(strip.background = element_blank(), strip.text = element_blank()) + xadj(),
        panel(s, "p", "p05", "p95", c(0, 1), cPaf, ylab = "Parasitized\nfraction") + facet_wrap(~paste0("p=", tag), nrow = 1) +
            theme(strip.background = element_blank(), strip.text = element_blank()) + xadj(),
        panel(s, "w", "w05", "w95", c(0, 1), cFit, ylab = "Host\nfitness", xlab = "Generations (thousands)") + facet_wrap(~paste0("p=", tag), nrow = 1) +
            theme(strip.background = element_blank(), strip.text = element_blank()) + xadj(),
        ncol = 1, align = "v", axis = "b", rel_heights = c(1.5, 1, 0.85, 1))
}


ss = load_parasites(len = "01", u = "01")   # 97 and up
ss = load_parasites(len = "01", u = "001")  # 97 and up
ss = load_parasites(len = "01", u = "0001") # 97 and up by 1e-06

ss = load_parasites(len = "05", u = "01")   # 96 and up
ss = load_parasites(len = "05", u = "001")  # 96 and up # THIS IS FIG. S2
ss = load_parasites(len = "05", u = "0001") # 96 and up

# All tragic
ss = load_parasites(len = "25", u = "01")   # nothing
ss = load_parasites(len = "25", u = "001")  # nothing
ss = load_parasites(len = "25", u = "0001") # nothing

s = load_parasites(len = "05", u = "001")
sp1 = paraplot(s[tag < 0.945], "", 50)
sp2 = paraplot(s[tag > 0.945], "", 50)

sp = cowplot::plot_grid(sp1, sp2, nrow = 2)
ggsave("./Figures/S-para.pdf", sp, width = 18, height = 20, units = "cm", useDingbats = FALSE)
ggsave("./Figures/S-para.png", sp, width = 18, height = 20, units = "cm")


#########################################
# SI - SUPPRESSION SENSITIVITY ANALYSIS #
#########################################

load_supps = function(type, cost, prefix = "csupp")
{
    ks = c("1", "10", "100", "1000", "10000");
    filenames = paste0("./TE/Runs/3-Suppression/", prefix, "_", type, "_co", cost, "_k", ks, ".ngd");
    su = list();
    for (i in seq_along(ks)) {
        if (file.exists(filenames[[i]])) {
            su[[i]] = load_ngd(filenames[[i]], tag = as.numeric(ks[[i]]));
        }
    }
    su = rbindlist(su);
    
    # Suppressed fraction
    su[, s   := S   / (10 + S)]
    su[, s05 := S05 / (10 + S05)]
    su[, s95 := S95 / (10 + S95)]
    
    return (su[])
}

supp_plot = function(s, title, nlim = c(0, 2000), alog = TRUE, alim = c(1e-3, 2e4), abreaks = c(1e-2, 1, 1e2, 1e4))
{
    xadj = function()
        list(
            scale_x_continuous(limits = c(0, 50), expand = expansion(c(0.02, 0), 0)),
            theme(panel.spacing = unit(0.7, "lines"), plot.margin = margin(4, 6, 4, 3, unit = "pt"))
        )
    
    plot_grid(
        panel(s, "n", "n05", "n95", 
            nlim, cTE1, ylab = "Copy\nnumber", title = title, xlim = c(0, 50)) + 
            facet_wrap(~paste0("italic(k)==",tag), nrow = 1, labeller = label_parsed) +
            theme(axis.text.x = element_blank()) + xadj(),
        panel(list(s, s), c("U", "u"), c("U05", "u05"), c("U95", "u95"), c(0, 1), c(cAct, cDup), ylab = "Duplication\nrate", xlim = c(0, 50)) + 
            facet_wrap(~tag, nrow = 1) +
            theme(strip.background = element_blank(), strip.text = element_blank(), axis.text.x = element_blank()) + xadj() +
            if (alog) {
                scale_y_log10(limits = alim, breaks = abreaks, labels = trans_format("log10", math_format(10^.x)))
            } else {
                scale_y_continuous(limits = alim, breaks = abreaks)
            },
        panel(s, "s", "s05", "s95", c(0, 1), cSuf, ylab = "Suppressed  \nfraction", xlim = c(0, 50)) + 
            facet_wrap(~tag, nrow = 1) +
            theme(strip.background = element_blank(), strip.text = element_blank(), axis.text.x = element_blank()) + xadj(),
        panel(s, "w", "w05", "w95", c(0, 1), cFit, ylab = "Host\nfitness", xlab = "Generations (thousands)", xlim = c(0, 50)) + 
            facet_wrap(~tag, nrow = 1) +
            theme(strip.background = element_blank(), strip.text = element_blank()) + xadj(),
        ncol = 1, align = "v", axis = "b", rel_heights = c(1.4, 0.9, 0.75, 1.1))
}

sN7 = load_supps("GN", "7")
sN6 = load_supps("GN", "6")
sN5 = load_supps("GN", "5")

sP7 = load_supps("GP", "7")
sP6 = load_supps("GP", "6")
sP5 = load_supps("GP", "5")

# For plotting
sN7[, U05 := pmax(U05, 1e-3)]
sN6[, U05 := pmax(U05, 1e-3)]
sN5[, U05 := pmax(U05, 1e-3)]

sP7[, U05 := pmax(U05, 1e-3)]
sP6[, U05 := pmax(U05, 1e-3)]
sP5[, U05 := pmax(U05, 1e-3)]

fsuppPr1 = plot_grid(
    supp_plot(sN7, expression("Transcriptional suppression, suppression cost = 10"^-7)),
    supp_plot(sN6, expression("Transcriptional suppression, suppression cost = 10"^-6)),
    supp_plot(sN5, expression("Transcriptional suppression, suppression cost = 10"^-5)),
    nrow = 3, ncol = 1, labels = LETTERS[1:3], label_size = 8
)

fsuppPr2 = plot_grid(
    supp_plot(sP7, expression("Post-transcriptional suppression, suppression cost = 10"^-7)),
    supp_plot(sP6, expression("Post-transcriptional suppression, suppression cost = 10"^-6)),
    supp_plot(sP5, expression("Post-transcriptional suppression, suppression cost = 10"^-5)),
    nrow = 3, ncol = 1, labels = LETTERS[1:3], label_size = 8
)

ggsave("./Figures/S-supp-private-1.pdf", fsuppPr1, width = 14, height = 22, units = "cm", useDingbats = FALSE)
ggsave("./Figures/S-supp-private-1.png", fsuppPr1, width = 14, height = 22, units = "cm")

ggsave("./Figures/S-supp-private-2.pdf", fsuppPr2, width = 14, height = 22, units = "cm", useDingbats = FALSE)
ggsave("./Figures/S-supp-private-2.png", fsuppPr2, width = 14, height = 22, units = "cm")


# Public elements suppression
sp0 = load_supps("0", "NA", "pubsupp")
# sp0 = load_supps("Null", "NA", "pubsupp")
spN = load_supps("GN", "5", "pubsupp")
spP = load_supps("GP", "5", "pubsupp")

fsuppPu = plot_grid(
    supp_plot(sp0, expression("No suppression"), nlim = c(0, 350), alog = FALSE, alim = c(0, 0.02), abreaks = waiver()),
    supp_plot(spN, expression("Transcriptional suppression, cost = 10"^-5), nlim = c(0, 350), alog = FALSE, alim = c(0, 0.02), abreaks = waiver()),
    supp_plot(spP, expression("Post-transcriptional suppression, cost = 10"^-5), nlim = c(0, 350), alog = FALSE, alim = c(0, 0.02), abreaks = waiver()),
    nrow = 3, ncol = 1, labels = LETTERS[1:6], label_size = 8
)

ggsave("./Figures/S-supp-public.pdf", fsuppPu, width = 14, height = 22, units = "cm", useDingbats = FALSE)
ggsave("./Figures/S-supp-public.png", fsuppPu, width = 14, height = 22, units = "cm")





##################
# DONG, ROTHKAMM #
##################

ddr = fread("./dongroth.csv")


pA = ggplot() + 
    geom_line(data = ddr[fig == "dong_line"], aes(x * 100, y), colour = "#0088ff", linewidth = 1) +
    geom_point(data = ddr[fig == "A"], aes(x * 100, y), colour = "red") +
    cowplot::theme_cowplot(font_size = 8) +
    labs(y = "Relative growth rate", x = "Amount of"~beta*"-galactosidases (%)") +
    coord_cartesian(xlim = c(0, 32), ylim = c(0, 1)) +
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30))


pB = ggplot() + 
    geom_line(data = ddr[fig == "dong_line"], aes(x * 100, y), colour = "#0088ff", linewidth = 1) +
    geom_point(data = ddr[fig == "B"], aes(x * 100, y), colour = "red") +
    cowplot::theme_cowplot(font_size = 8) +
    labs(y = "Relative growth rate", x = "Amount of"~Delta*"Ef-Tu (%)") +
    coord_cartesian(xlim = c(0, 32), ylim = c(0, 1)) +
    scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
    scale_x_continuous(breaks = c(0, 5, 10, 15, 20, 25, 30))

pC = ggplot() + 
    geom_line(data = ddr[fig == "rothkamm_line"], aes(x, y), colour = "#0088ff", linewidth = 1) +
    geom_point(data = ddr[fig == "C"], aes(x, y), colour = "red") +
    cowplot::theme_cowplot(font_size = 8) +
    labs(y = "Cell survival", x = "Double-strand breaks") +
    coord_cartesian(xlim = c(0, 300), ylim = c(0.03, 1)) +
    scale_y_log10(breaks = c(0.05, 0.1, 0.2, 0.5, 1.0)) +
    scale_x_continuous(breaks = c(0, 50, 100, 150, 200, 250, 300))

pp = cowplot::plot_grid(pA, pB, pC, ncol = 1, labels = LETTERS[1:3], label_size = 8)

ggsave("./Figures/S-dongroth.pdf", pp, width = 8, height = 16, units = "cm", useDingbats = FALSE)
ggsave("./Figures/S-dongroth.png", pp, width = 8, height = 16, units = "cm")

