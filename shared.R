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