# Helper functions for cost_benefit.R

# Load set of analyses from the deterministic simulation
load_analysis = function(mode, cost, dupl)
{
    # Load the analysis files with appropriate settings
    d = rbind(
        load_ngd(paste0("./Deterministic/Runs/4-Analysis/", mode, "_c", cost, "_u", dupl, "_up.ngd"), 
            thin = 1, gunit = 1, tag = "full"),
        load_ngd(paste0("./Deterministic/Runs/5-Analysis-Partial/", mode, "_c", cost, "_u", dupl, "_old.ngd"),
            thin = 1, gunit = 1, tag = "old"),
        load_ngd(paste0("./Deterministic/Runs/5-Analysis-Partial/", mode, "_c", cost, "_u", dupl, "_new.ngd"),
            thin = 1, gunit = 1, tag = "new")
    );
    
    # If file exists, load "par" too
    par_file = paste0("./Deterministic/Runs/5-Analysis-Partial/", mode, "_c", cost, "_u", dupl, "_par.ngd");
    if (file.exists(par_file)) {
        d = rbind(d, load_ngd(par_file, thin = 1, gunit = 1, tag = "par"));
    }
    
    # Decompose the 'g' variable into stage and generation 
    # (i.e. 1.3 becomes g = 1, stage = 3)
    d[, stage := round((g %% 1) * 10)];
    d[, g := floor(g)];
    
    # Add scenario information
    d[, mode := mode];
    d[, cost := as.numeric(paste0("0.", cost))];
    d[, dupl := as.numeric(paste0("0.", dupl))];
    
    return (d[])
}

# No-selection version of simulation, focusing on mutant lineage experience
# N     simulate NxN genotypes
# G     number of generations to run 
# GC    generation in which duplication event occurs
# u     duplication rate
# p     probability of additional duplication
# A_B   target experienced number of type-A TEs expected over type-B TEs at stage 3 (after recombination)
# mode  "pub" or "pri" for mode of replication
# focus "new", "old", or "none" for the duplication-event
simulate = function(N, G, GC, u, p, A_B, mode, focus)
{
    # Initialise genotype arrays
    xd = array(0, dim = c(N, N, G));
    xr = array(0, dim = c(N, N, G));
    xr[1, 2, 1] = 1;
    
    # Set omega
    omega = (1 - u) * A_B / 2;
    gametes = numeric(0);
    
    # Uncomment for empirical distribution
    # warning("Need to set cost through argument")
    # parents = eq_te(N, 10000, u, 5e-02, 5, 1e-12, TRUE)
    # gametes = rep(0, length(parents))
    # for (i in seq_along(parents)) {
    #     for (j in i:1) {
    #         gametes[j] = gametes[j] + parents[i] * dbinom(j - 1, i - 1, 0.5);
    #     }
    # }
    
    # Run each generation, including active generation
    for (g in 1:(G-1)) {
        if (g != GC) {
            x = mutant_lineage(xd, xr, g, u, 0, omega, mode, "none", gametes)
        } else {
            x = mutant_lineage(xd, xr, g, u, p, omega, mode, focus, gametes)
        }
        cat(".");
        xd = x[[1]];
        xr = x[[2]];
    }
    cat("\n");
    
    # Calculate statistics
    A = matrix(0:(N-1), N, N, byrow = FALSE)
    B = matrix(0:(N-1), N, N, byrow = TRUE)

    N_Bd = sapply(1:G, function(g) sum(xd[,,g] * B * (A+B)) / sum(xd[,,g] * B))
    A_Bd = sapply(1:G, function(g) sum(xd[,,g] * B *  A)    / sum(xd[,,g] * B))
    B_Bd = sapply(1:G, function(g) sum(xd[,,g] * B *  B)    / sum(xd[,,g] * B))
    N_Br = sapply(1:G, function(g) sum(xr[,,g] * B * (A+B)) / sum(xr[,,g] * B))
    A_Br = sapply(1:G, function(g) sum(xr[,,g] * B *  A)    / sum(xr[,,g] * B))
    B_Br = sapply(1:G, function(g) sum(xr[,,g] * B *  B)    / sum(xr[,,g] * B))
    
    # Package as data table
    # Slight oddities in generation/stage numbering here are to make this match
    # the deterministic simulations.
    ds = data.table(
        g = c(0:(G-1), 1:G),
        stage = c(rep(2, G), rep(1, G)),
        N_B = c(N_Bd, N_Br),
        A_B = c(A_Bd, A_Br),
        B_B = c(B_Bd, B_Br),
        u = u,
        p = p,
        omega = omega,
        mode = mode
    );
    ds = ds[order(g, stage)];
    
    return (ds)
}

# Helper functions to calculate dW/dn for fitness functions of the form
# w = (1 - UnitCost) ^ (0.5*n^2), where in this context capital W means
# relative fitness of the host (i.e. host fitness divided by average host
# fitness) and n is mean TE copy number over hosts.

# "Naive" dW/dn: point approximation (i.e. technically this is the slope
# of W evaluated at point n; other versions don't really calculate dW/dn but
# a least-squares slope of W across a distribution of host genotypes).
dWdn = function(n, UnitCost)
{
    log(1 - UnitCost) * n
}

# dWdn_poisson: Assumes that TEs are distributed across hosts according to a
# Poisson or a negative binomial distribution with mean n and variance factor vf 
# i.e. Var(X)/E(X) (vf = 1 is Poisson, vf > 1 is negative binomial) and 
# calculates the least-squares slope of W on n across hosts.
dWdn_poisson = function(n, UnitCost, vf = 1)
{
    mapply(dWdn_poisson_work, n, UnitCost, vf)
}

# Workhorse for dWdn_poisson
dWdn_poisson_work = function(n, UnitCost, vf)
{
    nmax = floor(3 * n + 20)
    nx = 0:nmax
    if (vf == 1) {
        nw = dpois(nx, n)
    } else {
        size = n / (vf - 1)
        nw = dnbinom(nx, size = size, mu = n)
    }
    ny = (1 - UnitCost) ^ (0.5 * nx^2)
    
    ww = weighted.mean(ny, nw)
    
    cm = cov.wt(cbind(nx, ny), nw, method = "ML")$cov
    (cm[2,1]/cm[1,1])/ww
}

# dWdn_mutant: Like dWdn_poisson, but weights according to mutant TE's
# perspective. Note n here is still equilibrium n across hosts.
dWdN_mutant = function(a, u, cost, vf = 1)
{
    mapply(dWdN_mutant_work, a, u, cost, vf)
}

# Workhorse for dWdN_mutant
dWdN_mutant_work = function(a, u, cost, vf)
{
    # Size of genotype matrix
    nmax = floor(6 * a + 20)
    mmax = floor(10 * 1/(1 - u))

    # Genotype matrix
    if (vf == 1) {
        X = dpois(0:(nmax - 1), a) %o% 
                dquasi1(1:mmax, u);
    } else {
        X = dnbinom(0:(nmax - 1), size = a / (vf - 1), mu = a) %o% 
                dquasi1(1:mmax, u);
    }
    
    A = matrix(0:(nmax-1), nrow = nmax, ncol = mmax)
    B = matrix(1:mmax, nrow = nmax, ncol = mmax, byrow = TRUE)

    abs_w = weighted.mean((1 - cost) ^ (0.5 * A^2), X);
    w = (1 - cost) ^ (0.5 * (A+B)^2) / abs_w;

    cm = cov.wt(cbind(c(A+B), c(w)), c(X*(A+B)), method = "ML")$cov
    (cm[2,1]/cm[1,1])
}

# This takes a different approach -- rather than assuming a given parametric
# distribution for TE copy number across hosts, simulates the duplication - 
# selection - recombination process for TEs with duplication rate u and the
# fitness function parameterised by UnitCost to get the distribution
# of TE number across hosts and then measures dW/dn across those hosts.
# Note that the n parameter here is just to estimate the number of genotypes
# needed to capture dynamics and does not affect the equilibrium distribution
# sought after.
dWdn_sim = function(n, UnitCost, u)
{
    # Create recycled args
    a = t(mapply(rbind, n, UnitCost, u));
    a = as.data.table(a);
    names(a) = c("nn", "UnitCost", "u");
    a[, id := 1:.N]
    
    # Apply to each unique UnitCost & u
    a[, dWdn := {
            nmax = floor(3 * max(n) + 20)
            nw = eq_te(nmax + 1, 10000, u[1], UnitCost[1], 10, 1e-10)
            dWdn_sim_work(UnitCost[1], nw = nw)
        }, by = .(UnitCost, u)]
    
    # Put back in correct order
    a = a[order(id)]
    a[, dWdn]
}

# Workhorse for dWdn_sim
dWdn_sim_work = function(UnitCost, nw)
{
    nmax = length(nw) - 1
    nx = 0:nmax
    ny = (1 - UnitCost) ^ (0.5 * nx^2)
    
    ww = weighted.mean(ny, nw)
    
    cm = cov.wt(cbind(nx, ny), nw, method = "ML")$cov
    (cm[2,1]/cm[1,1])/ww
}

# Helper functions to assemble a string giving expressions and their values
values = function(...)
{
    names = as.character(substitute(list(...))[-1])
    vals = list(...)
    paste(paste(names, "=", vals), collapse = ", ")
}

# Quasigeometric distribution density function, with support on 1..Inf
dquasi1 = function(x, u)
{
    dgeom(x - 1, 1 - u)
}

# Negative binomial distribution density function, parameterised by mean and variance factor
dnbinom1 = function(x, mu, vf)
{
    dnbinom(x, size = mu / (vf - 1), mu = mu)
}

choose2 = function(n, k)
{
    gamma(n + 1) / (gamma(k + 1) * gamma(n - k + 1))
}

dnbinom2 = function(x, mu, vf)
{
    choose2(x+mu/(vf-1)-1, mu/(vf-1)-1) * (1/vf)^(mu/(vf-1)) * (1 - 1/vf)^x
}

# Helper functions to calculate equilibrium mean copy number for fitness 
# functions of the form w = (1 - cost) ^ (0.5*n^2) at duplication rate u.

# Charlesworth point-estimate method.
eqn_cw = function(u, cost)
{
    # u == -log(1 - cost) * n
    -u / log(1 - cost)
}

eqn_vf = function(u, cost)
{
    # u == -log(1 - cost) * n * vf
    vf = 1 + 2.6 * u + 1.4 * cost - 1.9 * u * cost;
    -u / (log(1 - cost) * vf)
}

eqn_te = function(u, cost, precision = 2, epsilon = 1e-6)
{
    mapply(eqn_te_work, u, cost, precision, epsilon)
}

eqn_te_work = function(u, cost, precision = 2, epsilon = 1e-6)
{
    cat(".")
    n0 = eqn_vf(u, cost);
    nmax = max(20, n0 * precision);
    x = eq_te(nmax, 10000, u, cost, n0, epsilon, TRUE);
    weighted.mean(0:(length(x) - 1), x)
}

# Uses external ./determ executable to run simulation to equilibrium and returns 
# results. Note, the results are those after proliferation.
ext_determ = function(u, cost, gmax = 10000, nmax = 2000, epsilon = 1e-6, verbose = FALSE)
{
    n0 = eqn_vf(u, cost);
    vf = 1 + 2.6 * u + 1.4 * cost - 1.9 * u * cost;
    n_max = min(nmax, max(20, qnbinom(1e-8, size = n0 / (vf - 1), mu = n0, lower.tail = FALSE)));

    system2("./Deterministic/determ", args = c("-n_max", n_max, "-generations", as.integer(gmax),
        "-u", u, "-v", 0, "-w", paste0("\"n->(1-", cost, ")^(0.5*n^2)\""), "-epsilon", epsilon,
        "-fileout", "./Deterministic/output.ngd"),
        stdout = ifelse(verbose, "", FALSE))
    
    output = fread("./Deterministic/output.ngd", skip = "g\t")
    output = as.list(output[.N])
    output$n_max = n_max
    
    output
}

# Uses external ./determ executable to run invasion simulation and returns results.
ext_invade = function(u, cost, type, gens = 100, nmax = 100, nbmax = 15, epsilon = 1e-12, verbose = FALSE)
{
    n0 = eqn_vf(u, cost);
    vf = 1 + 2.6 * u + 1.4 * cost - 1.9 * u * cost;
    n_max = min(nmax, max(20, qnbinom(1e-8, size = n0 / (vf - 1), mu = n0, lower.tail = FALSE)));
    
    if (n_max == nmax) {
        warning("n_max == nmax for u ", u, " cost ", cost)
    }

    system2("./Deterministic/determ", args = c("-generations", 0L, "-epsilon", epsilon,
        "-gamete_threshold", 0L, "-transition_threshold", 0L, "-marginal_threshold", 0L, 
        "-analysis_gens0", -1L, "-analysis_gens", gens, "-analysis_seed", 1e-6,
        "-analysis_nb_max", nbmax, "-n_max", n_max, "-mode", "analysis2", 
        "-type", type, "-n_0", 0.1, "-w", paste0("\"n->(1-", cost, ")^(0.5*n^2)\""),
        "-u", u, "-u_B", u + 0.0001, "-v", 0.0, "-invasion_seed", 0.001, "-report_freq", 1000L,
        "-fileout", "./Deterministic/output.ngd"),
        stdout = ifelse(verbose, "", FALSE))
    
    output = fread("./Deterministic/output.ngd", skip = "g\t")
    
    b0 = output[g == 1.1, b]
    b1 = output[g == 100.1, b]
    
    (b1 - b0) / b0
}


# If simulate == FALSE, get a rough equlibrium n and equilibrium dW/dn
# through analytic approximations. If simulate == TRUE, actually simulate
# a 1-type transposon model to get equilibrium n and dW/dn. Note that 
# simulate == TRUE is still part of the "analytic" approach, as it then 
# calculates R, B and C in the normal way.
# Then, with:
# how == "product" use product definitions of B and C
# how == "sum"     use sum approximations to product definitions of B and C
# how == "simple"  use simple approximations to B and C (B = 1, C = 2u basically)
get_RBC = function(u, cost, simulate = FALSE, how = "product", mode = "private")
{
    # Estimate equilibrium n
    n_est = eqn_vf(u, cost);
    
    # Definition of fitness function
    wf = function(n) (1 - cost) ^ (0.5*n^2);
    
    if (simulate) {
        # Try out a value for n_max
        n_max = qpois(-1e-10, n_est, lower.tail = TRUE, log.p = TRUE);
        n_max = max(n_max, 25)
        if (n_max > 500) {
            stop(n_max)
        }
        # Get equilibrium TE distribution
        X = eq_te(n_max = n_max, g_max = 5000, u = u, cost = cost, 
            n_0 = 10, epsilon = 1e-16, after_selection = FALSE);

        # Measure average n, average w, and slope of host relative fitness on transposon number
        nn = 0:(n_max - 1);
        n_avg = weighted.mean(nn, X); # Note: this is average n after proliferation.
        w_avg = weighted.mean(wf(nn), X);
        cov_wn = cov.wt(cbind(n = nn, w = wf(nn) / w_avg), wt = X, method = "ML")$cov;
        beta = cov_wn["n", "w"] / cov_wn["n", "n"];
    } else {
        # Just do everything through analytics
        n_avg = n_est;
        w_avg = wf(n_avg);
        beta = log(1 - cost) * n_avg;
    }
    n = n_avg; # This name is more convenient below...

    # Now, calculate various key quantities...
    W0 = 1 / (1 + u);
    
    if (mode == "private") {
        N1Off = 1/(1+u)
        R = 1
    } else {
        # get around issues with exp(n) overflowing or with the division underflowing
        if (n < 100) { 
            dN1 = (exp(n) * (n * (1 + 4*u) * (1 + n + n*u) - 4*u) - n + 4*u) / 
                (exp(n) * (n + 4 * (n - 1) * u) - n + 4*u);
        } else {
            dN1 = (n * (1 + 4*u) * (1 + n + n*u) - 4*u) / 
                (n + 4 * (n - 1) * u);
        }

        N1Off = dN1 - 
            n * (1 + u) - 
            (1 + u * (1 + u) * (7 + 4*u)) / 
            ((1 + u) * (1 + 2 * u));
        R = (1 - exp(-n) + 4*u*(n - 1 + exp(-n))/n)/(n + 2*u*n)
    }
    N1Par = (1 + u) / (1 - u)
    
    if (how == "product") {
        t = 1:100
        B = prod(1 + beta/W0 * N1Off * (exp(beta) * (1+u)/2)^(t - 1))
        C = 1 - prod(1 + beta/W0 * N1Par * (exp(beta) * (1+u)/2)^(t - 1))
    } else if (how == "sum") {
        B = 1 + beta/W0 * N1Off * 2/(2 - exp(beta) * (1 + u))
        C = -beta/W0 * N1Par * 2/(2 - exp(beta) * (1 + u))
    } else if (how == "simple") {
        if (mode == "private") {
            R = 1
            B = 1
            C = 2*u
        } else {
            R = (1 - exp(-n_avg)) / n_avg
            # Can use the below, but the improvement is small.
            # R = (1 - exp(-n) + 4*u*(n - 1 + exp(-n))/n)/(n + 2*u*n)
            B = 1
            C = 2 * u
        }
    }

    list(R = R, B = B, C = C)
}
