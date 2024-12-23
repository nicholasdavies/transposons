# The tragedy of the transposons

Code and data for manuscript by Nicholas G. Davies, J. Arvid Ågren, Kevin R. Foster.

## Guide to the repository

`Cost-Benefit/` contains intermediate files processed by `cost_benefit_load.R` from data files in `Deterministic/Runs/4-Analysis/`.

`Deterministic/` contains the one-type and two-type deterministic transposon simulation code.

`TE/` contains the stochastic individual-based transposon simulation code.

The code in `Deterministic/` and `TE/` is written in C++ and can be compiled with `make`.

Other files in the top level are R scripts and intermediate data files generated by those scripts, which mainly have to do with generating figures for the manuscript. These were run using R 4.4.1.

### System requirements

To run the scripts, R (version 4.0.0 or higher) is required with the packages `cmocean`, `cowplot`, `data.table`, `ggarchery`, `ggpattern`, `ggplot2`, `Rcpp`, `rlang`, `scales`, and `stringr`. A C++ compiler (C++14) is required with the Lua library installed.

The software has been tested with R 4.4.1 and the clang compiler running on macOS 14.2.1.

### Installation guide

To run the R scripts, download the repository (or use `git clone`) and run the R scripts in the top level directory. There is an RStudio project file that can be used to navigate the code. Downloading R and RStudio and installing the packages can be done in around 20 minutes depending on the speed of your internet connection.

The code in the folders `Deterministic/` and `TE/` is written in C++ and can be compiled with `make`. Instructions for installing a C++ compiler such as gcc can be found online.

### Demo

After building executables in the folders `Deterministic/` and `TE/` using `make`, the resulting executables can be run from the command line with no parameters for some demonstration simulations. The parameters used by default for the simulations will be printed to the screen and results will be written to the default output file in each directory. The expected run time of each demo is less than a second.

### Instructions for use

Config files for running the analyses reported in the manuscript can be found in `Deterministic/Runs/` and `TE/Runs/`. For example, from the `TE/` directory, you can run from the command line
```
te ./Runs/1-Tragedy/config.cfg
```
to run the analyses in the `TE/Runs/1-Tragedy` directory.

For the R scripts in the top level directory, these can be run from start to finish to reproduce the figures.
