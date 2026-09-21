# repMeta
R package for quantifying replicability of multiple studies in a meta-analysis

# Requirements

- R (> 3.5)
- `metafor (>= 3.8.1)`
- `evd (>= 2.3.3)`
- `parallel (>= 3.5.3)`

# Download
```R
install.packages("devtools") # comment out if you have the devtools package
devtools::install_github("menglix/repMeta")
```

> **Note (written by Claude):** repMeta 0.1.1+ works with current `metafor` releases (tested against 3.8.1 and 5.2.1). If you are stuck on repMeta 0.1.0 and cannot upgrade, that version only works with the older `metafor 2.4.0`, bundled in this repo as `metafor_2.4-0.tar.gz` for reference; otherwise ignore that file and just install the latest `metafor` from CRAN.

# Minimum toy example to reproduce the case study in the paper

```R
## Load required packages
library(metafor) 
library(evd)
library(parallel)
library(repMeta)

## Transform the data
### The data should contain the columns with effect size named by "y" and within-study variance named by "s2" 
data.case <- to.dat.repMeta(data=moller12,ai = r1, n1i = n1, ci = r2, n2i = n2,measure="OR")

## Calculate $R_1$ for the given data
calR(data.case,m=1)

## Replicability test
### Replicability test based on Gumbel approximation
pRmasym(data.case,m=1)

### Replicability test based on the bootstrap approximation
pRmboot(data.case,m=1)

## Identify studies with non-replicable results
### Identify studies using Gumbel approximation
Rm.func.iterative(data.case,m=1)

### identify studies using bootstrap approximation
Rm.func.iterative.boot(data.case,m=1)

```
