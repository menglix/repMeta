# repMeta
R package for quantifying replicability of multiple studies in a meta-analysis

# Requirements

- R (> 3.5)
- `metafor`
- `evd`
- `parallel`

# Download
```R
install.packages("devtools") # comment out if you have the devtools package
devtools::install_github("menglix/repMeta")
```

# Minimum toy example to reproduce the case study in the paper

```R
## Load required packages
library(metafor) 
library(evd)
library(parallel)

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
