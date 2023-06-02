# repMeta
R package for quantifying replicability of multiple studies in a meta-analysis

# Download
> install.packages("devtools") # comment out if you have the devtools package
>
> devtools::install_github("menglix/repMeta")

# Required R packages

> metafor

> evd

> parallel

# Minimum toy example to reproduce the case study in the paper

## load required packages

> library(metafor)
> 
> library(evd)
>
> library(parallel)

## transform the data

### the data should contain the columns with effect size named by "y" and within-study variance named by "s2" 

> data.case <- data.trans.bin(data=moller12,ai = r1, n1i = n1, ci = r2, n2i = n2,measure="OR")

## calculate $R_1$ for the given data
> calR(data.case,m=1)

## replicability test

### replicability test based on Gumbel approximation
> pRmasym(data.case,m=1)

### replicability test based on the bootstrap approximation
> pRmboot(data.case,m=1)

## identify studies with non-replicable results

### identify studies using Gumbel approximation
> Rm.func.iterative(data.case,m=1)

### identify studies using bootstrap approximation
> Rm.func.iterative.boot(data.case,m=1)

