# repMeta
R package for quantifying replicability of multiple studies in a meta-analysis

# Download
install.packages("devtools") # comment out if you have the devtools package

devtools::install_github("menglix/repMeta")

# Minimum toy example to reproduce the case study in the paper
library(metafor)

library(evd)

library(parallel)

data.case <- data.trans.bin(data=moller12,ai = r1, n1i = n1, ci = r2, n2i = n2,measure="OR")

calR(data.case,m=1)

pRmasym(data.case,m=1)

pRmboot(data.case,m=1)

Rm.func.iterative(data.case,m=1)

Rm.func.iterative.boot(data.case,m=1)

