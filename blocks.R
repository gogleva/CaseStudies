# Blocks

library(minfi)
path="/home/anna/anna/study/DNA_methylation/tcgaMethylationSubset-master"
targets=read.delim(file.path (path,"targets.txt"),as.is=TRUE)
index = which( targets$Tissue=="colon")
targets = targets[index,]
dat = read.metharray.exp(base=path,targets = targets, verbose=TRUE)

#transform dat to get methylation avlues and locations of CpGs

dat = preprocessIllumina(dat)
dat = mapToGenome(dat)

# collapse the data:

cdat = cpgCollapse(dat)
nrow(dat)

#Q: How many regions are represented in the collapsed object?
nrow(cdat$object)

# A = 223497