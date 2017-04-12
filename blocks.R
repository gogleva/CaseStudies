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

#-----Q1: How many regions are represented in the collapsed object?
nrow(cdat$object)

# A = 223497

#-----Q2: What proportion of the regions are OpenSea regions?

OpenSeq_ratio <- length(which(granges(cdat$obj)$type == 'OpenSea'))/length(granges(cdat$obj)$type)


#-----Q3: blockfinder to find differentially methylated regions between cancer and normal:

status = factor(pData(cdat$obj)$Status,
                level=c("normal","cancer"))
X=model.matrix(~status)
res = blockFinder(cdat$obj,X,cutoff=0.05)

#blockFinder calls bumphunter and returns a similar object. We can see the blocks:

head(res$table)

#Q: What proportion of the blocks reported in res$table are hypomethyated (lower methylation in cancer versus normal)?
hypo <- length(which(res$table$value < 0))/nrow(res$table)

