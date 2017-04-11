#Data analysis and integration, practice

library(minfi)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

grset=getGenomicRatioSetFromGEO("GSE32148")
##after donwload save the object like this:
save(grset,file="grset.rda")
##then to load later
load("grset.rda")

#Q1:
#specify path to TCGA data:
path="/home/anna/anna/study/DNA_methylation/tcgaMethylationSubset-master"
targets=read.delim(file.path (path,"targets.txt"),as.is=TRUE)
#How many samples are represented in this table?
dim(targets)
#[1] 98 65
# A1 = 98


