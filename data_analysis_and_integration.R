#Data analysis and integration, practice

library(minfi)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

grset=getGenomicRatioSetFromGEO("GSE32148")
##after donwload save the object like this:
save(grset,file="grset.rda")
##then to load later
load("grset.rda")

#-----Q1:
#specify path to TCGA data:
path="/home/anna/anna/study/DNA_methylation/tcgaMethylationSubset-master"
targets=read.delim(file.path (path,"targets.txt"),as.is=TRUE)
#How many samples are represented in this table?
dim(targets)
#[1] 98 65
# A1 => 98


#-----Q2:
#How many samples are from normal colon samples?
length(which(targets$Tissue == 'colon' & targets$Status == 'normal'))
# A2 => 17


#----Q3:
#read in breast and colon normal samples:
index = which( targets$Status=="normal" & targets$Tissue%in%c("colon","breast") )
targets = targets[index,]
#Now we are ready to read in the data (this will take about 2 minutes):
dat = read.metharray.exp(base=path,targets = targets, verbose=TRUE)

#dat includes the raw data. To convert this into an object that includes methylation values, 
#as well as the location of CpGs, we do the following (we show you the class of dat as we transform it):

class(dat)
## preprocess the data
dat = preprocessIllumina(dat)
class(dat)
## assign locations to each CpG
dat = mapToGenome(dat)
class(dat)
## precompute methylation values from U and M values
dat = ratioConvert(dat,type="Illumina")
class(dat)

#create some quality assessment plots. First look at the distribution of each sample:

library(rafalib)
mypar(1,1)
##extract methylation values
y = getBeta(dat)
shist(y)

# => distributions seem similar. Nothing stands out
mds = cmdscale( dist(t(y)))
tissue = as.factor(pData(dat)$Tissue)
plot(mds,col=tissue)
#create an MDS plot to search for outlier samples. 

# => The first PC splits the data by tissue as expected and no sample stands out as an outlier.


# we are ready to use statistical inference to find differentially methylated regions.
# start by using the limma package to perform a site-by-site analysis.

library(limma)
##create design matrix
tissue = as.factor(pData(dat)$Tissue)
X = model.matrix(~tissue)
##extract methylation values
y = getBeta(dat)
## obtain effect sizes and pvals with limma
fit = lmFit(y,X)

#Q3.1: Find the CpG with the largest effect size when comparing the two tissues. What chromosome is it on?
which.max(fit$coefficients[,2])
#to get the CpG name:
names(which.max(fit$coefficients[,2]))
#get the chr and position for the CpG
granges(dat["cg22365276",])
#1-liner:
granges(dat[names(which.max(fit$coefficients[,2])),])







