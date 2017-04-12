#Data analysis and integration, practice

library(minfi)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

grset=getGenomicRatioSetFromGEO("GSE32148")
##after donwload save the object like this:
save(grset,file="grset.rda")
##then to load later
load("grset.rda")

#-----HW08
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


#Q4: use the qvalue function to determine the q-value for the CpG found in the previous question.

library(qvalue)
##create design matrix
tissue = as.factor(pData(dat)$Tissue)
X = model.matrix(~tissue)
##extract methylation values
y = getBeta(dat)
## obtain effect sizes and pvals with limma
fit = lmFit(y,X)
eb = ebayes(fit)
## obtain q-values
qvals = qvalue(eb$p.value[,2])$qvalue

#Q: What is the q-value for this CpG?
qvals[360649]
#A => 1.269936e-27 


#Q5: 
#use the location of the CpG discussed in the previous two questions. 
#Find all the CpGs within 5000 basepairs of the location of this CpG. 

i <- 360649
gr <- granges(dat)
middle <- gr[i, ]
Index <- gr%over%(middle + 5000)
cols <- ifelse(pData(dat)$Tissue == 'breast', 'red', 'green')
chr <- as.factor(seqnames(gr))
pos <- start(gr)

#Create a plot showing the methylation values for all samples for these CpGs.
#Use color to distinguish breast from colon
#methylation values are stored in y

matplot(pos[Index], y[Index,], col=cols, pch=16, xlab="genomic location", ylab="methylation")
legend("bottomright", inset=.05, legend=c("breast", "colon"), pch=1, col=c("red", "green"), horiz=TRUE)

#Plot the estimated effect size and the -log10 of the q-value in two separate plots for a total of three plots.
splot(fit$coef[,2],-log10(eb$p.value[,2]),xlab="Effect size",ylab="-log10 p-value")
plot(pos[Index],fit$coef[Index,2],type="b",xlab="genomic location",ylab="difference")

