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

# => A: A region of about 1000 base pairs appears to be different.


#-----Q6. Repeat the above exercise, but now make the same plots for the top 10 CpGs ranked by absolute value of effect size. 
#get the order like this:
o = order(abs(fit$coef[,2]), decreasing = TRUE)[1:10]

par(mfrow = c(2,5))

for (i in o) {
  middle <- gr[i, ]
  Index <- gr%over%(middle + 5000)
  cols <- ifelse(pData(dat)$Tissue == 'breast', 'red', 'green')
  chr <- as.factor(seqnames(gr))
  pos <- start(gr)
  matplot(pos[Index], y[Index,], col=cols, pch=16, xlab="genomic location", ylab="methylation")
#  legend("bottomright", inset=.05, legend=c("breast", "colon"), pch=1, col=c("red", "green"), horiz=TRUE)
}

# => A: For most plots we see groups of CpGs that are differentially methylated.

#-----Q7: 

# explicitly search for regions using the bumphunter function.
#We will use permutation to assess statistical significance.
#Because the function is slow, we will restrict our analysis to chromosome 15.

index= which(seqnames(dat)=="chr15")
dat2 = dat[index,]

library(doParallel)
ncores = detectCores()
registerDoParallel(cores = ncores)

#We can now run the bumphunter function to find differentially methylated regions (DMR). 
#we will use 100 permutations, although we recommend more in practice. 
#Here we will use a cutoff of 0.1. 
#The permutations are random so make sure you set seed to 1 to obtain exact results in the assessment question.

##create design matrix
tissue = as.factor(pData(dat)$Tissue)
X = model.matrix(~tissue)
##extract methylation values
set.seed(1)
res = bumphunter(dat2,X,cutoff=0.1,B=100)
head(res$tab)

#Q:  how many regions achieve an FWER lower than 0.05?

length(which(res$tab[12] < 0.05))
# A => 296



#-----Q8:

#Previously we performed a CpG by CpG analysis and obtained qvalues.
#Create an index for the CpGs that achieve qvalues smaller than 0.05 and a large effect size larger than 0.5 (in absolute value):

##fit and qvals were defined in a previous answer
index = which(qvals < 0.05 & abs(fit$coef[,2]) > 0.5 & seqnames(dat)=="chr15")

#create a table of the DMRs returned by bumphunter that had 3 or more probes and convert the table into GRanges:
tab = res$tab[ res$tab$L >= 3,]
tab = makeGRangesFromDataFrame(tab)

#Q: What proportion of the CpGs indexed by index are inside regions found in tab
Ind <- gr[index, ]
findOverlaps(Ind, tab)
# => 12 hits
A = 12/length(index)

# => A = 0.5714286


#-----Q9:
#download the table of CGI using AnnotationHub:
library(AnnotationHub)
cgi = AnnotationHub()[["AH5086"]]

# create a GRanges object from the list of DMRs we computed in the previous questions:

tab = res$tab[res$tab$fwer <= 0.05,]
tab = makeGRangesFromDataFrame(tab)


# Q: What proportion of the regions represented in tab do not overlap islands, 
# but overall CpG islands shores (within 2000 basepairs)?
# Hint: use the distanceToNearest

overlap_islands <-findOverlaps(tab, cgi)
map = distanceToNearest(tab,cgi)
d = mcols(map)$dist
tab_shore_and_island <- ind_2000 <- which(d <= 2000)

A = (length(tab_shore_and_island) - length(overlap_islands))/length(tab)

# => A = 0.1892744

#-----Q10
#study the relationship between gene expression and DNA methylation by integrating gene expression and DNA methylation high throughput data
#compare colon and lung samples

path="/home/anna/anna/study/DNA_methylation/tcgaMethylationSubset-master"
targets=read.delim(file.path (path,"targets.txt"),as.is=TRUE)
index = which( targets$Status=="normal" & targets$Tissue%in%c("colon","lung") )
targets = targets[index,]

#read in the data (this will take about 2 minutes):

library(minfi)
dat = read.metharray.exp(base=path,targets = targets, verbose=TRUE)

## preprocess the data
dat = preprocessIllumina(dat)
dat = mapToGenome(dat)
dat = ratioConvert(dat,type="Illumina")

# run the bumphunter function with cutoff=0.25 and default parameters.

tissue=pData(dat)$Tissue
X = model.matrix(~tissue)
res = bumphunter(dat,X,cutoff=0.25)
nrow(res$tab)

# Q: What proportion of regions are just one CpG?
length(which(res$tab$L == 1)) / nrow(res$tab)

# => A = 0.833483


#Q11: 

#match these regions to genes.
#First load the expression related data. 
#this will load an object called tcgaLungColonExpLM which is the result of a differential expression analysis using limma on TCGA rawdata:

path="/home/anna/anna/study/DNA_methylation/tcgaMethylationSubset-master"
load(file.path(path,"tcgaLungColonExpLM.rda"))
class(tcgaLungColonExpLM)

# we saved the annotation of the gene expression array in this object:
print(annotation)

# obtain q-values using the qvalue package:

library(limma)
library(qvalue)
eb=ebayes(tcgaLungColonExpLM)
qvals=qvalue(eb$p.value[,2])$qvalue

# obtain locations for these genes. 

source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("hgu133plus2.db")
library("hgu133plus2.db")

map=select(hgu133plus2.db,keys=rownames(tcgaLungColonExpLM$coef),columns=c("ENTREZID"),keytype="PROBEID")
biocLite("Homo.sapiens")
library(Homo.sapiens)

Genes=genes(Homo.sapiens)
Genes=resize(Genes,1) ## we want the tss

index1=match(as.character(mcols(Genes)$GENEID),map$ENTREZID)
index2 = match(map$PROBEID[index1],rownames(tcgaLungColonExpLM$coef))
M = tcgaLungColonExpLM$coef[index2,2]

#M is now in the same order as Genes. We can now find the closest gene to each DMR.

tab = makeGRangesFromDataFrame(res$tab,keep.extra.columns = TRUE)
map2=distanceToNearest(tab,Genes)

# make plots comparing the methylation differences to the gene expression differences.
# We consider DMRs of different size

index1=subjectHits(map2)
dist = mcols(map2)$dist

library(rafalib)
mypar(2,2)
for(i in c(0,1,2,3)){
  keep = dist< 10000 & tab$L>i
  plot(tab$value[keep],M[index1][keep],main=paste("cor=",signif(cor(tab$value[keep],M[index1][keep],use="complete"),2)))
}

# A: There is a negative correlation between gene expression and DNA methylation and it is stronger for larger DMRs


