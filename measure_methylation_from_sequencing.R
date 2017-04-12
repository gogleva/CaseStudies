# measure methylation from sequencing

#Reduced Representation Bisulfite Sequencing or RRBS is an experimental technique widely used to manipulate the regions of the genome we measure. An enzyme is used to cut DNA at CCGG and the general idea is to filter out small or large molecules once DNA is cut. We can use Bioconductor tools to predict the size of these regions. Load the genome package and create an object with the sequence for chr22:

library("BSgenome.Hsapiens.UCSC.hg19")
chr22 = Hsapiens[["chr22"]]

#use the matchPattern function to find all the locations in which CCGG occurs on chr22.

#-----Q1: How many CCGG do we find on chr22?
res = matchPattern('CCGG', chr22)
length(res)

# A = 58102

#----Q2: Plot a histogram of the DNA fragment sizes after we cut with CCGG. How would you describe this distribution?

size=diff(start(res))
hist(size)
hist(log10(size))

#A : The distribution has a long right tail with most values from 0-1000, but some very large values.

#-----Q3: What proportion of the fragments created for chr22 are between 40 and 220 basepairs?

up <- size[which(size <= 220)]
both <- which(up >= 40)

length(both)/length(size)

#=> A = 0.3309754

#-----Q4: If we sequence all of chromosome 22 we need to sequence 51,304,566 bases. However, if instead we focus only on fragments of size between 40 and 220 basepairs, how much sequence would we need?

target_sizes <- up[both]
res =  sum(target_sizes)

# a bit more elegant:

sum( size[size<=220 & size>=40] )

# => A = 2203342

#-----Q5: 

path = "/home/anna/anna/study/DNA_methylation/colonCancerWGBS-master"
targets = read.table(file.path(path,"targets.txt"), header = TRUE, sep = "\t")
targets

biocLite("bsseq")
library("bsseq")
cov.files = list.files(path=path,pattern="*chr22.cov",full.names=TRUE)  # coverage files
colonCancerWGBS =read.bismark(files=cov.files, rmZeroCov=TRUE, sampleNames = as.character(targets$Run), strandCollapse = FALSE)
# add sample information to object
colData(colonCancerWGBS) = DataFrame(targets)
###Note you might see a warning message here. You can ignore.

#To view the bsseq object and the phenotypic information about each sample:

colonCancerWGBS
# phenotypic information
pData(colonCancerWGBS)
# granges object
granges(colonCancerWGBS)


#Now we can extract the coverage and the number of reads with evidence from methylation:

cov=getCoverage(colonCancerWGBS,type = "Cov")
m=getCoverage(colonCancerWGBS,type = "M")

filter <- apply(cov, 1, function(x) length(x[x != 0])>= 6)
some_coverage <- cov[filter,]

#Q: What proportion of the reported CpGs have some coverage in all sample?

dim(some_coverage)[1]/nrow(m)

# A = 0.7743644

#-----Q6: Compute the total coverage (across all samples) for each CpG. Plot it against location.

#total coverage:

tot_cov <- rowSums(cov)
pos <- start(granges(colonCancerWGBS))

plot(pos, tot_cov)

# => A: Has some very large values (>200) as well as general varaibility


#-----Q7
#we can get coverage and the number of reads including evidence for methylation like this:
cov=getCoverage(colonCancerWGBS,type = "Cov")
m=getCoverage(colonCancerWGBS,type = "M")

#make plot of a selected region:
gr = GRanges(seqnames="22",ranges=IRanges(start=43793678,end= 45022550))
index=granges(colonCancerWGBS)%over%gr
library(rafalib)
i=1
index2=which(index & cov[,i]>=5 & cov[,i]<=50)
x=start(colonCancerWGBS)[index2]
y=m[index2,i]/cov[index2,i]
w=sqrt(cov[index2,i])/7
plot(x,y,cex=w)



