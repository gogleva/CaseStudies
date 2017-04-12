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

#A: 0.3309754

#-----Q4: If we sequence all of chromosome 22 we need to sequence 51,304,566 bases. However, if instead we focus only on fragments of size between 40 and 220 basepairs, how much sequence would we need?

target_sizes <- up[both]
res =  sum(target_sizes)

# => A = 2203342




