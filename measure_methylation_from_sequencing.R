# measure methylation from sequencing

#Reduced Representation Bisulfite Sequencing or RRBS is an experimental technique widely used to manipulate the regions of the genome we measure. An enzyme is used to cut DNA at CCGG and the general idea is to filter out small or large molecules once DNA is cut. We can use Bioconductor tools to predict the size of these regions. Load the genome package and create an object with the sequence for chr22:

library("BSgenome.Hsapiens.UCSC.hg19")
chr22 = Hsapiens[["chr22"]]

#use the matchPattern function to find all the locations in which CCGG occurs on chr22.

#-----Q1: How many CCGG do we find on chr22?
res = matchPattern('CCGG', chr22)
length(res)

# A = 58102

#----Q2: 

size=diff(start(res))
hist(log10(size))
