# Cell composition

#To examine the importance of accounting for cellular composition in DNA methylation analysis, we are going to download a GEO dataset used in an analysis of whole blood data.

#The minfi package has a function to read data directly from GEO. Run the following commands. Note that this command 
#downloads 64.7 MB of data and can take several minutes, depending on your download speed.

library(minfi)
grset=getGenomicRatioSetFromGEO("GSE32148")

#This creates an object of class:
class(grset)

# use pData to examine the sample information table
pData(grset)

#Q1: which column includes the age of the individual?