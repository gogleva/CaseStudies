# Cell composition

#To examine the importance of accounting for cellular composition in DNA methylation analysis, we are going to download a GEO dataset used in an analysis of whole blood data.

#The minfi package has a function to read data directly from GEO. Run the following commands. Note that this command 
#downloads 64.7 MB of data and can take several minutes, depending on your download speed.

library(minfi)
grset = getGenomicRatioSetFromGEO("GSE32148")
save(grset,file="grset.rda")
load("grset.rda")


#This creates an object of class:
class(grset)

# use pData to examine the sample information table
pData(grset)

#-----Q1: which column includes the age of the individual?
View(pData(grset))

# A = characteristics_ch1.1


#-----Q2
#extract the age as a number of the individuals sampled here:
age=pData(grset)$characteristics_ch1.1
age=as.character(age)
age[grep("N/A",age)] = NA
age=as.numeric(gsub("age \\(y\\):\ ","",age))

# This experiment was performed to find DMRs between 
# individuals with Crohn's disease and controls. We can
# extract this information like this:


group = rep("normal",nrow(pData(grset)))
group[grepl("ulcerative",pData(grset)[,1])]="ulcerative"
group[grepl("Crohn",pData(grset)[,1])]="crohn"
group = factor(group,levels=c("normal","ulcerative","crohn"))

#create some exploratory plots based on distance between sample. Before doing this, we will remove CpGs with NA calls as well as the sex chromosomes:

keep = which(rowSums(is.na(getBeta(grset)))==0 & 
               !seqnames(grset)%in%c("chrX","chrY"))
##create a new object 
grset2=grset[keep,]

y = getBeta(grset2)
mds = cmdscale( dist(t(y)))
#age_cols <- ifelse(pData(grset2)$Tissue == 'breast', 'red', 'green')

age_cols = ifelse(age >= 40, 'red', 'green')
cod_pch = ifelse(group == 'normal', 6, ifelse(group == 'crohn', 8, 19))
plot(mds, col = age_cols, pch = cod_pch)

legend('topright', col= c('red', 'green'), legend= c("old", "young"), pch = 16, cex = 0.7)
legend('bottomleft', pch = c(6,8,19), legend = c('normal', 'crohn', 'ulcer'))

# A: The individuals that are older than 40 form a cluster, with the rest showing larger variability.


