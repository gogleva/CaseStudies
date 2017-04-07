#Data analysis and integration, practice

library(minfi)
grset=getGenomicRatioSetFromGEO("GSE32148")

##after donwload save the object like this:
save(grset,file="grset.rda")
##then to load later
load("grset.rda")



