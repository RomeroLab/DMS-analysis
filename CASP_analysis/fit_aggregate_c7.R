remove(list=ls())
library(pudms)
source("fn_fit_aggregate.R")
#### example rocker dataset

dat1 = pudms::create_protein_dat(path_l = "Casp7_Lib1_sorted_sequences.txt",
                                path_u = "Casp7_Lib1_unsorted_sequences.txt")
dat2 = pudms::create_protein_dat(path_l = "Casp7_Lib2_sorted_sequences.txt",
                                path_u = "Casp7_Lib2_unsorted_sequences.txt")
dat3 = pudms::create_protein_dat(path_l = "Casp7_Lib3_sorted_sequences.txt",
                                path_u = "Casp7_Lib3_unsorted_sequences.txt")

### for illustration, divide example dataset into three subsamples
### pretend they are from three replicated experiments
#shuffle_idx=cut(sample(1:nrow(dat)),breaks = 3,labels = F)
rep1=dat1
rep2=dat2
rep3=dat3

rep_all = list(rep1,rep2,rep3) # put them in a list

### obtain coefficients and inverse fisher information matrices
py1_rocker  = 0.01
fit_all<- lapply(1:3, fit_repi ,rep_all = rep_all, py1_rep = py1_rocker)

# collect coefficients and inv matrices
coef_all <- lapply(1:3, function(i) fit_all[[i]]$coef)
invI_all <- lapply(1:3, function(i) fit_all[[i]]$invI)

# aggreate results from reps
pvalues = aggregate_reps(coef_all = coef_all,invI_all = invI_all)
blockidx=gsub(sapply(1:length(pvalues$coef[-1]), function(i) strsplit(names(pvalues$coef[-1]),"\\.")[[i]][1]),pattern = "[^0-9]",replacement = "")

# summarise results from reps
res_table<-return_table_for_aggregated_dat(pvalues = pvalues, blockidx = blockidx)

# export result
write.csv(res_table[-1,], file="C7_aggregated_result.csv")
