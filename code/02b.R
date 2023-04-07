library(data.table); library(ggplot2); library(patchwork)
load("./DATA/fst_estimates_subsampled.RData")
n_nearby <- 200; min_maf <- 0.05; rsid_status <- 'keep_good_only'
info <- info[!is.na(info$maf.ML_same),]

######## Split off neutral sites, then filter candidates based on minor allele frequency ########
info_neut <- info[info$type == "neut",]; info_neut$maf_rank <- rank(info_neut$maf.ML_same)
info <- info[info$maf.ML_same >= min_maf,]
########

######## calculate p-values ########
library(parallel); pvals <- info[,1:5]

calc_pval_L13 <- function(site) {
  input_maf <- info$maf.ML_same[site]; input_fst <- info$london.post.pre.fst.ML_same[site]
  n_less <- sum(info_neut$maf.ML_same <= input_maf)
  
  tmp_fst <- info_neut$london.post.pre.fst.ML_same[info_neut$maf_rank <= (n_less + n_nearby/2) & info_neut$maf_rank >= n_less - n_nearby/2]
  if (length(tmp_fst) < n_nearby) {tmp_fst <- info_neut$london.post.pre.fst.ML_same[info_neut$maf_rank >= length(info_neut$maf_rank) - n_nearby]}
  
  tmp_percentiles <- ecdf(tmp_fst)
  tmp_p <- 1-tmp_percentiles(input_fst)
  return(tmp_p)
}
pvals$L13.pval <- do.call("c",mclapply(1:nrow(info), FUN = calc_pval_L13))
calc_pval_L12 <- function(site) {
  input_maf <- info$maf.ML_same[site]; input_fst <- info$london.during.pre.fst.ML_same[site]
  n_less <- sum(info_neut$maf.ML_same <= input_maf)
  
  tmp_fst <- info_neut$london.during.pre.fst.ML_same[info_neut$maf_rank <= (n_less + n_nearby/2) & info_neut$maf_rank >= n_less - n_nearby/2]
  if (length(tmp_fst) < n_nearby) {tmp_fst <- info_neut$london.during.pre.fst.ML_same[info_neut$maf_rank >= length(info_neut$maf_rank) - n_nearby]}
  
  tmp_percentiles <- ecdf(tmp_fst)
  tmp_p <- 1-tmp_percentiles(input_fst)
  return(tmp_p)
}
pvals$L12.pval <- do.call("c",mclapply(1:nrow(info), FUN = calc_pval_L12))

calc_pval_D13 <- function(site) {
  input_maf <- info$maf.ML_same[site]; input_fst <- info$denmark.post.pre.fst.ML_same[site]
  n_less <- sum(info_neut$maf.ML_same <= input_maf)
  
  tmp_fst <- info_neut$denmark.post.pre.fst.ML_same[info_neut$maf_rank <= (n_less + n_nearby/2) & info_neut$maf_rank >= n_less - n_nearby/2]
  if (length(tmp_fst) < n_nearby) {tmp_fst <- info_neut$denmark.post.pre.fst.ML_same[info_neut$maf_rank >= length(info_neut$maf_rank) - n_nearby]}
  
  tmp_percentiles <- ecdf(tmp_fst)
  tmp_p <- 1-tmp_percentiles(input_fst)
  return(tmp_p)
}
pvals$D13.pval <- do.call("c",mclapply(1:nrow(info), FUN = calc_pval_D13))

rm(list=ls(pattern="calc_pval_"))
rm(drop_samples_missing, design, design2, info_neut)
########


pvals <- pvals[,-c(2:5)]; colnames(pvals)[2:4] <- paste(colnames(pvals)[2:4],"ML_same", sep=".")
info <- merge(pvals, info, by="site"); rm(pvals)
########

save.image("./DATA/pvalues_200neutralsites_subsampled.RData")

