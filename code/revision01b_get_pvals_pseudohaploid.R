library(data.table); library(ggplot2); library(patchwork)
load("../DATA/fst_estimates_pseudohaploid.RData")
n_nearby <- 250; min_maf <- 0.05; rsid_status <- 'keep_good_only'

# report correlation between ML and pseudohaploid allele frequencies
sqrt(summary(lm(info$maf.ML ~ info$maf))$r.squared)

######## Split off neutral sites, then filter candidates based on minor allele frequency ########
info_neut <- info[info$type == "neut",]; info_neut$maf_rank <- rank(info_neut$maf)
info <- info[info$maf >= min_maf,]
########

######## calculate p-values ########
library(parallel); pvals <- info[,1:5]

calc_pval_L13 <- function(site) {
  input_maf <- info$maf.ML[site]; input_fst <- info$L13.fst[site]
  n_less <- sum(info_neut$maf.ML <= input_maf)
  
  tmp_fst <- info_neut$L13.fst[info_neut$maf_rank <= (n_less + n_nearby/2) & info_neut$maf_rank >= n_less - n_nearby/2]
  if (length(tmp_fst) < n_nearby) {tmp_fst <- info_neut$L13.fst[info_neut$maf_rank >= length(info_neut$maf_rank) - n_nearby]}
  
  tmp_percentiles <- ecdf(tmp_fst)
  tmp_p <- 1-tmp_percentiles(input_fst)
  return(tmp_p)
}
pvals$L13.pval <- do.call("c",mclapply(1:nrow(info), FUN = calc_pval_L13))
calc_pval_L12 <- function(site) {
  input_maf <- info$maf.ML[site]; input_fst <- info$L12.fst[site]
  n_less <- sum(info_neut$maf.ML <= input_maf)
  
  tmp_fst <- info_neut$L12.fst[info_neut$maf_rank <= (n_less + n_nearby/2) & info_neut$maf_rank >= n_less - n_nearby/2]
  if (length(tmp_fst) < n_nearby) {tmp_fst <- info_neut$L12.fst[info_neut$maf_rank >= length(info_neut$maf_rank) - n_nearby]}
  
  tmp_percentiles <- ecdf(tmp_fst)
  tmp_p <- 1-tmp_percentiles(input_fst)
  return(tmp_p)
}
pvals$L12.pval <- do.call("c",mclapply(1:nrow(info), FUN = calc_pval_L12))

calc_pval_D13 <- function(site) {
  input_maf <- info$maf.ML[site]; input_fst <- info$D13.fst[site]
  n_less <- sum(info_neut$maf.ML <= input_maf)
  
  tmp_fst <- info_neut$D13.fst[info_neut$maf_rank <= (n_less + n_nearby/2) & info_neut$maf_rank >= n_less - n_nearby/2]
  if (length(tmp_fst) < n_nearby) {tmp_fst <- info_neut$D13.fst[info_neut$maf_rank >= length(info_neut$maf_rank) - n_nearby]}
  
  tmp_percentiles <- ecdf(tmp_fst)
  tmp_p <- 1-tmp_percentiles(input_fst)
  return(tmp_p)
}
pvals$D13.pval <- do.call("c",mclapply(1:nrow(info), FUN = calc_pval_D13))

rm(list=ls(pattern="calc_pval_"))
rm(drop_samples_missing, design, design2, info_neut)
########


pvals <- pvals[,-c(1:4)]
info <- merge(pvals, info, by="site"); rm(pvals)
########


# save.image("../DATA/pvalues_200neutralsites_pseudohaploid.RData")

save.image("../DATA/pvalues.pseudohaploid.n_neutral_250.RData")


