library(data.table); library(ggplot2); library(patchwork); library(parallel)

# n_nearby sets the number of neutral sites
# num_lowfreq sets the number of populations allowed to not contain the alternate allele

for (n_nearby in c(50, 100, 150, 200, 250, 300, 400, 2000)) { # 50, 100, 2000
  for (num_lowfreq in 1) {
  load("../DATA/fst_estimates.RData")
  
  min_maf <- 0.05; rsid_status <- 'keep_good_only'

  info <- as.data.frame(info)
  info <- info[rowSums(info[,which(colnames(info) %like% "alternate.ML")] < 0.01) <= num_lowfreq,]
  
  # report correlation between ML and GL allele frequencies
  sqrt(summary(lm(info$maf.ML ~ info$maf.GL))$r.squared)
  
  ######## Split off neutral sites, then filter candidates based on minor allele frequency ########
  info_neut <- info[info$type == "neut",]
  info_neut$maf_rank <- rank(info_neut$maf.ML)
  info <- info[info$maf.ML >= min_maf,]
  ########
  
  ######## calculate p-values ########
  pvals <- info[,1:5]
  
  calc_pval_L13 <- function(site) {
    input_maf <- info$maf.ML[site]
    input_fst <- info$london.post.pre.fst.ML[site]
    n_less <- sum(info_neut$maf.ML <= input_maf)
    
    tmp_fst <- info_neut$london.post.pre.fst.ML[info_neut$maf_rank <= (n_less + n_nearby/2) & info_neut$maf_rank >= n_less - n_nearby/2]
    if (length(tmp_fst) < n_nearby) {tmp_fst <- info_neut$london.post.pre.fst.ML[info_neut$maf_rank >= length(info_neut$maf_rank) - n_nearby]}
    
    tmp_percentiles <- ecdf(tmp_fst)
    tmp_p <- 1-tmp_percentiles(input_fst)
    return(tmp_p)
  }
  pvals$L13.pval <- do.call("c",mclapply(1:nrow(info), FUN = calc_pval_L13))
  calc_pval_L12 <- function(site) {
    input_maf <- info$maf.ML[site]
    input_fst <- info$london.during.pre.fst.ML[site]
    n_less <- sum(info_neut$maf.ML <= input_maf)
    
    tmp_fst <- info_neut$london.during.pre.fst.ML[info_neut$maf_rank <= (n_less + n_nearby/2) & info_neut$maf_rank >= n_less - n_nearby/2]
    if (length(tmp_fst) < n_nearby) {tmp_fst <- info_neut$london.during.pre.fst.ML[info_neut$maf_rank >= length(info_neut$maf_rank) - n_nearby]}
    
    tmp_percentiles <- ecdf(tmp_fst)
    tmp_p <- 1-tmp_percentiles(input_fst)
    return(tmp_p)
  }
  pvals$L12.pval <- do.call("c",mclapply(1:nrow(info), FUN = calc_pval_L12))
  
  calc_pval_D13 <- function(site) {
    input_maf <- info$maf.ML[site]
    input_fst <- info$denmark.post.pre.fst.ML[site]
    n_less <- sum(info_neut$maf.ML <= input_maf)
    
    tmp_fst <- info_neut$denmark.post.pre.fst.ML[info_neut$maf_rank <= (n_less + n_nearby/2) & info_neut$maf_rank >= n_less - n_nearby/2]
    if (length(tmp_fst) < n_nearby) {tmp_fst <- info_neut$denmark.post.pre.fst.ML[info_neut$maf_rank >= length(info_neut$maf_rank) - n_nearby]}
    
    tmp_percentiles <- ecdf(tmp_fst)
    tmp_p <- 1-tmp_percentiles(input_fst)
    return(tmp_p)
  }
  pvals$D13.pval <- do.call("c",mclapply(1:nrow(info), FUN = calc_pval_D13))
  
  rm(list=ls(pattern="calc_pval_"))
  rm(drop_samples_missing, design, design2, info_neut)
  ########
  
  pvals <- pvals[,-c(1:4)]; colnames(pvals)[2:4] <- paste(colnames(pvals)[2:4],"ML", sep=".")
  info <- merge(pvals, info, by="site"); rm(pvals)
  ########
  
  candidates <- info$rsid[info$L13.pval.ML < 0.05 & info$D13.pval.ML < 0.1 & 
                       sign(info$delta_L13.ML) == sign(info$delta_D13.ML) &
                       sign(info$delta_L13.ML) != sign(info$delta_L12.ML)]
  
  save.image(paste0("../DATA/pvalues.n_neutral_",n_nearby, ".max_low_freq_pops_", num_lowfreq,".RData"))
  
  print(paste0("looking at ", n_nearby, " neutral sites"))
  print(candidates)
}}

# info[info$rsid %like% "rs2549794",c(1:4,64:66,68)]
