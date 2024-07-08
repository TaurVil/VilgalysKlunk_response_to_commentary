library(data.table); library(ggplot2); library(patchwork); library(parallel)

n_nearby <- 250; min_maf <- 0.05; rsid_status <- 'keep_good_only'
num_lowfreq=1

get_perm <- function(fff) {
  
  load("../DATA/fst_estimates.RData")
  
  info <- as.data.frame(info)
  info <- info[rowSums(info[,which(colnames(info) %like% "alternate.ML")] < 0.01) <= num_lowfreq,]
  
  ######## Randomly sample neutral sites, then split them off and filter candidates by MAF ########
  info_neut <- info[sample(size = sum(info$type == "neut"), x = 1:length(info$type), replace = F),]
  info <- subset(info, ! info$site %in% info_neut$site)
  info_neut$maf_rank <- rank(info_neut$maf.ML)
  
  info <- info[info$maf.ML >= min_maf,]
  ########
  
  pvals <- info[,1:5]
  ########
  # calculate p-values
  ########
  library(parallel)
  
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
  ########
  pvals <- pvals[,-c(1:4)]
  colnames(pvals)[2:4] <- paste(colnames(pvals)[2:4],"ML", sep=".")
  info <- merge(pvals, info, by="site")
  rm(pvals)
  ########???
  
  info$perm <- fff
  info$test <- paste0("ML_n",n_nearby)
  keep <- c(1,5:8,68:ncol(info),which(colnames(info) == "type"),which(colnames(info) %like% "ML"))
  info <- info[,keep]
  info <- info[,c(1:12,21,29,32,35,36:38)]
  rm(keep)
  
  # info_perm_neutral <- rbind(info_perm_neutral, info)
  return(info)
  print(fff)
  rm(info)
}

set.seed(42)
library(parallel) # split up so it runs without crashing anything and we have intermediate files saved
fff <- 1:1000; info_perm_neutral1 <- do.call("rbind", mclapply(fff, get_perm))
fff <- 1001:2000; info_perm_neutral2 <- do.call("rbind", mclapply(fff, get_perm)); save.image("../DATA/fst_estimates_perm_sites.RData")
fff <- 2001:3000; info_perm_neutral3 <- do.call("rbind", mclapply(fff, get_perm))
fff <- 3001:4000; info_perm_neutral4 <- do.call("rbind", mclapply(fff, get_perm)); save.image("../DATA/fst_estimates_perm_sites.RData")
fff <- 4001:5000; info_perm_neutral5 <- do.call("rbind", mclapply(fff, get_perm))
fff <- 5001:6000; info_perm_neutral6 <- do.call("rbind", mclapply(fff, get_perm)); save.image("../DATA/fst_estimates_perm_sites.RData")
fff <- 6001:7000; info_perm_neutral7 <- do.call("rbind", mclapply(fff, get_perm))
fff <- 7001:8000; info_perm_neutral8 <- do.call("rbind", mclapply(fff, get_perm)); save.image("../DATA/fst_estimates_perm_sites.RData")
fff <- 8001:9000; info_perm_neutral9 <- do.call("rbind", mclapply(fff, get_perm))
fff <- 9001:10000; info_perm_neutral10 <- do.call("rbind", mclapply(fff, get_perm)); save.image("../DATA/fst_estimates_perm_sites.RData")
rm(fff)

info_perm_neutral <- rbind(info_perm_neutral1, info_perm_neutral2, info_perm_neutral3, info_perm_neutral4, info_perm_neutral5, 
                           info_perm_neutral6, info_perm_neutral7, info_perm_neutral8, info_perm_neutral9, info_perm_neutral10)
rm(info_perm_neutral1, info_perm_neutral2, info_perm_neutral3, info_perm_neutral4, info_perm_neutral5, 
   info_perm_neutral6, info_perm_neutral7, info_perm_neutral8, info_perm_neutral9, info_perm_neutral10)

save.image("../DATA/fst_estimates_perm_sites.RData")
