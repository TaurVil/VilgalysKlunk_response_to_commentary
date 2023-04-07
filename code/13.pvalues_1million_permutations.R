library(data.table); library(parallel)

########## Read in frequencies without permutation ########## 
load("./fst_estimates.RData")
info <- info[(info$rsid_number == 1),]
info <- info[!duplicated(info$site),] 
info <- info[info$maf.ML >= 0.05,]
##########  Rename info ########## 
info -> info_real; rm(info, design, design2, min_maf, min_n, drop_samples_missing)
info_real <- as.data.frame(info_real)

method="ML"

get_ps <- function(tmp_n) {
  tmp_snp <- info_real$site[tmp_n]
  tmp_file <- paste0("./perms.",tmp_snp,".txt.gz")
  
  res <- info_real[info_real$site == tmp_snp,]
  
  if (file.exists(tmp_file)) {
    extra_perms <- fread(tmp_file)
    
    info_real$number.neutralML <- nrow(extra_perms)
    
    res$double.pval.ML <- sum(
      extra_perms$site == tmp_snp & 
        extra_perms$london.post.pre.fst >= res$london.post.pre.fst.ML &
        extra_perms$denmark.post.pre.fst >= res$denmark.post.pre.fst.ML 
    )/sum(extra_perms$site == tmp_snp)
    
    res$double_plus_sign.pval.ML <- sum(
      extra_perms$site == tmp_snp & 
        extra_perms$london.post.pre.fst >= res$london.post.pre.fst.ML &
        extra_perms$denmark.post.pre.fst >= res$denmark.post.pre.fst.ML &
        sign(extra_perms$delta_L13) == sign(extra_perms$delta_D13)
    )/sum(extra_perms$site == tmp_snp)
    
    res$double_sign_during.pval.ML <- sum(
      extra_perms$site == tmp_snp & 
        extra_perms$london.post.pre.fst >= res$london.post.pre.fst.ML &
        extra_perms$denmark.post.pre.fst >= res$denmark.post.pre.fst.ML &
        sign(extra_perms$delta_L13) == sign(extra_perms$delta_D13) & 
        sign(extra_perms$delta_L13) != sign(extra_perms$delta_L12)
    )/sum(extra_perms$site == tmp_snp)
    
    res$D13.pval.ML <- sum(
      extra_perms$site == tmp_snp & 
        extra_perms$denmark.post.pre.fst >= res$denmark.post.pre.fst.ML
    )/sum(extra_perms$site == tmp_snp)
    
    res$L13.pval.ML <- sum(
      extra_perms$site == tmp_snp & 
        extra_perms$london.post.pre.fst >= res$london.post.pre.fst.ML
    )/sum(extra_perms$site == tmp_snp)
    
    return(res)
  }
}

res <- NULL
for (i in 1:nrow(info_real)) {
  d <- get_ps(i)
  res <- rbind(res,d)
  print(i); rm(d)
}

save.image("./pvalues_1mil_perms.RData")
