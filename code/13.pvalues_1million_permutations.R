library(data.table); library(parallel)

########## Read in frequencies without permutation ########## 
load("./fst_estimates.RData")
info <- info[(info$rsid_number == 1),]
info <- info[!duplicated(info$site),] 
info <- info[info$maf.ML >= 0.05,]
##########  Rename info ########## 
info -> info_real; rm(info, design, design2, min_maf, min_n, drop_samples_missing)
info_real <- as.data.frame(info_real)

info_real <- info_real[sign(info_real$delta_L13.ML) == sign(info_real$delta_D13.ML) & sign(info_real$delta_L13.ML) != sign(info_real$delta_L12.ML),]
table(info_real$type)

method="ML"

get_ps <- function(tmp_n) {
  tmp_snp <- info_real$site[tmp_n]
  tmp_file <- paste0("./perms.",tmp_snp,".txt.gz")
  
  res <- info_real[info_real$site == tmp_snp,]
  
  if (file.exists(tmp_file)) {
    extra_perms <- fread(tmp_file)
    
    info_real$number.neutralML <- nrow(extra_perms)
    
    res$double_sign_during.pval.ML <- sum(
      extra_perms$site == tmp_snp & 
        extra_perms$london.post.pre.fst >= res$london.post.pre.fst.ML &
        extra_perms$denmark.post.pre.fst >= res$denmark.post.pre.fst.ML &
        sign(extra_perms$delta_L13) == sign(extra_perms$delta_D13) & 
        sign(extra_perms$delta_L13) != sign(extra_perms$delta_L12)
    )/sum(extra_perms$site == tmp_snp)

    # Condition for appropriate changes
    extra_perms <- extra_perms[sign(extra_perms$delta_L13) == sign(extra_perms$delta_D13) & 
                                 sign(extra_perms$delta_L13) != sign(extra_perms$delta_L12),]
    res$number.perms.conditioned <- nrow(extra_perms)
    
    res$L_1v3.pval.ML.conditioned <- sum(
      extra_perms$london.post.pre.fst >= res$london.post.pre.fst.ML
    )/sum(extra_perms$site == tmp_snp)
    res$D_1v3.pval.ML.conditioned <- sum(
      extra_perms$denmark.post.pre.fst >= res$denmark.post.pre.fst.ML
    )/sum(extra_perms$site == tmp_snp)
    res$L_12v3.pval.ML.conditioned <- sum(
      extra_perms$london.combined.post.fst >= res$london.combined.post.fst.ML
    )/sum(extra_perms$site == tmp_snp)
    
    res$abs_L_1v3.pval.ML.conditioned <- sum(
      abs(extra_perms$delta_L13) >= abs(res$delta_L13.ML)
    )/sum(extra_perms$site == tmp_snp)
    res$abs_D_1v3.pval.ML.conditioned <- sum(
      abs(extra_perms$delta_D13) >= abs(res$delta_D13.ML)
    )/sum(extra_perms$site == tmp_snp)
    res$abs_L_12v3.pval.ML.conditioned <- sum(
      abs(extra_perms$delta_L123) >= abs(res$delta_L123.ML)
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

save.image("./pvalues.RData")
