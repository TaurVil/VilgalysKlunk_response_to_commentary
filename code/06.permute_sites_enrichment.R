load("../DATA/fst_estimates_perm_sites.RData")

# set bins to consider
tmp_bin_breaks <- c(0.1,.2,.3,.4,0.5)
tmp_bins <-  matrix(ncol=2, nrow=length(tmp_bin_breaks)); tmp_bins[,1] <- c(0,tmp_bin_breaks[-length(tmp_bin_breaks)]); tmp_bins[,2] <- c(tmp_bin_breaks); rm(tmp_bin_breaks)

res_perm <- NULL
for (fff in unique(info_perm_neutral$perm)) {
  info_candidate <- info_perm_neutral[info_perm_neutral$perm == fff,]
  
  res <- NULL
  for (tmp_enrich in c(0.01)) {
    for (tmp_bin in 1:nrow(tmp_bins)) {
      for (tmp_pop in c("L13", "D13", "L12")) {
        tmp_data <- info_candidate[info_candidate$maf.ML > tmp_bins[tmp_bin,1] & info_candidate$maf.ML <= tmp_bins[tmp_bin,2],]
        tmp_res <- as.data.frame(matrix(ncol=6, nrow=1)); colnames(tmp_res) <- c('maf_low', 'maf_high', 'population', 'enrichment', 'observed', 'expected')
        tmp_res$population <- tmp_pop
        tmp_res$maf_low <- tmp_bins[tmp_bin,1]; tmp_res$maf_high <- tmp_bins[tmp_bin,2]
        tmp_res$enrichment <- tmp_enrich  
        if (tmp_pop == "L13") {tmp_res$observed <- sum(tmp_data$L13.pval.ML < tmp_enrich, na.rm=T)}
        if (tmp_pop == "D13") {tmp_res$observed <- sum(tmp_data$D13.pval.ML < tmp_enrich, na.rm=T)}
        if (tmp_pop == "L12") {tmp_res$observed <- sum(tmp_data$L12.pval.ML < tmp_enrich, na.rm=T)}
        tmp_res$expected <- nrow(tmp_data)*tmp_enrich
        rbind(res, tmp_res) -> res; rm(tmp_res, tmp_data)
      }; rm(tmp_pop)
    }; rm(tmp_bin)
  }
  res$maf_bin <- paste(100*res$maf_low, "% to ", 100*res$maf_high, "%", sep="")
  res$fc <- log2(res$observed/res$expected)
  
  tmp2 <- as.data.frame(matrix(ncol=8, nrow=1))
  colnames(tmp2) <- c("enr", "maf", "fc_L13", "p_L13", "fc_L12", "p_L12", "fc_D13", "p_D13")
  tmp2$enr <- c(0.01); tmp2$maf <- "all"
  
  tmp <- res[res$population == "L13" ,]
  tmp2$fc_L13[1] <- binom.test(sum(tmp$observed), n = (1/0.01)*sum(tmp$expected), p=0.01)$estimate/0.01
  tmp2$p_L13[1] <- binom.test(sum(tmp$observed), n = (1/0.01)*sum(tmp$expected), p=0.01)$p.value
  tmp <- res[res$population == "L12" ,]
  tmp2$fc_L12[1] <- binom.test(sum(tmp$observed), n = (1/0.01)*sum(tmp$expected), p=0.01)$estimate/0.01
  tmp2$p_L12[1] <- binom.test(sum(tmp$observed), n = (1/0.01)*sum(tmp$expected), p=0.05)$p.value
  tmp <- res[res$population == "D13" ,]
  tmp2$fc_D13[1] <- binom.test(sum(tmp$observed), n = (1/0.01)*sum(tmp$expected), p=0.01)$estimate/0.01
  tmp2$p_D13[1] <- binom.test(sum(tmp$observed), n = (1/0.01)*sum(tmp$expected), p=0.01)$p.value
  
  res_perm <- rbind(res_perm,tmp2); rm(tmp2, tmp, res)
  print(fff)
}; rm(fff); rm(tmp_enrich)


# add in the real enrichment as well
info_candidate <- info_candidate[info_candidate$type != "neut",]
res <- NULL
for (tmp_enrich in c(0.01)) {
  for (tmp_bin in 1:nrow(tmp_bins)) {
    for (tmp_pop in c("L13", "D13", "L12")) {
      tmp_data <- info_candidate[info_candidate$maf.ML > tmp_bins[tmp_bin,1] & info_candidate$maf.ML <= tmp_bins[tmp_bin,2],]
      tmp_res <- as.data.frame(matrix(ncol=6, nrow=1)); colnames(tmp_res) <- c('maf_low', 'maf_high', 'population', 'enrichment', 'observed', 'expected')
      tmp_res$population <- tmp_pop
      tmp_res$maf_low <- tmp_bins[tmp_bin,1]; tmp_res$maf_high <- tmp_bins[tmp_bin,2]
      tmp_res$enrichment <- tmp_enrich  
      if (tmp_pop == "L13") {tmp_res$observed <- sum(tmp_data$L13.pval.ML < tmp_enrich)}
      if (tmp_pop == "D13") {tmp_res$observed <- sum(tmp_data$D13.pval.ML < tmp_enrich)}
      if (tmp_pop == "L12") {tmp_res$observed <- sum(tmp_data$L12.pval.ML < tmp_enrich)}
      tmp_res$expected <- nrow(tmp_data)*tmp_enrich
      rbind(res, tmp_res) -> res; rm(tmp_res, tmp_data)
    }; rm(tmp_pop)
  }; rm(tmp_bin)
}
res$maf_bin <- paste(100*res$maf_low, "% to ", 100*res$maf_high, "%", sep="")
res$fc <- log2(res$observed/res$expected)

tmp2 <- as.data.frame(matrix(ncol=8, nrow=1))
colnames(tmp2) <- c("enr", "maf", "fc_L13", "p_L13", "fc_L12", "p_L12", "fc_D13", "p_D13")
tmp2$enr <- c(0.01); tmp2$maf <- "all"

tmp <- res[res$population == "L13" & res$enrichment == 0.01,]
tmp2$fc_L13[1] <- binom.test(sum(tmp$observed), n = (1/0.01)*sum(tmp$expected), p=0.01)$estimate/0.01
tmp2$p_L13[1] <- binom.test(sum(tmp$observed), n = (1/0.01)*sum(tmp$expected), p=0.01)$p.value
tmp <- res[res$population == "L12" & res$enrichment == 0.01,]
tmp2$fc_L12[1] <- binom.test(sum(tmp$observed), n = (1/0.01)*sum(tmp$expected), p=0.01)$estimate/0.01
tmp2$p_L12[1] <- binom.test(sum(tmp$observed), n = (1/0.01)*sum(tmp$expected), p=0.01)$p.value
tmp <- res[res$population == "D13" & res$enrichment == 0.01,]
tmp2$fc_D13[1] <- binom.test(sum(tmp$observed), n = (1/0.01)*sum(tmp$expected), p=0.01)$estimate/0.01
tmp2$p_D13[1] <- binom.test(sum(tmp$observed), n = (1/0.01)*sum(tmp$expected), p=0.01)$p.value

res_real <- tmp2; rm(tmp2, tmp, res)

rm(x,y,tmp_fst, rsid_status, n_nearby, min_n, min_maf)
rm(info_candidate, get_perm, info, info_perm_neutral)

save.image("../DATA/perm_sites_results.RData")
