library(ggplot2); library(patchwork)

method="ML"
load("./perms/permuted_ITERATION/pvals.RData")

# let's not consider the neutral sites for enrichment
info_candidate <- info[!info$type == "neut",]
info_neutral <- info[info$type == "neut",]

  ######### enrichment by bin #########
  tmp_bin_breaks <- c(0.1,.2,.3,.4,0.5)
  tmp_bins <-  matrix(ncol=2, nrow=length(tmp_bin_breaks)); tmp_bins[,1] <- c(0,tmp_bin_breaks[-length(tmp_bin_breaks)]); tmp_bins[,2] <- c(tmp_bin_breaks); rm(tmp_bin_breaks)
  res <- NULL
  for (tmp_enrich in c(0.005, seq(0.01,0.1,0.01), seq(0.1,0.2,0.05))) {
    for (tmp_bin in 1:nrow(tmp_bins)) {
      for (tmp_pop in c("L13", "D13", "L12")) {
        tmp_data <- info_candidate[info_candidate[[paste0("maf.",method)]] > tmp_bins[tmp_bin,1] & info_candidate[[paste0("maf.",method)]] <= tmp_bins[tmp_bin,2],]
        tmp_res <- as.data.frame(matrix(ncol=6, nrow=1)); colnames(tmp_res) <- c('maf_low', 'maf_high', 'population', 'enrichment', 'observed', 'expected')
        tmp_res$population <- tmp_pop
        tmp_res$maf_low <- tmp_bins[tmp_bin,1]; tmp_res$maf_high <- tmp_bins[tmp_bin,2]
        tmp_res$enrichment <- tmp_enrich
        if (tmp_pop == "L13") {tmp_res$observed <- sum(tmp_data[[paste0("L13.pval.",method)]] < tmp_enrich, na.rm=T)}
        if (tmp_pop == "D13") {tmp_res$observed <- sum(tmp_data[[paste0("D13.pval.",method)]] < tmp_enrich, na.rm=T)}
        if (tmp_pop == "L12") {tmp_res$observed <- sum(tmp_data[[paste0("L12.pval.",method)]] < tmp_enrich, na.rm=T)}
        tmp_res$expected <- nrow(tmp_data)*tmp_enrich
        rbind(res, tmp_res) -> res; rm(tmp_res, tmp_data)
      }; rm(tmp_pop)
    }; rm(tmp_bin)
  }; rm(tmp_enrich)

  res$maf_bin <- paste(100*res$maf_low, "% to ", 100*res$maf_high, "%", sep="")
  res$fc <- log2(res$observed/res$expected)

  ###### summary by binomial test #########
  res_binom <- as.data.frame(matrix(ncol=8, nrow=1))
  colnames(res_binom) <- c("enr", "maf", "fc_L13", "p_L13", "fc_L12", "p_L12", "fc_D13", "p_D13")
  res_binom$enr <- c(0.01)
  res_binom$maf <- "all"
  res_binom$method <- method

  tmp <- res[res$population == "L13" ,]
  res_binom$fc_L13[1] <- binom.test(sum(tmp$observed), n = (1/.01)*sum(tmp$expected), p=0.01)$estimate/0.01
  res_binom$p_L13[1] <- binom.test(sum(tmp$observed), n = (1/.01)*sum(tmp$expected), p=0.01)$p.value
  tmp <- res[res$population == "L12" ,]
  res_binom$fc_L12[1] <- binom.test(sum(tmp$observed), n = (1/.01)*sum(tmp$expected), p=0.01)$estimate/0.01
  res_binom$p_L12[1] <- binom.test(sum(tmp$observed), n = (1/.01)*sum(tmp$expected), p=0.01)$p.value
  tmp <- res[res$population == "D13" ,]
  res_binom$fc_D13[1] <- binom.test(sum(tmp$observed), n = (1/.01)*sum(tmp$expected), p=0.01)$estimate/0.01
  res_binom$p_D13[1] <- binom.test(sum(tmp$observed), n = (1/.01)*sum(tmp$expected), p=0.01)$p.value

  res_binom$method <- method
  res_binom$iteration <- ITERATION

rm(method, min_n, min_maf)

write.table(res_binom, "enrichment_res/ITERATION.txt", row.names=F, col.names=F, sep="\t", quote=F)
