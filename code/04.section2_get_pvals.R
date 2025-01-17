
library(data.table); library(ggplot2); library(parallel); library(patchwork)

load("./perms/permuted_ITERATION/fst.RData")
num_lowfreq=1 # maximum number of population-time points which we'll allow to not observe the alternate allele

info <- as.data.frame(info)
info <- info[rowSums(info[,which(colnames(info) %like% "alternate.ML")] < 0.01) <= num_lowfreq,]

info -> info_all

for (method in c("sliding")) {
  if (method == "bins") (ns <- 200)
  if (method == "sliding") (ns <- c(50, 100, 150, 200, 250, 300, 400, 2000))
  
  for (n_nearby in ns) {
    for (min_maf in c(0.05)) {
      for (rsid_status in "1kg") {
        
        ######## Split off neutral sites, then filter candidates based on minor allele frequency ########
        info_neut <- info_all[info_all$type == "neut",]; info_neut$maf_rank <- rank(info_neut$maf.ML)
        info <- info_all[info_all$maf.ML >= min_maf,]
        
        ## filter neutral sites as well if we're comparing all candidate vs all neutral
        ## not introducing a floor where low MAF sites are only compared against higher maf sites
        if (n_nearby == 2000) {info_neut <- info_neut[info_neut$maf.ML > min_maf,]}
        
        ########
        
        ######## keep sites in the rsid set ########
        # already done in the first step
        ########
        
        ######## Set bin sizes if using bins ##########
        # Bins were removed in favor of sliding windows. For an example, see code for the original publication. 
        ##########
        
        ######## calculate p-values ########
        pvals <- info[,c(1:5,21,29,61:63)]
        
        get_tmp_fst_table <- function(input_maf) {
          if (method == "sliding") {
            n_less <- sum(info_neut$maf.ML <= input_maf)
            tmp_fst <- info_neut[info_neut$maf_rank <= (n_less + n_nearby/2) & info_neut$maf_rank >= n_less - n_nearby/2,]
            if (nrow(tmp_fst) < n_nearby) {tmp_fst <- info_neut[info_neut$maf_rank >= length(info_neut$maf_rank) - n_nearby,]}
          }
          if (method == "bins") {
            tmp_bin <- bins[input_maf > bins[,1] & input_maf <= bins[,2],]
            tmp_fst <- info_neut[info_neut$maf.ML > tmp_bin[3] & info_neut$maf.ML <= tmp_bin[4]]
          }
          return(tmp_fst)
        }
        
        calc_pvals <- function(site, test) {
          input_maf <- info$maf.ML[site]
          tmp_fst <- get_tmp_fst_table(input_maf)
          
          if (test == "L13") {
            input_fst <- info$london.post.pre.fst.ML[site]
            tmp_percentiles <- ecdf(tmp_fst$london.post.pre.fst.ML)
          }
          if (test == "L12") {
            input_fst <- info$london.during.pre.fst.ML[site]
            tmp_percentiles <- ecdf(tmp_fst$london.during.pre.fst.ML)
          }
          if (test == "D13") {
            input_fst <- info$denmark.post.pre.fst.ML[site]
            tmp_percentiles <- ecdf(tmp_fst$denmark.post.pre.fst.ML)
          }
          
          tmp_p <- 1-tmp_percentiles(input_fst)
          return(tmp_p)
        }
        
        pvals$L13.pval <- do.call("c",mclapply(1:nrow(info), FUN = calc_pvals, test="L13"))
        pvals$L12.pval <- do.call("c",mclapply(1:nrow(info), FUN = calc_pvals, test="L12"))
        pvals$D13.pval <- do.call("c",mclapply(1:nrow(info), FUN = calc_pvals, test="D13"))
        ########

        write.table(pvals, paste("./perms/permuted_ITERATION/pvalues.with_during", rsid_status, min_maf, method, n_nearby, "txt", sep="."), row.names = F, col.names=T, sep="\t", quote=F)
        
      }
    }
  }
}


