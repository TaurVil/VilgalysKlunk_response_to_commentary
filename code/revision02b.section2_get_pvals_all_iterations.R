
library(data.table); library(ggplot2); library(patchwork)

load("./perms/permuted_ITERATION/fst_allsites.with_during.RData")
info -> info_all

for (method in c("bins", "sliding")) {
  if (method == "bins") (ns <- 200)
  if (method == "sliding") (ns <- c(100, 200, 300, 2000))
  
  for (n_nearby in ns) {
    for (min_maf in c(0.05, 0.1)) {
      for (rsid_status in c("biobank", "1kg")) {
        
        ######## Split off neutral sites, then filter candidates based on minor allele frequency ########
        info_neut <- info_all[info_all$type == "neut",]; info_neut$maf_rank <- rank(info_neut$maf.ML)
        info <- info_all[info_all$maf.ML >= min_maf,]
        
        ## filter neutral sites as well if we're comparing all candidate vs all neutral
        ## not introducing a floor where low MAF sites are only compared against higher maf sites
        if (n_nearby == 2000) {info_neut <- info_neut[info_neut$maf.ML > min_maf,]}
        
        ########
        
        ######## keep sites in the rsid set ########
        if (rsid_status == "biobank") {
          sites <- fread("./UKB_variants.txt")
          info_neut <- info_neut[paste(gsub("chr", "", info_neut$site), info_neut$ref, info_neut$alt, sep="_") %in% sites$V1,]
          info <- info[paste(gsub("chr", "", info$site), info$ref, info$alt, sep="_") %in% sites$V1,]
        }
        if (rsid_status == "1kg") {
          targets <- read.delim("./DATA/sites_for_genotyping.txt", header=F)
          colnames(targets) <- c("chr", "loc", "ref", "alt")
          targets$site <- paste(targets$chr, targets$loc, sep="_")
          info_neut <- info_neut[info_neut$site %in% targets$site,]
          info <- info[info$site %in% targets$site,]
          rm(targets)
        }
        ########
        
        ######## Set bin sizes if using bins ##########
        if (method == "bins") {
          mafs <- seq(0.005,0.40, 0.01); window_size <- 0.05*rep(1, length(mafs))
          window_size[mafs>=0.135] <- 0.06; window_size[mafs>=0.145] <- 0.07; window_size[mafs>=0.155] <- 0.08
          window_size[mafs>=0.165] <- 0.09; window_size[mafs>=0.175] <- 0.11; window_size[mafs>=0.185] <- 0.12
          window_size[mafs>=0.205] <- 0.13; window_size[mafs>=0.225] <- 0.15; window_size[mafs>=0.255] <- 0.16
          window_size[mafs>=0.275] <- 0.17; window_size[mafs>=0.305] <- 0.18; window_size[mafs>=0.335] <- 0.19
          window_size[mafs>=0.345] <- 0.2; window_size[mafs>=0.365] <- 0.21; window_size[mafs>=0.375] <- 0.22
          window_size[mafs>=0.395] <- 0.24
          bins <- cbind(mafs-0.005, mafs+0.005, mafs-window_size/2, mafs+window_size/2)
          bins[nrow(bins),2] <- 0.5
        }
        ##########
        
        ######## calculate p-values ########
        library(parallel); pvals <- info[,c(1:5,21,29,61:63)]
        
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
        
        write.table(pvals, paste0("./perms/permuted_ITERATION/pvalues", rsid_status, min_maf, method, n_nearby, "txt", sep="."), row.names = F, col.names=T, sep="\t", quote=F)
        
      }
    }
  }
  
}

