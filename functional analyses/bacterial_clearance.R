cfu <- read.delim("./DATA/cfu_data.txt")
cfu$perc_kill <- 1-cfu$X24mean/cfu$X2mean
cfu <- cfu[cfu$perc_kill > 0,]

# plot(cfu$perc_kill ~ cfu$X2mean, pch=16, cex=3, xlab="bacteria at 2h", ylab="percent of bacteria killed")
cfu$id <- as.factor(cfu$id)

cfu$Genotype <- NA
cfu$Genotype[cfu$allele == "Allele 1"] <- 0
cfu$Genotype[cfu$allele == "Allele 2"] <- 2
cfu$Genotype[cfu$allele == "HE"] <- 1




# Get sites and genotype data
################
load(("../DATA/pvalues_200neutralsites.RData"))
info <- info[,c(5:8, 24, 1)]
rm(min_maf, min_n, n_nearby, rsid_status)
################
# gt <- as.data.frame(t(read.delim("./targets.012", header=F)[,-1]))
# colnames(gt) <- read.delim("targets.012.indv", header=F)$V1
# sites <- read.delim("targets.012.pos", header=F)
# sites$loc <- paste0("chr", sites$V1, "_", sites$V2)
# 
# gt <- gt[!duplicated(sites$loc),]
# sites <- sites[!duplicated(sites$loc), ]
# row.names(gt) <- sites$loc 
# rm(sites)
# gt <- gt[row.names(gt) %in% c(info$site, "chr5_96235896"),]

gt <- read.delim("./DATA/genotypes_for_functional_analyses.txt")
colnames(gt) <- gsub("X","", colnames(gt))

###########################

gt <- t(gt)

res <- NULL
for (site in colnames(gt)) {
  
  tmp_gt <- as.data.frame(gt[, which(colnames(gt) == site)])
  colnames(tmp_gt) <- "gt"
  tmp_gt$name <- row.names(tmp_gt)
  
  if (length(unique(tmp_gt$gt)) > 1) {
    tmp <- merge(cfu[,-1], by.x="id", tmp_gt, by.y="name")
    
    tmp_res <- as.data.frame(matrix(ncol=4, nrow=1))
    colnames(tmp_res) <- c("snp", "estimate", "pval", "method")
    tmp_res$snp <- site
    tmp_res$method <- "raw_ratio"
    
    tmp_res$estimate <- cor.test(tmp$perc_kill, tmp$gt, method="spearman")$estimate
    tmp_res$pval <- cor.test(tmp$perc_kill, tmp$gt, method="spearman")$p.value
    
    rbind(res, tmp_res) -> res
    
    rm(tmp_res, tmp)
  }
  rm(tmp_gt)
}; rm(site)

res -> res_cfu; rm(res)

save.image("./cfu_results.RData")




## ERAP2 results 

res_cfu[res_cfu$snp == "chr5_96235896",]
res_cfu[res_cfu$snp == "chr5_96244549",]

sum(res_cfu$pval <= res_cfu$pval[res_cfu$snp == "chr5_96235896"]) #
sum(res_cfu$pval <= res_cfu$pval[res_cfu$snp == "chr5_96235896"])/nrow(res_cfu)

# original value which is without dropping one sample with failed genome-wide genotype data
sum(res_cfu$pval <= 0.0155) # 
sum(res_cfu$pval <= 0.0155)/nrow(res_cfu)
