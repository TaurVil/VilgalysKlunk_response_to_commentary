# module load plink 


pacman::p_load(vroom, tidyverse, ggforestplot, TwoSampleMR, patchwork, gt, gwasvcf, gwasglue)
set_bcftools('/project/lbarreiro/USERS/tauras/Programs/bcftools/bcftools/bcftools')

library(ggplot2); library(plyr); library(patchwork)
load("cyto_results.RData")
load("cfu_results.RData")

## Create table of all test statistics
# colnames(res_cyto)[2] <- c("trait") # (cytokine or bacterial clearance)
# res_cfu$trait <- "bacterial_clearance"
# table_s1 <- rbind(res_cfu[,c(1,5,2:3)],res_cyto[,1:4])
# write.table(table_s1, "./Table_S1.txt", row.names=F, col.names=T, sep="\t", quote=F)

# Let's try a meta-analysis of the cytokine data
res_cfu$p_meta <- res_cfu$num_cyto <- NA
for (site in (unique(res_cyto$snp))) {
  tmp <- res_cyto[res_cyto$snp == site,]
  xi <- -2*sum(log(tmp$pval))
  res_cfu$p_meta[res_cfu$snp == site] <- 1-pchisq(xi, 2*nrow(tmp))
  rm(xi)
  
  res_cfu$num_cyto[res_cfu$snp == site] <- sum(tmp$pval < 0.05)
  
  rm(tmp)
}; rm(site)

r <- res_cfu[!is.na(res_cfu$p_meta) & res_cfu$snp %in% c("chr5_96235896", info$site[info$type != "neut"]),]

head(r)
sum(r$num_cyto >= r$num_cyto[r$snp == "chr5_96235896"])
sum(r$num_cyto >= r$num_cyto[r$snp == "chr5_96235896"])/nrow(r)
sum(r$pval <= r$pval[r$snp == "chr5_96235896"])
sum(r$pval <= r$pval[r$snp == "chr5_96235896"])/nrow(r)
sum(r$pval <= 0.0155)
sum(r$pval <= 0.0155)/nrow(r)

# Correlation between CFU data and number of cytokines
cor.test(-log10(r$pval), r$num_cyto, method="spearman")
cor.test(-log10(r$pval), r$num_cyto, method="pearson")

# number of loci with similar evidence 
sum(r$pval[r$num_cyto >= 4] <= 0.05)
r[r$pval <= 0.05 & r$num_cyto >= 4,]

table1_snps <- r[r$pval <= 0.05 & r$num_cyto >= 4,c(1:3,5)]
table1_snps$estimate <- abs(table1_snps$estimate)
colnames(table1_snps)[2:3] <- c("estimate_bacterialclearance", "p_bacterialclearance")

table1_snps$rsid <- c("rs2248374", "rs12949531", "rs9380739", "rs9394492", "rs6651252")
table1_snps$notes <- c("candidate and splice variant", "", "in linkage with rs9394492", "in linkage with rs9380739", "")

table1_snps <- table1_snps[,c(1,5:6,2:4)]

for (iii in table1_snps$snp) {
  t <- res_cyto[res_cyto$snp == iii,]
  for (j in t$cytokine) {
    table1_snps[[paste0("estimate_", j)]][table1_snps$snp == iii] <- t$estimate[t$cytokine == j]
    table1_snps[[paste0("p_", j)]][table1_snps$snp == iii] <- t$pval[t$cytokine == j]
  }; rm(j,t)
}; rm(iii)


table1_snps$rsid[table1_snps$rsid == "rs2248374"] <- "rs2549794"

# Pneumonia (death)
expd2 <- gwasvcf::query_gwas("ieu-b-4979.vcf.gz", rsid=c(table1_snps$rsid, "rs2549794"))
expd3 <- gwasglue::gwasvcf_to_TwoSampleMR(expd2, type="exposure")
tmp <- expd3[,c(11,5,7)]
colnames(tmp) <- c("rsid", "estimate_pneumonia (death)", "p_pneumonia (death)")
table1_snps <- merge(table1_snps, tmp, by="rsid")

# Pneumonia (death in critical care)
expd2 <- gwasvcf::query_gwas("ieu-b-4977.vcf.gz", rsid=c(table1_snps$rsid, "rs2549794"))
expd3 <- gwasglue::gwasvcf_to_TwoSampleMR(expd2, type="exposure")
tmp <- expd3[,c(11,5,7)]
colnames(tmp) <- c("rsid", "estimate_pneumonia (death critical care)", "p_pneumonia (death critical care)")
table1_snps <- merge(table1_snps, tmp, by="rsid")

# Pneumonia (critical care)
expd2 <- gwasvcf::query_gwas("ieu-b-4978.vcf.gz", rsid=c(table1_snps$rsid, "rs2549794"))
expd3 <- gwasglue::gwasvcf_to_TwoSampleMR(expd2, type="exposure")
tmp <- expd3[,c(11,5,7)]
colnames(tmp) <- c("rsid", "estimate_pneumonia (critical care)", "p_pneumonia (critical care)")
table1_snps <- merge(table1_snps, tmp, by="rsid")

# Pneumonia
expd2 <- gwasvcf::query_gwas("ieu-b-4976.vcf.gz", rsid=c(table1_snps$rsid, "rs2549794"))
expd3 <- gwasglue::gwasvcf_to_TwoSampleMR(expd2, type="exposure")
tmp <- expd3[,c(11,5,7)]
colnames(tmp) <- c("rsid", "estimate_pneumonia", "p_pneumonia")
table1_snps <- merge(table1_snps, tmp, by="rsid")


table1_snps$rsid[table1_snps$rsid == "rs2549794"] <- "rs2549794/rs2248374"
write.table(table1_snps, "./Table_S1.txt", row.names=F, col.names=T, sep="\t", quote=F)


panel1 <- ggplot(data=d, aes(x=Freq)) + 
  geom_bar(fill="darkgray") + 
  theme_classic() + xlab("number of cytokines") + ylab("number of variants") + 
  ggtitle("number of cytokines (p < 0.05)") +
  theme(text=element_text(size=14), axis.text = element_text(color="black"))  ; panel1
  

panel2 <- ggplot(data=r, aes(x=-log10(p_meta))) +
  geom_density() + 
  geom_vline(xintercept=abs(r$p_meta[res_cfu$snp == "chr5_96235896"]), col='red') +
  # geom_hline(yintercept=0, col='gray15') +
  coord_cartesian(ylim=c(0,2)) +
  theme_classic() + ggtitle("multivariate p-value for all 10 cytokines") +
  ylab("density") +
  expand_limits(x = 0, y = 0) + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
  theme(text=element_text(size=14), axis.text = element_text(color="black"))

panel3 <- ggplot(data=res_cfu, aes(x=abs(estimate))) +
  geom_density() + 
  geom_vline(xintercept=abs(res_cfu$estimate[res_cfu$snp == "chr5_96235896"]), col='red') +
  # geom_hline(yintercept=0, col='gray15') +
  coord_cartesian(ylim=c(0,4.5)) +
  theme_classic() + ggtitle("correlation coefficient for bacterial clearance") +
  ylab("density") +
  expand_limits(x = 0, y = 0) + scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 0)) + 
  theme(text=element_text(size=14), axis.text = element_text(color="black"))


panel1 + panel2 + panel3


sum(r$p_meta <= r$p_meta[r$snp == "chr5_96235896"]) # top 1.96% of hits
sum(r$p_meta <= r$p_meta[r$snp == "chr5_96235896"])/nrow(r)

sum(r$p_meta <= r$p_meta[r$snp == "chr5_96244549"]) # top 11% of hits
sum(r$p_meta <= r$p_meta[r$snp == "chr5_96244549"])/nrow(r)


r$r_cfu <- rank(r$pval)
r$r_cyto <- rank(r$p_meta)

r[r$snp == "chr5_96235896",]
View(r[r$r_cfu < 50 & r$r_cyto < 50,])





sum(r$num_cyto >= 3)
sum(r$num_cyto >= 3 & r$pval < 0.05)


sum(r$pval[r$num_cyto >= 4] <= 0.0155)/length(r$pval[r$num_cyto >= 4])

sum(r$pval[r$num_cyto >= 3] <= 0.0155)
sum(r$pval[r$num_cyto >= 3] <= 0.0155)/length(r$pval[r$num_cyto >= 3])


d <- as.data.frame(table(res_cyto$snp[res_cyto$pval > 0.05 & res_cyto$snp %in% c("chr5_96235896", info$site[info$type != "neut"])]))
colnames(d) <- c("site", "Freq")
d$Freq <- 10-d$Freq
