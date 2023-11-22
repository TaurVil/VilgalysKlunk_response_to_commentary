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
