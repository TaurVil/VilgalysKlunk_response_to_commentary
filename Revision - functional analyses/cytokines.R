load(("../DATA/pvalues_200neutralsites.RData"))
info <- info[,c(5:8, 24, 1)]
rm(min_maf, min_n, n_nearby, rsid_status)

## Read in data
# each row contains the cytokine level for pestis-stimulated samples versus the null control samples at the same time point. we'll calculate the difference between them below, and call it `value`
###########################
data <- read.delim("./DATA/cytokines.txt")
data$value <- data$Pestis - data$Null
# Effects become more obvious over time. We therefore focus on the 24h time point rather than the 5h timepoint. 
data <- subset(data, data$time == '24h') 
# get the list of cytokines analyzed
cks <- unique(data$cytokine)
###########################

# Genotype data
###########################
gt <- read.delim("./DATA/genotypes_for_functional_analyses.txt")
colnames(gt) <- gsub("X","", colnames(gt))

gt <- t(gt); data$gt <- NA

# correlation in ERAP2 genotypes
tmp_gt <- as.data.frame(gt[, which(colnames(gt) == "chr5_96244549")])
colnames(tmp_gt) <- "gt"
tmp_gt$name <- row.names(tmp_gt)
tmp <- data[data$cytokine == "CCL3",]
tmp <- merge(tmp[,-1], by.x="sample", tmp_gt, by.y="name")
plot(tmp$gt ~ factor(tmp$genotype))
table(tmp$gt, factor(tmp$genotype))
cor.test(tmp$gt, as.numeric(factor(tmp$genotype, levels=c("Allele 1", "HE", "Allele 2"))))

tmp_gt <- as.data.frame(gt[, which(colnames(gt) == "chr5_96235896")])
colnames(tmp_gt) <- "gt"
tmp_gt$name <- row.names(tmp_gt)
tmp <- data[data$cytokine == "CCL3",]
tmp <- merge(tmp[,-1], by.x="sample", tmp_gt, by.y="name")
plot(tmp$gt ~ factor(tmp$genotype))
table(tmp$gt, factor(tmp$genotype))
cor.test(tmp$gt, as.numeric(factor(tmp$genotype, levels=c("Allele 1", "HE", "Allele 2"))))

###########################


## For each cytokine, measure cytokine level as a function of genotype and plot
## Normalize data for each cytokine using qqnorm
###########################


res <- NULL
for (site in colnames(gt)) {
  
  tmp_gt <- as.data.frame(gt[, which(colnames(gt) == site)])
  colnames(tmp_gt) <- "gt"
  tmp_gt$name <- row.names(tmp_gt)
  
  if (length(unique(tmp_gt$gt)) > 1) {
    for (cyto in cks) {
      tmp <- data[data$cytokine == cyto,]
      tmp <- merge(tmp[,-1], by.x="sample", tmp_gt, by.y="name")
      tmp$value <- qqnorm(tmp$value, plot.it = F)$x
      
      tmp_res <- as.data.frame(matrix(ncol=6, nrow=1))
      colnames(tmp_res) <- c("snp", "cytokine", "estimate", "pval", "r.squared", "method")
      tmp_res$snp <- site
      tmp_res$cytokine <- cyto
      tmp_res$method <- "qqnorm"
      
      if (length(unique(tmp$gt)) > 1) {
        model <- lm(tmp$value ~ tmp$gt)
        
        tmp_res$estimate <- summary(model)$coefficients[2,1]
        tmp_res$pval <- summary(model)$coefficients[2,4]
        tmp_res$r.squared <- summary(model)$r.squared
        
        rbind(res, tmp_res) -> res
        rm(model)
      }
      rm(tmp_res, tmp)
    }; rm(cyto)
  }
  rm(tmp_gt)
}; rm(site)


res -> res_cyto; rm(res)

save.image("./cyto_results.RData")


res_cyto[res_cyto$snp == "chr5_96235896",]
res_cyto[res_cyto$snp == "chr5_96244549",]

table(table(res_cyto$snp[res_cyto$pval < 0.05 & res_cyto$snp %in% info$site]))
barplot(table(table(res_cyto$snp[res_cyto$pval < 0.05 & res_cyto$snp %in% info$site])), xlab="number of cytokines with p < 0.05", ylab="number of sites")

