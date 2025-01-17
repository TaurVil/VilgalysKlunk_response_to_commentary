# Originated from FINAL_FIGURE1.R
library(ggplot2); library(patchwork)
cbPalette <- c("#999999", "#D55E00", "#CC79A7", "#56B4E9", "#0072B2")

########  Enrichment for original data ######## 
load("../DATA/pvalues.n_neutral_200.max_low_freq_pops_1.RData")

info[info$rsid == "rs2549794",]

table(info$type[info$L13.pval.ML <= 0.05 & info$D13.pval.ML <= 0.1 & sign(info$delta_L13.ML) == sign(info$delta_D13.ML) & sign(info$delta_L13.ML) != sign(info$delta_L12.ML)])
info$rsid[info$L13.pval.ML <= 0.05 & info$D13.pval.ML <= 0.1 & sign(info$delta_L13.ML) == sign(info$delta_D13.ML) & sign(info$delta_L13.ML) != sign(info$delta_L12.ML)]

method="ML"
# exclude neutral sites
info_candidate <- info[!info$type == "neut",]
info_neutral <- info[info$type == "neut",]
######### by bin #########
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
      if (tmp_pop == "L13") {tmp_res$observed <- sum(tmp_data[[paste0("L13.pval.",method)]] < tmp_enrich)}
      if (tmp_pop == "D13") {tmp_res$observed <- sum(tmp_data[[paste0("D13.pval.",method)]] < tmp_enrich)}
      if (tmp_pop == "L12") {tmp_res$observed <- sum(tmp_data[[paste0("L12.pval.",method)]] < tmp_enrich)}
      tmp_res$expected <- nrow(tmp_data)*tmp_enrich
      rbind(res, tmp_res) -> res; rm(tmp_res, tmp_data)
    }; rm(tmp_pop)
  }; rm(tmp_bin)
}; rm(tmp_enrich)
res$maf_bin <- paste(100*res$maf_low, "% to ", 100*res$maf_high, "%", sep="")
res$fc <- log2(res$observed/res$expected)

log2(sum(res$observed[res$enrichment == 0.01 & res$population == "L13"])/sum(res$expected[res$enrichment == 0.01 & res$population == "L13"]))
binom.test(sum(res$observed[res$enrichment == 0.01 & res$population == "L13"]), n = sum(res$expected[res$enrichment == 0.01 & res$population == "L13"])*100, p=0.01)$p.value

log2(sum(res$observed[res$enrichment == 0.01 & res$population == "D13"])/sum(res$expected[res$enrichment == 0.01 & res$population == "D13"]))
binom.test(sum(res$observed[res$enrichment == 0.01 & res$population == "D13"]), n = sum(res$expected[res$enrichment == 0.01 & res$population == "D13"])*100, p=0.01)$p.value


res -> res_real

######### plot panel 1A #########
p1A <- ggplot(res_real[res_real$population == "L13",], aes(x=1-enrichment , y=fc, color=maf_bin)) +
  geom_point(size=2) + geom_line() +
  theme_classic() + ylab("log2(enrichment)") +
  ggtitle("London: pre vs post") + coord_cartesian(xlim=c(0.8,1), ylim=c(-0.5,4.2)) +
  geom_abline(slope=0, intercept=0, color='dark gray') + theme(legend.position = "bottom") + 
  guides(fill = guide_legend(label.position = "bottom")) + xlab("percentile of neutral variants") + scale_color_manual(values=cbPalette)
rm(info, res, tmp_bins, method, info_neutral, info_candidate)
######## 

######### read in enrichment results for real data and for permuting sites #########
load("../DATA/perm_sites_results.RData")
######### 

######### calculate enrichment statistics for subsampled data #########
load("../DATA/pvalues_200neutralsites_subsampled.RData")
rm(calc_pval_D13, calc_pval_L12, calc_pval_L13, site, design, design2, input_fst, input_maf, drop_samples_missing, n_less, info_neut)
info <- info[!info$type == "neut",]
res_subsample <- res_real[1,]; res_subsample[,3:8] <- NA

percentile=0.01
res_subsample$fc_L13 <- sum(info$L13.pval.ML_same <= percentile)/(length(info$L13.pval.ML_same)*percentile)
res_subsample$fc_D13 <- sum(info$D13.pval.ML_same <= percentile)/(length(info$D13.pval.ML_same)*percentile)
res_subsample$p_L13 <- binom.test(sum(info$L13.pval.ML_same <= percentile), n = length(info$L13.pval.ML_same), p=percentile)$p.value
res_subsample$p_D13 <- binom.test(sum(info$D13.pval.ML_same <= percentile), n = length(info$D13.pval.ML_same), p=percentile)$p.value

res_subsample

rm(percentile, info)
######### lollipop plot for 1B #########
res_real$method <- "original"
res_subsample$method <- "subsample"
res_binom_format2 <- rbind(res_subsample, res_real)

tmp <- res_binom_format2[,c(1:2,9,3:4)]
colnames(tmp)[4:5] <- c("FC", "pval")
tmp$pop <- "London"
to_plot <- tmp
tmp <- res_binom_format2[,c(1:2,9,7:8)]
colnames(tmp)[4:5] <- c("FC", "pval")
tmp$pop <- "Denmark"
to_plot <- rbind(to_plot,tmp); rm(tmp, res_binom_format2)

p1B <- ggplot(to_plot[to_plot$enr == 0.01,], aes(y=FC, x=factor(method), col=method)) + 
  facet_grid( pop ~ . , switch = "y") +
  geom_segment(aes(xend=factor(method), y=0,yend=FC,),
               color = "gray", lwd = 1.5) + 
  scale_color_manual(values = c("red", "orange")) + 
  geom_point(size = 8) + ylim(c(0,5.5)) +
  coord_flip() + 
  xlab("") + 
  ylab("degree of enrichment") +
  geom_text(aes(label = "*", x=factor(method)), color = "black", size = 6) +
  theme_classic() + theme(legend.position="none", axis.text.y = element_blank())
p1B
rm(to_plot)
######### 


######### plot 1C #########
p1C <- ggplot(data=res_perm, aes(x=fc_D13, y=fc_L13)) + 
  stat_density_2d(aes(fill = ..level..), geom = "polygon") + 
  ggtitle("ML, 200 sites: 99th percentile") + 
  xlab("Denmark enrichment") + ylab("London enrichment") + 
  scale_fill_distiller(palette=4, direction=1) +
  coord_cartesian(xlim=c(0,5), ylim=c(0,5.5)) +
  geom_point(x = res_real$fc_D13, y=res_real$fc_L13, col="red", size=7) +
  geom_point(x = res_subsample$fc_D13 , y=res_subsample$fc_L13, col="orange", size=7) +
  theme_bw() + theme(legend.position = 'none'); p1C

#########



(p1A + p1B + p1C) + 
  plot_layout(guides="collect", width=c(1,1,1.5)) + plot_annotation(tag_levels = 'A') & 
  theme(legend.position = 'bottom')


