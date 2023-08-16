library(ggplot2); library(data.table); library(patchwork)
################## Read in data, set offset ########################
load("../DATA/pvalues.RData")

offset=9e-7; method="ML"
res <- res[apply(abs(res[,which(colnames(res) %like% "delta" & colnames(res) %like% method)]),1,min) > 0,]

table(res$type)

original_sites <- c("chr5_96244549",
                    "chr2_204738938",
                    "chr5_114915460", 
                    "chr18_77287776")

# calculate p-value from L & D independently
xi <- -2*(log(res$L_1v3.pval.ML.conditioned) + log(res$D_1v3.pval.ML.conditioned))
res$p.conditioned <- 1 - pchisq(xi, 4)


####################

res2 <- res

res2$logp <- -log10(res2$p.conditioned + offset)
res2$rank <- rank(res2$logp)/length(res2$logp)

res2$sh <- NA; 
res2$sh[sign(res2$delta_L13.ML) != sign(res2$delta_D13.ML)] <- 1
res2$sh[sign(res2$delta_L13.ML) == sign(res2$delta_L12.ML) & 
          sign(res2$delta_L13.ML) == sign(res2$delta_D13.ML)] <- 21
res2$sh[sign(res2$delta_L13.ML) != sign(res2$delta_L12.ML) & 
          sign(res2$delta_L13.ML) == sign(res2$delta_D13.ML)] <- 16
res2$sh <- factor(res2$sh)

res2$str <- NA; 
res2$str[sign(res2$delta_L13.ML) != sign(res2$delta_D13.ML)] <- .1
res2$str[sign(res2$delta_L13.ML) == sign(res2$delta_L12.ML) & 
           sign(res2$delta_L13.ML) == sign(res2$delta_D13.ML)] <- .1
res2$str[sign(res2$delta_L13.ML) != sign(res2$delta_L12.ML) & 
           sign(res2$delta_L13.ML) == sign(res2$delta_D13.ML)] <- 1
res2$str[1] <- 0
# res2$fill[1] <- 0
library(patchwork); library(data.table)

orig <- res2[res2$site %in% original_sites,]
res2 <- res2[! res2$site %in% original_sites,]

# shape=sh, 
cbPalette <- c("darkolivegreen", "darkolivegreen", "darkorange3", "darkorange3", "darkolivegreen", "#0072B2", "#D55E00", "#CC79A7")
fig2a <- ggplot(data=res2[res2$type != "neut",], aes(x=rank, y=logp, alpha=str, color=factor(type))) +
  # facet_grid( . ~ sh) +
  geom_jitter(width = 0.04, size=3, alpha=0.8) + 
  ylab("-log10( p-value )") + 
  xlab("ranked evidence") +
  ylim(c(0,3)) +
  geom_point(x = orig$rank[orig$site == "chr5_96244549"],
             y = orig$logp[orig$site == "chr5_96244549"],
             col="firebrick1", size=8) +
  # geom_point(x = res$rank[res$site == "chr2_204738938"],
  #            y = res$logp[res$site == "chr2_204738938"],
  #            col="#D55E00", size=4, shape=16) +
  # geom_point(x = res$rank[res$site == "chr5_114915460"],
  #            y = res$logp[res$site == "chr5_114915460"],
  #            col="#D55E00", size=4, shape=16) +
  geom_point(x = orig$rank[orig$site == "chr18_77287776"],
             y = orig$logp[orig$site == "chr18_77287776"],
             col="firebrick1", size=8) +
  geom_hline(yintercept = -log10(0.05), col='orange', lty=2) +
  geom_hline(yintercept = -log10(0.01), col='red', lty=2) +
  # geom_hline(yintercept = -log10(0.05/1723), col='black', lty=2) +
  theme_classic() + theme(legend.position="none") +
  scale_color_manual(values=cbPalette); fig2a


saveRDS(fig2a, "./fig2a_revision.RData")
