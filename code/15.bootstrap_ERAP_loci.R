library(data.table); library(ggplot2)

load("./20Feb2022_all3_estimates_no_pvalues.RData")
rm(method, min_maf, min_n, min_neutral_nsites, tmp_n, window_size, drop_samples_missing, mafs, bins, design2)

sites <- c("chr5_96244549", "chr5_96244585")

likelihood <- function(p, data){
  gt.freq <- c((1-p)^2, 2*p*(1-p), p^2)
  ll <- sum(log(rowSums(t(t(data)*gt.freq))))
  return(ll)
}
estimate_af_ml<-function(d){
  n<-NROW(d)
  af<-rep(NA, n)
  for(i in 1:n){
    gl<-d[i,,drop=FALSE]
    gl <- gl[!is.na(gl)]
    gl <- matrix(gl, ncol=3, byrow=T)
    gl <- gl[!is.na(gl[,1]),,drop=FALSE]
    opt <- optimize(likelihood, interval=c(0,1), maximum=TRUE, data=gl)
    af[i] <- opt$maximum
  }
  return(af)
}

info <- info[info$site %in% sites,]

plot_info <- NULL
for (ttt in c(4:nrow(design),1:3)) {
  name=paste("./DATA/genoliks/genolik.gwas_",design$pop[ttt],"_", design$time[ttt],".genolik",sep="")
  n2=paste(design$pop[ttt],design$time[ttt],sep="_")
  d <- as.data.frame(fread(name))
  row.names(d) <- paste(d$V1, d$V2, sep="_")
  d <- d[,-c(1:2)]
  d <- d[row.names(d) %in% sites,]
  d[d == -1] <- NA; n <- ncol(d)/3
  
  #### bootstrap step 
  boot_data <- NULL
  for (i in 1:10000) {
    k <- sample(1:n, n, replace=T)
    new_d <- as.data.frame(matrix(ncol=1, nrow=nrow(d)))
    for (j in k) {
      new_d <- cbind(new_d, d[,c((j*3 - 2):(j*3))])
    }; rm(j)
    new_d[,1] <- row.names(d)
    colnames(new_d) <- c("site", paste0("V",1:(3*n)))
    boot_data <- rbind(boot_data, new_d); rm(new_d, k)
  }; rm(i)
  rm(d)
  
  boot_info <- as.data.frame(boot_data$site)
  colnames(boot_info)[1] <- "site"
  boot_info$freq <- estimate_af_ml(boot_data[,-1])
  
  tmp_info <- info[,c(1:5)]
  tmp_info$pop <- design$pop[ttt]
  tmp_info$time <- design$time[ttt]
  tmp_info$real_af <- info[[paste0(design$pop[ttt],".",design$time[ttt],".alternate.ML")]]
  
  for (s in sites) {
    mean(boot_info$freq[boot_info$site == s]) -> tmp_info$boot.mean[tmp_info$site == s]
    quantile(boot_info$freq[boot_info$site == s], 0.05) -> tmp_info$boot.lower95[tmp_info$site == s]
    quantile(boot_info$freq[boot_info$site == s], 0.95) -> tmp_info$boot.upper95[tmp_info$site == s]
    sd(boot_info$freq[boot_info$site == s]) -> tmp_info$boot.sd[tmp_info$site == s]
  }; rm(s)
  print(design[ttt,])
  
  rbind(plot_info, tmp_info) -> plot_info
  
  rm(n, n2, name, tmp_info, boot_data, boot_info)
  
}

rm(ttt)

plot_info$time <- factor(plot_info$time, levels=c("pre", "during", "post"))
plot_info$xpos <- as.numeric(plot_info$time)
plot_info$xpos[plot_info$pop == 'denmark'] <- plot_info$xpos[plot_info$pop == 'denmark'] + 0.1

plot_info[,c(8:11)] <- 1-plot_info[,c(8:11)]
plot(data=plot_info, boot.mean ~ real_af); abline(a=0,b=1,col='red')

plot1 <- ggplot(data=plot_info, aes(x=xpos, y=boot.mean, col=pop)) +
  facet_grid( site ~ .) + 
  geom_point(size=4) + 
  geom_line(linetype="dashed") + 
  coord_cartesian(ylim=c(0.2,1)) +
  # geom_segment(aes(x=xpos, xend=xpos, y=boot.lower95, yend=boot.upper95), size=1.1) +
  # geom_segment(aes(x=xpos, xend=xpos, y=real_af - boot.sd, yend=real_af + boot.sd), size=1.1) +
  geom_segment(aes(x=xpos, xend=xpos, y=boot.mean - boot.sd, yend=boot.mean + boot.sd), size=1.1) +
  # ggtitle("linked ERAP2 variants") +
  scale_color_manual(values=c("#56B4E9", 'red', "#56B4E9", 'red')) + 
  theme_classic() + theme(legend.position = 'inset') + xlab("time") + ylab("derived allele frequency") 
plot1

plot1 -> plot_erap_loci; rm(plot1)
rm(estimate_af_ml, likelihood, design)
info -> erap_info
plot_info -> erap_plot_info; rm(plot_info, info)
sites -> erap_sites; rm(sites)
save.image("./bootstrap_ERAP_loci.RData")
