library(data.table); library(ggplot2); library(patchwork)
load("../DATA/fst_estimates.RData")
n_nearby <- 200; min_maf <- 0.05

######## 
# Load in ML data - original
########
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


for (pop in c("london", "denmark")) {
  if (pop == "london") {pop2 <- "L"}
  if (pop == "denmark") {pop2 <- "D"}
  
  if (pop == "london") {times <- c("pre", "post", "during")}
  if (pop == "denmark") {times <- c("pre", "post")}
  
  for (time in times) {
    # get genolik file
    d <- as.data.frame(fread(paste("DATA/subsample_results/GL/genolik.",pop,"_", time, ".genolik",sep="")))
    # remove site. change missing data to NA. get number of samples.
    row.names(d) <- paste(d[,1], d[,2], sep="_")
    d <- d[,-(1:2)] ; d[d == -1] <- NA; tmp_n <- ncol(d)/3
    # report data on missingness per sample
    tmp_missing=(colSums(is.na(d))/nrow(d))[seq(1,3*tmp_n,3)]
    # remove samples missing too much data
    keep <- which(tmp_missing <= drop_samples_missing)
    k2 <- as.data.frame(matrix(nrow=nrow(d), ncol=1)); for (k in keep) {k2 <- cbind(k2, d[,(k*3-2):(k*3)])}
    d <- k2[,-1];
    n <- ncol(d)/3
    tmp_info <- info[info$site %in% row.names(d),1]
    d <- d[order(row.names(d)),]
    tmp_info$site <- row.names(d)
    # tmp_info[[paste(pop,time,"alternate.GL_same",sep=".")]] <- rowMeans(d[,seq(3,n*3,3)],na.rm=T) + rowMeans(d[,seq(2,n*3,3)],na.rm=T)/2
    tmp_info[[paste(pop,time,"alternate.ML_same",sep=".")]] <- estimate_af_ml(d)
    tmp_info[[paste(pop,time,"n_called.ML_same",sep=".")]] <- rowSums(!is.na(d))/3
    tmp_info <- tmp_info[tmp_info[[paste(pop,time,"n_called.ML_same",sep=".")]] >= 10,]
    
    rm(d, n, k2,keep, tmp_missing,k)
    info <- merge(info, tmp_info, by="site", all=T); rm(tmp_info)
  }
}; rm(pop, pop2, times, time)

########

#######
# calculate alternate and minor allele frequency
#######
info$london.alternate.ML_same <- rowMeans(cbind(info$london.pre.alternate.ML_same, info$london.post.alternate.ML_same,info$london.during.alternate.ML_same))
info$denmark.alternate.ML_same <- rowMeans(cbind(info$denmark.pre.alternate.ML_same, info$denmark.post.alternate.ML_same))
info$alternate.ML_same <- rowMeans(cbind(info$london.alternate.ML_same, info$denmark.alternate.ML_same))
info$maf.ML_same <- apply(cbind(info$alternate.ML_same, 1-info$alternate.ML_same), 1, min)

#######

#######
# calculate fst
#######
for (pop in c("london", "denmark")) {
  if (pop == "london") {times <- c("pre", "post", "during")}
  if (pop == "denmark") {times <- c("pre", "post")}
  
  for (time in times) {
    # expected heterozygosity
    info[[paste(pop,time,"ehet.ML_same",sep=".")]] <- 2*info[[paste(pop,time,"alternate.ML_same",sep=".")]]*(1-info[[paste(pop,time,"alternate.ML_same",sep=".")]])
  }
  rm(time, times)
}; rm(pop)
# pairwise means, pairwise expectations, and Fst
for (i in 1:3) {
  n2=paste(design2$pop[i],design2$time1[i],design2$time2[i],sep=".")
  
  info[[paste(n2,"mean.ML_same",sep=".")]] <- rowMeans(cbind(info[[paste(design2$pop[i],".",design2$time1[i],".alternate.ML_same",sep="")]], info[[paste(design2$pop[i],".",design2$time2[i],".alternate.ML_same",sep="")]]))
  info[[paste(n2,"ehet.ML_same",sep=".")]] <- 2*info[[paste(n2,"mean.ML_same",sep=".")]]*(1-info[[paste(n2,"mean.ML_same",sep=".")]])
  info[[paste(n2,"fst.ML_same",sep=".")]] <- (info[[paste(n2,"ehet.ML_same",sep=".")]] - (info[[paste(design2$pop[i],".",design2$time1[i],".ehet.ML_same",sep="")]] + info[[paste(design2$pop[i],".",design2$time2[i],".ehet.ML_same",sep="")]])/2)/info[[paste(n2,"ehet.ML_same",sep=".")]]
  
}; rm(i,n2)
########

info$delta_L13.ML_same <- info$london.post.alternate.ML_same-info$london.pre.alternate.ML_same
info$delta_L12.ML_same <- info$london.during.alternate.ML_same-info$london.pre.alternate.ML_same
info$delta_D13.ML_same <- info$denmark.post.alternate.ML_same-info$denmark.pre.alternate.ML_same

rm(n_nearby, tmp_n)
rm(estimate_af_ml, likelihood)
save.image("./DATA/fst_estimates_subsampled.RData")
