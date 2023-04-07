sites=read.delim("./DATA/sites_maf_to_permute.txt", header=F)

min_n=10; min_maf=0.025
drop_samples_missing=.5
######## get list of all candidate snames per population ########
# get list of all candidate snames per population
L_snames_1 <- unique(c(
  read.delim(paste("./DATA/SampleNames/keep_london_pre_exon",sep=""), header=F)$V1, 
  read.delim(paste("./DATA/SampleNames/keep_london_pre_gwas",sep=""), header=F)$V1, 
  read.delim(paste("./DATA/SampleNames/keep_london_pre_neutral",sep=""), header=F)$V1
))
L_snames_2 <- unique(c(
  read.delim(paste("./DATA/SampleNames/keep_london_during_exon",sep=""), header=F)$V1, 
  read.delim(paste("./DATA/SampleNames/keep_london_during_gwas",sep=""), header=F)$V1, 
  read.delim(paste("./DATA/SampleNames/keep_london_during_neutral",sep=""), header=F)$V1
))
L_snames_3 <- unique(c(
  read.delim(paste("./DATA/SampleNames/keep_london_post_exon",sep=""), header=F)$V1, 
  read.delim(paste("./DATA/SampleNames/keep_london_post_gwas",sep=""), header=F)$V1, 
  read.delim(paste("./DATA/SampleNames/keep_london_post_neutral",sep=""), header=F)$V1
))
L_snames_joined <- c(L_snames_1, L_snames_2, L_snames_3)
# same for denmark
D_snames_1 <- unique(c(
  read.delim(paste("./DATA/SampleNames/keep_denmark_pre_exon",sep=""), header=F)$V1, 
  read.delim(paste("./DATA/SampleNames/keep_denmark_pre_gwas",sep=""), header=F)$V1, 
  read.delim(paste("./DATA/SampleNames/keep_denmark_pre_neutral",sep=""), header=F)$V1
))
D_snames_3 <- unique(c(
  read.delim(paste("./DATA/SampleNames/keep_denmark_post_exon",sep=""), header=F)$V1, 
  read.delim(paste("./DATA/SampleNames/keep_denmark_post_gwas",sep=""), header=F)$V1, 
  read.delim(paste("./DATA/SampleNames/keep_denmark_post_neutral",sep=""), header=F)$V1
))
D_snames_joined <- c(D_snames_1, D_snames_3)
########
######## Get set of sites to consider ########
# Create design matrix for which we'll pull in data
design <- expand.grid(time = gl(3, 1, labels = c("pre", "post", "during")), pop = gl(2, 1, labels = c("denmark", "london")))
design <- design[-3,]
#type = gl(3, 1, labels = c("exons", "neut", "neutral"))
design2 <- expand.grid(time1 = gl(3, 1, labels = c("pre", "post", "during")), time2 = gl(3, 1, labels = c("pre", "post", "during")), pop = gl(2, 1, labels = c("denmark", "london")))
design2<-subset(design2, !(design2$time1 == design2$time2)); design2 <- design2[c(1,7:8,10),]
######## Read in info file ########
info <- fread("./DATA/genoliks/neutral.frq")[,-c(3:4)]
colnames(info) <- c("chr", "pos", "ref", "alt")
info$ref <- gsub(":.*", "",info$ref); info$alt <- gsub(":.*", "",info$alt)
paste(info$chr,info$pos,sep="_") -> info$site
########
######## Functions for maximum likelihood estimates of allele frequency ######## 
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

#######
####### set type ###########
type="neut"
if (type == "gwas") {type2 <- "gwas"}
if (type == "exon") {type2 <- "exons"}
if (type == "neut") {type2 <- "neutral"}
####### read in data for those sites, then permute ###########
for (pop in c("london", "denmark")) {
  if (pop == "london") {pop2 <- "L"}
  if (pop == "denmark") {pop2 <- "D"}
  
  if (pop == "london") {times <- c("pre", "post", "during")}
  if (pop == "denmark") {times <- c("pre", "post")}
  
  # p - by - 3*n matrix of genotype likelihoods, and the corresponding list of sample names 
  NAMES <- NULL
  DATA <- as.data.frame(matrix(nrow=nrow(get("info")), ncol=1))
  for (time in times) {
    
    # read in sample names 
    tmp_names=read.delim(paste("./DATA/SampleNames/", pop, "_", time, "_",type2,".012.indv",sep=""), header=F)
    # get genolik file 
    d <- as.data.frame(fread(paste("./DATA/genoliks/genolik.",type2,"_",pop,"_", time, ".genolik",sep="")))
    # remove site. change missing data to NA. get number of samples.
    row.names(d) <- paste(d[,1], d[,2], sep="_")
    d <- d[,-(1:2)] ; d[d == -1] <- NA; tmp_n <- ncol(d)/3
    # report data on missingness per sample
    tmp_missing=(colSums(is.na(d))/nrow(d))[seq(1,3*tmp_n,3)]
    # remove samples missing too much data
    keep <- which(tmp_missing <= drop_samples_missing)
    k2 <- as.data.frame(matrix(nrow=nrow(d), ncol=1)); for (k in keep) {k2 <- cbind(k2, d[,(k*3-2):(k*3)])}
    d <- k2[,-1];
    tmp_names <- tmp_names$V1[keep]; n_pre <- ncol(d)/3; rm(n_pre)
    
    # SAVE TO FILE FOR pop and type 
    NAMES <- c(NAMES, tmp_names)
    DATA <- cbind(DATA, d)
    rm(k2, k, tmp_missing,keep)
  }; rm(time, times, d, tmp_n, tmp_names)
  # Remove buffer row from data
  DATA <- DATA[,-1]
  
  DATA <- DATA[row.names(DATA) %in% sites$V1,]
  
  assign(value = DATA, x=paste0("DATA.",pop2))
  assign(value = NAMES, x=paste0("NAMES.",pop2))
  rm(NAMES,DATA)
}; rm(pop, pop2, type2)

info <- info[info$site %in% sites$V1,]

info <- info[!duplicated(info$site),]

save.image("./setup_for_permutations.neut.RData")


