# We call haploid data using angsd

# ls rescaled_bams/*bam > haploid.bamlist
# ../../Programs/angsd/angsd/angsd -bam haploid.bamlist -dohaplocall 1 -doCounts 1 -minMinor 1 -out iteration1 -seed 1
# # -dohaplocall 1 : sample a random base 
# # -r 1:  would say use just chromosome 1


# library(data.table)
# onek <- fread('./DATA/sites_for_genotyping.txt')
# haplo <- fread("./angsdput.haplo.gz"); noms <- fread("./haploid.bamlist", header=F)
# 
# noms$V1 <- gsub("rescaled_bams/merged.", "", noms$V1)
# noms$V1 <- gsub(".rescaled.bam", "", noms$V1)
# colnames(haplo)[4:ncol(haplo)] <- noms$V1 
# 
# onek$site <- paste(onek$V1, onek$V2, sep="_")
# haplo$site <- paste(haplo$chr, haplo$pos, sep="_")
# h2 <- haplo[haplo$site %in% onek$site,]
# rm(haplo, onek, noms)
# save.image("../DATA/pseudohaploid.1.RData")





library(data.table); library(ggplot2); library(patchwork)
load("../DATA/fst_estimates.RData")
n_nearby <- 200; min_maf <- 0.05

######## Initiation, parameters, and functions ######## 
# parameters
min_n=10
min_maf=0.0
drop_samples_missing=.5
########

######## get sample names per population and time point ########
L_snames_1 <- unique(c(
  read.delim(paste("../DATA/SampleNames/keep_london_pre_exon",sep=""), header=F)$V1, 
  read.delim(paste("../DATA/SampleNames/keep_london_pre_gwas",sep=""), header=F)$V1, 
  read.delim(paste("../DATA/SampleNames/keep_london_pre_neutral",sep=""), header=F)$V1
))
L_snames_2 <- unique(c(
  read.delim(paste("../DATA/SampleNames/keep_london_during_exon",sep=""), header=F)$V1, 
  read.delim(paste("../DATA/SampleNames/keep_london_during_gwas",sep=""), header=F)$V1, 
  read.delim(paste("../DATA/SampleNames/keep_london_during_neutral",sep=""), header=F)$V1
))
L_snames_3 <- unique(c(
  read.delim(paste("../DATA/SampleNames/keep_london_post_exon",sep=""), header=F)$V1, 
  read.delim(paste("../DATA/SampleNames/keep_london_post_gwas",sep=""), header=F)$V1, 
  read.delim(paste("../DATA/SampleNames/keep_london_post_neutral",sep=""), header=F)$V1
))
L_snames_joined <- c(L_snames_1, L_snames_2, L_snames_3)
# same for denmark
D_snames_1 <- unique(c(
  read.delim(paste("../DATA/SampleNames/keep_denmark_pre_exon",sep=""), header=F)$V1, 
  read.delim(paste("../DATA/SampleNames/keep_denmark_pre_gwas",sep=""), header=F)$V1, 
  read.delim(paste("../DATA/SampleNames/keep_denmark_pre_neutral",sep=""), header=F)$V1
))
D_snames_3 <- unique(c(
  read.delim(paste("../DATA/SampleNames/keep_denmark_post_exon",sep=""), header=F)$V1, 
  read.delim(paste("../DATA/SampleNames/keep_denmark_post_gwas",sep=""), header=F)$V1, 
  read.delim(paste("../DATA/SampleNames/keep_denmark_post_neutral",sep=""), header=F)$V1
))
D_snames_joined <- c(D_snames_1, D_snames_3)
######## PERMUTATION STEP -- EMPTY FOR THIS SCRIPT WITH THE OBSERVED DATA ########
# Empty for real data
# For permutations, see other examples in this directory
# Convert new list of names to L_snames_[1:3] and D_snames_[1:2]
rm(L_snames_joined, D_snames_joined)
########

load("../DATA/pseudohaploid.1.RData"); # rm(haplo, onek, noms); save.image("../DATA/pseudohaploid.1.RData")
h2 <- as.data.frame(h2[h2$site %in% info$site,])
info <- as.data.frame(info[info$site %in% h2$site, c(1:5,21,65,29)])
h2 <- h2[order(h2$site),]
info <- info[order(info$site),]
sum(h2$site == info$site)

for (pop in c("london", "denmark")) {
  if (pop == "london") {pop2 <- "L"}; if (pop == "denmark") {pop2 <- "D"}
  if (pop == "london") {times <- c("1", "2", "3")}; if (pop == "denmark") {times <- c("1", "3")}
  
  for (time in times) { 
    t_names <- get(paste0(pop2,"_snames_",time)) 
    t_data <- h2[,which(colnames(h2) %in% t_names)]
    row.names(t_data) <- h2$site
    
    t_data[t_data == "N"] <- NA
    info[[paste0(pop2,time,".called")]] <- ncol(t_data) - rowSums(is.na(t_data))
    info[[paste0(pop2,time,".alternate")]] <- NA
    for (i in 1:nrow(info)) {
      t_ref <- info$ref[i]
      t_alt <- info$alt[i]
      n_ref <- sum(t_data[i,] == t_ref, na.rm=T)
      n_alt <- sum(t_data[i,] == t_alt, na.rm=T)
      info[[paste0(pop2,time,".alternate")]][i] <- n_alt/(n_ref + n_alt)
      rm(n_ref, n_alt, t_ref, t_alt)
    }; rm(i)
    rm(t_data, t_names)
  }
  rm(time, times)
}; rm(pop, pop2)


####### filter for sites with 10 samples per population in each time point #######
attach(info)
info <- subset(info, apply(X = cbind(L1.called, L2.called, L3.called, D1.called, D3.called), 1, min) >= min_n)
detach(info)
#######

####### calculate alternate and minor allele frequency #######
# do for both ML and GL calls, even though we only use ML for this response
info$london.alternate <- rowMeans(cbind(info$L1.alternate, info$L2.alternate, info$L3.alternate))
info$denmark.alternate <- rowMeans(cbind(info$D1.alternate, info$D3.alternate))
info$alternate <- rowMeans(cbind(info$london.alternate, info$denmark.alternate))
info$maf <- apply(cbind(info$alternate, 1-info$alternate), 1, min)

plot(info$maf ~ info$maf.ML)
#######

####### calculate fst #######
# calculate expected heterozygosity
for (pop in c("L", "D")) {
  if (pop == "L") {times <- c("1", "2", "3")}
  if (pop == "D") {times <- c("1", "3")}
  
  for (time in times) {
    # expected heterozygosity
    info[[paste(pop,time,"ehet",sep=".")]] <- 2*info[[paste(pop,time,".alternate",sep="")]]*(1-info[[paste(pop,time,".alternate",sep="")]])
  }
  rm(time, times)
}; rm(pop)
# pairwise means, pairwise expectations, and Fst
design2$pop <- gsub("denmark", "D", design2$pop)
design2$pop <- gsub("london", "L", design2$pop)
design2$time2 <- gsub("pre", "1", design2$time2)
design2$time1 <- gsub("post", "3", design2$time1)
design2$time1 <- gsub("during", "2", design2$time1)

for (i in 1:3) {
  n2=paste(design2$pop[i],design2$time2[i],design2$time1[i],sep="")
  info[[paste(n2,"mean",sep=".")]] <- rowMeans(cbind(info[[paste(design2$pop[i],design2$time1[i],".alternate",sep="")]], info[[paste(design2$pop[i],design2$time2[i],".alternate",sep="")]]))
  info[[paste(n2,"ehet",sep=".")]] <- 2*info[[paste(n2,"mean",sep=".")]]*(1-info[[paste(n2,"mean",sep=".")]])
  info[[paste(n2,"fst",sep=".")]] <- (info[[paste(n2,"ehet",sep=".")]] - (info[[paste(design2$pop[i],".",design2$time1[i],".ehet",sep="")]] + info[[paste(design2$pop[i],".",design2$time2[i],".ehet",sep="")]])/2)/info[[paste(n2,"ehet",sep=".")]]

}; rm(i,n2)
########

######## Get difference between time points ########
info$delta_L13 <- info$L3.alternate - info$L1.alternate
info$delta_L12 <- info$L2.alternate - info$L1.alternate
info$delta_D13 <- info$D3.alternate - info$D1.alternate

########


save.image("../DATA/fst_estimates_pseudohaploid.RData")

