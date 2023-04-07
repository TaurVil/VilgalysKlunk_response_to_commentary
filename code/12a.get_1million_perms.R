library(parallel); library(data.table)

# Load in the right data set which has our variant to permute
load("./setup_for_permutations.exon.RData")
if (! "MYSITE" %in% info$site) {load("./setup_for_permutations.gwas.RData")}
if (!"MYSITE" %in% info$site) {load("./setup_for_permutations.neut.RData")}

info <- info[info$site %in% "MYSITE",]
perm_individuals <- function(perm) {
  # Sample individuals with replacement for each population and time point
  D_snames_1 <- sample(D_snames_joined, length(D_snames_1), replace = T)
  D_snames_3 <- sample(D_snames_joined, length(D_snames_3), replace = T)
  L_snames_1 <- sample(L_snames_joined, length(L_snames_1), replace = T)
  L_snames_2 <- sample(L_snames_joined, length(L_snames_2), replace = T)
  L_snames_3 <- sample(L_snames_joined, length(L_snames_3), replace = T)

  info <- get("info")

  for (pop in c("london", "denmark")) {
    if (pop == "london") {pop2 <- "L"}
    if (pop == "denmark") {pop2 <- "D"}

    if (pop == "london") {times <- c("pre", "post", "during")}
    if (pop == "denmark") {times <- c("pre", "post")}

    NAMES <- get(paste0("NAMES.",pop2))
    DATA <- get(paste0("DATA.",pop2))
    DATA <- DATA[row.names(DATA) %in% info$site,]

    for (time in times) {
      if (time == "pre") {t2 <- "1"}
      if (time == "during") {t2 <- "2"}
      if (time == "post") {t2 <- "3"}

      keep <- get(paste0(pop2, "_snames_", t2))
      keep <- which(NAMES %in% keep)
      k2 <- as.data.frame(matrix(nrow=nrow(DATA), ncol=1))
      for (k in keep) {k2 <- cbind(k2, DATA[,(k*3-2):(k*3)])};
      tmp_data <- k2[,-1]; rm(k, keep, k2)
      n <- ncol(tmp_data)/3

      info[[paste(pop, time, "called", sep=".")]] <- (ncol(tmp_data)-rowSums(is.na(tmp_data)))/3
      info[[paste(pop,time,"alternate",sep=".")]] <- estimate_af_ml(tmp_data)
    }; rm(time)
  }
  info$perm <- perm
  return(info)
}

# Get allele frequency from 10k permutations
X=1:10000
perm_sites <- do.call("rbind", mclapply(X, perm_individuals))


######## get fst and allele frequency ########
perm_sites$london.alternate <- rowMeans(cbind(perm_sites$london.pre.alternate, perm_sites$london.post.alternate, perm_sites$london.during.alternate))
perm_sites$denmark.alternate <- rowMeans(cbind(perm_sites$denmark.pre.alternate, perm_sites$denmark.post.alternate))
perm_sites$alternate <- rowMeans(cbind(perm_sites$london.alternate, perm_sites$denmark.alternate))
perm_sites$maf <- apply(cbind(perm_sites$alternate, 1-perm_sites$alternate), 1, min)

for (pop in c("london", "denmark")) {
  if (pop == "london") {times <- c("pre", "post", "during")}
  if (pop == "denmark") {times <- c("pre", "post")}

  for (time in times) {
    # expected heterozygosity
    perm_sites[[paste(pop,time,"ehet",sep=".")]] <- 2*perm_sites[[paste(pop,time,"alternate",sep=".")]]*(1-perm_sites[[paste(pop,time,"alternate",sep=".")]])
  }
  rm(time, times)
}; rm(pop)
# pairwise means, pairwise expectations, and Fst
for (i in 1:3) {
  n2=paste(design2$pop[i],design2$time1[i],design2$time2[i],sep=".")
  perm_sites[[paste(n2,"mean",sep=".")]] <- rowMeans(cbind(perm_sites[[paste(design2$pop[i],".",design2$time1[i],".alternate",sep="")]], perm_sites[[paste(design2$pop[i],".",design2$time2[i],".alternate",sep="")]]))
  perm_sites[[paste(n2,"ehet",sep=".")]] <- 2*perm_sites[[paste(n2,"mean",sep=".")]]*(1-perm_sites[[paste(n2,"mean",sep=".")]])
  perm_sites[[paste(n2,"fst",sep=".")]] <- (perm_sites[[paste(n2,"ehet",sep=".")]] - (perm_sites[[paste(design2$pop[i],".",design2$time1[i],".ehet",sep="")]] + perm_sites[[paste(design2$pop[i],".",design2$time2[i],".ehet",sep="")]])/2)/perm_sites[[paste(n2,"ehet",sep=".")]]
}; rm(i,n2)
########
perm_sites$delta_L13 <- perm_sites$london.post.alternate - perm_sites$london.pre.alternate
perm_sites$delta_L12 <- perm_sites$london.during.alternate - perm_sites$london.pre.alternate
perm_sites$delta_D13 <- perm_sites$denmark.post.alternate - perm_sites$denmark.pre.alternate
########

# If the file already exists, write to it without column headers. If the output file doesn't exist (ie this is the first set of 10k sites), write the file header as well. 
if (file.exists("./perms.MYSITE.txt")) {write.table(perm_sites, "./perms.MYSITE.txt", row.names=F, col.names=F, sep="\t", quote=F, append=T)} else {write.table(perm_sites, "./perms.MYSITE.txt", row.names=F, col.names=T, sep="\t", quote=F)}

