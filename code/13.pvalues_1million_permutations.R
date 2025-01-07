library(data.table); library(parallel)

########## Read in frequencies without permutation ########## 
load("./fst_estimates.RData")
info <- info[(info$rsid_number == 1),]
info <- info[!duplicated(info$site),] 
info <- info[info$maf.ML >= 0.05,]
##########  Rename info ########## 
info -> info_real; rm(info, design, design2, min_maf, min_n, drop_samples_missing)
info_real <- as.data.frame(info_real)

info_real <- info_real[sign(info_real$delta_L13.ML) == sign(info_real$delta_D13.ML) & sign(info_real$delta_L13.ML) != sign(info_real$delta_L12.ML),]


method="ML"

get_ps <- function(tmp_n) {
  tmp_snp <- info_real$site[tmp_n]
  tmp_file <- paste0("./perms.no_replacement.",tmp_snp,".txt.gz")
  
  res <- info_real[info_real$site == tmp_snp,]
  
  if (file.exists(tmp_file)) {
    extra_perms <- fread(tmp_file)
    
    res$number.neutral.all <- nrow(extra_perms)
    
    res$D13.pval.ML.all <- sum(
      extra_perms$site == tmp_snp & 
        extra_perms$denmark.post.pre.fst >= res$denmark.post.pre.fst.ML
    )/sum(extra_perms$site == tmp_snp)
    
    res$L13.pval.ML.all <- sum(
      extra_perms$site == tmp_snp & 
        extra_perms$london.post.pre.fst >= res$london.post.pre.fst.ML
    )/sum(extra_perms$site == tmp_snp)
    
    extra_perms <- extra_perms[sign(extra_perms$delta_L13) == sign(extra_perms$delta_D13) & 
                                 sign(extra_perms$delta_L13) != sign(extra_perms$delta_L12),]
    
    res$number.neutral.dx <- nrow(extra_perms)
    
    res$D13.pval.ML.dx <- sum(
      extra_perms$site == tmp_snp & 
        extra_perms$denmark.post.pre.fst >= res$denmark.post.pre.fst.ML
    )/sum(extra_perms$site == tmp_snp)
    
    res$L13.pval.ML.dx <- sum(
      extra_perms$site == tmp_snp & 
        extra_perms$london.post.pre.fst >= res$london.post.pre.fst.ML
    )/sum(extra_perms$site == tmp_snp)
    
    # 
    # res$double.pval.ML <- sum(
    #   extra_perms$site == tmp_snp & 
    #     extra_perms$london.post.pre.fst >= res$london.post.pre.fst.ML &
    #     extra_perms$denmark.post.pre.fst >= res$denmark.post.pre.fst.ML 
    # )/sum(extra_perms$site == tmp_snp)
    # 
    # res$double_plus_sign.pval.ML <- sum(
    #   extra_perms$site == tmp_snp & 
    #     extra_perms$london.post.pre.fst >= res$london.post.pre.fst.ML &
    #     extra_perms$denmark.post.pre.fst >= res$denmark.post.pre.fst.ML &
    #     sign(extra_perms$delta_L13) == sign(extra_perms$delta_D13)
    # )/sum(extra_perms$site == tmp_snp)
    # 
    # res$double_sign_during.pval.ML <- sum(
    #   extra_perms$site == tmp_snp & 
    #     extra_perms$london.post.pre.fst >= res$london.post.pre.fst.ML &
    #     extra_perms$denmark.post.pre.fst >= res$denmark.post.pre.fst.ML &
    #     sign(extra_perms$delta_L13) == sign(extra_perms$delta_D13) & 
    #     sign(extra_perms$delta_L13) != sign(extra_perms$delta_L12)
    # )/sum(extra_perms$site == tmp_snp)
    
    return(res)
  }
}

res <- NULL
for (i in 1:nrow(info_real)) {
  d <- get_ps(i)
  res <- rbind(res,d)
  print(i); rm(d)
}

save.image("./pvalues_1mil_perms_new_approaches.RData")



attach(res)
# par(mfrow=c(1,2))
# qqplot(-log10(L13.pval.ML.all), -log10(L13.pval.ML.dx)); abline(a=0,b=1,col='red')
# qqplot(-log10(D13.pval.ML.all), -log10(D13.pval.ML.dx)); abline(a=0,b=1,col='red')

res$xi <- -2*(log(L13.pval.ML.dx) + log(D13.pval.ML.dx))
res$p.dx <- 1 - pchisq(res$xi, 4)

res$xi <- -2*(log(L13.pval.ML.all) + log(D13.pval.ML.all))
res$p.all <- 1 - pchisq(res$xi, 4)

attach(res)
par(mfrow=c(1,1))
# qqplot(-log10(p.all), -log10(p.dx)); abline(a=0,b=1,col='red')


t1 <- -log10(res$p.dx[res$type == "neut"])
t2 <- -log10(res$p.dx[res$type != "neut"])



library(patchwork)
p1 <- ggplot(data=res, aes(sample=-log10(p.dx), col=as.factor(type == "neut"))) + geom_qq() + theme_classic() + ggtitle("conditioned")
p2 <- ggplot(data=res, aes(sample=-log10(p.all), col=as.factor(type == "neut"))) + geom_qq() + theme_classic() + ggtitle("all perms")
p1 + p2 + plot_layout(guides = "collect")


qvalue(res$p.dx[res$type != "neut"]) -> d
quantile(d$qvalues)
