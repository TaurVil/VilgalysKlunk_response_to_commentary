# To be run on cluster. Requires permutations of each loci in their own file. 

# cd ./2022_KlunkVilgalys/Mathieson_Response/extra_permutations
# module load R; R
library(data.table); library(parallel)

########## Load observed data ########## 
load("./fst_estimates.RData")
method="ML"
# Filter and clean up the dataset
info <- info[(info$rsid_number == 1), ]              # Filter for known SNPs (by rsID)
info <- info[!duplicated(info$site), ]               # Remove duplicates by site
info <- info[info$maf.ML >= 0.05, ]                  # Retain sites with maf >= 0.05

# Select candidate sites
info_real <- as.data.frame(info)                     # Convert to data frame
info_real <- info_real[                              # identify candidates based on allele frequency changes consistent with selection
  sign(info_real$delta_L13.ML) == sign(info_real$delta_D13.ML) & 
    sign(info_real$delta_L13.ML) != sign(info_real$delta_L12.ML), ]
# Clean up environment
rm(info, design, design2, min_maf, min_n, drop_samples_missing)

# Check the distribution of the 'type' column
cat("Types of candidate variants:\n")
print(table(info_real$type))

colnames(info_real) <- gsub(".ML", "", colnames(info_real))

########## Calculate permutation-based p-values ########## 
# Define the main function to calculate p-values
# P-value calculations
# abs() is added for the estimates using AF differences rather than Fst
calc_pvals <- function(column, results, perms) {
  sum(abs(perms[[column]]) >= abs(res[[column]])) / sum(perms$site == tmp_snp)
}

perm_type="" # or ".no_replacement"

get_ps <- function(tmp_n, perm_type) {
  tmp_snp <- info_real$site[tmp_n]    # Extract site name
  tmp_file <- paste0("./perms",perm_type,".",tmp_snp,".txt.gz")
  
  res <- info_real[info_real$site == tmp_snp, ]       # Subset for current SNP
  
  # Check if the permutation file exists
  if (!file.exists(tmp_file)) {
    warning(paste("File not found:", tmp_file))
    return(NULL)
  }
  
  # Load permutation data
  extra_perms <- fread(tmp_file)
  # Add permutation statistics
  res$number.perms.all <- nrow(extra_perms)
  
  # p-value calculations
  
  # non-conditioned (compared against all permutations)
  # using Fst
  res$L_1v3.pval.ML.all <- calc_pvals("london.post.pre.fst", res, extra_perms)
  res$D_1v3.pval.ML.all <- calc_pvals("denmark.post.pre.fst", res, extra_perms)
  res$L_12v3.pval.ML.all <- calc_pvals("london.combined.post.fst", res, extra_perms)
  # using the difference in allele frequency instead of Fst
  res$abs_L_1v3.pval.ML.all <- calc_pvals("delta_L13", res, extra_perms)
  res$abs_D_1v3.pval.ML.all <- calc_pvals("delta_D13", res, extra_perms)
  res$abs_L_12v3.pval.ML.all <- calc_pvals("delta_L123", res, extra_perms)

  # conditioned (right patterns, used in the MA)
  conditioned_perms <- extra_perms[
    sign(extra_perms$delta_L13) == sign(extra_perms$delta_D13) & 
      sign(extra_perms$delta_L13) != sign(extra_perms$delta_L12), ]
  res$number.perms.conditioned <- nrow(extra_perms)
  
  # using Fst
  res$L_1v3.pval.ML.conditioned <- calc_pvals("london.post.pre.fst", res, conditioned_perms)
  res$D_1v3.pval.ML.conditioned <- calc_pvals("denmark.post.pre.fst", res, conditioned_perms)
  res$L_12v3.pval.ML.conditioned <- calc_pvals("london.combined.post.fst", res, conditioned_perms)
  # using the difference in allele frequency instead of Fst
  res$abs_L_1v3.pval.ML.conditioned <- calc_pvals("delta_L13", res, conditioned_perms)
  res$abs_D_1v3.pval.ML.conditioned <- calc_pvals("delta_D13", res, conditioned_perms)
  res$abs_L_12v3.pval.ML.conditioned <- calc_pvals("delta_L123", res, conditioned_perms)
  
  # combine effect sizes and directions
  res$L_1v3_with_London_direction <- sum (
    extra_perms$london.post.pre.fst >= res$london.post.pre.fst & 
      sign(extra_perms$delta_L13) != sign(extra_perms$delta_L12)
  )/ nrow(extra_perms)
  
  res$D_1v3_with_same_direction_conditioned_on_london_value <- sum (
    extra_perms$denmark.post.pre.fst >= res$denmark.post.pre.fst.ML & 
      sign(extra_perms$delta_L13) != sign(extra_perms$delta_L12) & 
      extra_perms$london.post.pre.fst >= res$london.post.pre.fst.ML
  )/ sum(extra_perms$london.post.pre.fst >= res$london.post.pre.fst.ML)
  
  res$L_1v3_conditioning_upon_london <- sum (
    extra_perms$london.post.pre.fst >= res$london.post.pre.fst.ML & 
      sign(extra_perms$delta_L13) != sign(extra_perms$delta_L12)
  )/ sum(sign(extra_perms$delta_L13) != sign(extra_perms$delta_L12))
  
  res$D_1v3_conditioning_direction_only <- sum( 
    extra_perms$denmark.post.pre.fst >= res$denmark.post.pre.fst.ML &
      sign(extra_perms$delta_L13) == sign(extra_perms$delta_D13)
  )/ sum(sign(extra_perms$delta_L13) == sign(extra_perms$delta_D13))
  
  # for SI?
  res$all_conditions <- sum (
    extra_perms$london.post.pre.fst >= res$london.post.pre.fst & 
    extra_perms$denmark.post.pre.fst >= res$denmark.post.pre.fst & 
    sign(extra_perms$delta_L13) != sign(extra_perms$delta_L12) &
    sign(extra_perms$delta_L13) == sign(extra_perms$delta_D13)
  ) / nrow(extra_perms)
  
  # Return results
  return(res)
}


# Parallelized execution
results <- mclapply(1:nrow(info_real), function(i) get_ps(i, "ML"), mc.cores = detectCores())

# Combine results
results <- rbindlist(results, fill = TRUE)

# Calculate chi-sq statistics as well
xi <- -2*(log(L_1v3.pval.ML.all) + log(D_1v3.pval.ML.all))
all_res$p.all <- 1 - pchisq(xi, 4)

## THESE ARE FOR THE MAIN FIGURE
xi <- -2*(log(L_1v3.pval.ML.conditioned) + log(D_1v3.pval.ML.conditioned))
all_res$p.conditioned <- 1 - pchisq(xi, 4)

xi <- -2*(log(L_12v3.pval.ML.all) + log(D_1v3.pval.ML.all))
all_res$p.123.all <- 1 - pchisq(xi, 4)

xi <- -2*(log(L_12v3.pval.ML.conditioned) + log(D_1v3.pval.ML.conditioned))
all_res$p.123.conditioned <- 1 - pchisq(xi, 4)

# xi <- -2*(log(L_1v3.pval.ML.conditioned) + log(D_1v3.pval.ML.conditioned) + log(L_2v3.conditioned))
# all_res$p.conditioned.13_and_23 <- 1 - pchisq(xi, 6)


# Save the output
save(results, file = "./DATA/pvalues.RData")


