# Vilgalys, et al. 2023
## Response to Commentary
Code for the response to the commentary by Barton & Santander et al. 

## Enrichment of highly differentiated immune loci

First, we will get allele frequencies and Fst values for each site, using the ML likelihood approach and filtering for known variants. We will then get p-values by comparing candidate immune loci against 200 MAF-matched neutral SNPs. 

```console
01_get_fst_ML.R # calculate fst
02_get_pvals.R # identify outliers by the proportion of 200 maf-matched SNPs their fst level exceeds
```

Second, we will consider the dataset subsampled such that the mean coverage is the same between the neutral, immune GWAS, and immune exon enrichments. 

```console
03_get_subsampled_genoliks.sh
# modified versions of 01 and 02 which use the subsampled genotype likelihoods instead
01b.R
02b.R
```

Third, we consider 10,000 permutations using the strategy proposed by Barton et al. Specifically, we randomly draw individuals for each time point and estimate the enrichment of immune loci which exceed the 99th percentile of neutral variants. 

```console
sbatch --array=1-10000 04.permute_samples.sh
# This saves the Fst values and an enrichment result for each permutation, which can be compared against the real data. 
# It runs modified versions of 01 and 02, with an extra script to get the enrichment for candidate sites which exceed the 99th percentile of neutral loci. 
```

Fourth, we consider 10,000 permutations where site labels are permuted between neutral and candidate immune loci (maintaining the same number of neutral sites). We believe this is a more appropriate null distribution than the time-based permutation suggested by Barton et al. as permutation across time points assumes Fst is 0 for all sites rather than that immune loci are systematically different than putatively neutral loci.  

```console 
# get permutations
05.permute_sites.R
# get enrichment for each permutation
06.permute_sites_enrichment.R
# saves "./DATA/perm_sites_results.RData" which contains a table of 10k results for enrichment statistics beyond the 99th percentile of neutral sites, and the observed enrichment stats
```

Finally, we combine results from these chunks to create Figure 1. 
```console
07.plot_fig1.R
```


## Permutation-based p-values

Although we do not think it is an appropriate null for the overall enrichment of highly differentiated immune loci, the time-based permutation suggested by Barton et al. does provide us with the opportunity to calculate a per site p-value that is independent of MAF-matched neutral loci. Therefore, for each site in our data set with MAF > 5% (n=2264, sites_maf_to_permute.txt), we generate 1 million permutations and calculate the proportion for which the Fst values exceed those of the observed data (requiring that qualitative patterns of genetic differentiation match those consistent with natural selection). 

```console
# Get data files for permutations. Repeat to obtain RData files for neutral, immune gwas, and immune exonic SNPs. 
# Lines 44, 72, and 123 need to be changed for each enrichment
11.get_data_files.R

# Run jobs to get one permutation per site
sbatch --array=1-2264 --mem=4G 12.get_all_perms.sh

# Calculate p-values and attach to fst_estimates.RData from script 01
13.pvalues_1million_permutations.R

# Plot Fig 2A
14.plot_2A.R

# Bootstrap ERAP2 loci for Fig 2B
15.bootstrap_ERAP_loci.R

# Plot Figure 2
panel2A <- readRDS("./fig2a.RData") # print as 8x9"
load("./bootstrap_ERAP_loci.RData") # print plot_erap_loci as 4x9"
```
