# Get list of samples to keep 
mkdir kept_samples
cp VilgalysKlunk_yersinia_pestis-main/part1_aDNA/DATA/SampleNames/keep* ./kept_samples/ ; dos2unix kept_samples/*
ls kept_samples/ > list_of_samples_kept.txt; sed -i s/keep_//g list_of_samples_kept.txt; grep -v '3' list_of_samples_kept.txt > tmp; mv tmp list_of_samples_kept.txt
cat kept_samples/keep_* | sort | uniq > 206_kept_names.txt

# Pair with list of samples that are kept in all 3 conditions
cat kept_samples/keep_london_during_* | sort | uniq -c | grep ' 3 ' | sed 's/.* //g' > kept_samples/keep_london_during_3
cat kept_samples/keep_london_post_* | sort | uniq -c | grep ' 3 ' | sed 's/.* //g' > kept_samples/keep_london_post_3
cat kept_samples/keep_london_pre_* | sort | uniq -c | grep ' 3 ' | sed 's/.* //g' > kept_samples/keep_london_pre_3
cat kept_samples/keep_denmark_post_* | sort | uniq -c | grep ' 3 ' | sed 's/.* //g' > kept_samples/keep_denmark_post_3
cat kept_samples/keep_denmark_pre_* | sort | uniq -c | grep ' 3 ' | sed 's/.* //g' > kept_samples/keep_denmark_pre_3
cat kept_samples/keep_*_3 > kept_samples/samples_all_3.txt

# In addition to the samples kept across all contexts, get the samples that are only kept in some contexts
for f in `cat list_of_samples_kept.txt`; do grep -v -x -f kept_samples/samples_all_3.txt kept_samples/keep_${f} > kept_samples/only_${f} ; done 

# Get bam file for each sample and each enrichment
for f in `cat 206_kept_names.txt`; do 
	samtools view -b -o rescaled_bams_split/exon.$f.bam rescaled_bams/merged.$f.rescaled.bam -L exon.bed; 
	samtools view -b -o rescaled_bams_split/gwas.$f.bam rescaled_bams/merged.$f.rescaled.bam -L immune.bed; 
	samtools view -b -o rescaled_bams_split/neut.$f.bam rescaled_bams/merged.$f.rescaled.bam -L neutral.bed; 
echo $f; done

# Get lists of bam files for individuals kept in all 3 enrichments
cat kept_samples/keep_*_3 > kept_samples/samples_all_3.txt
for f in `cat kept_samples/samples_all_3.txt`; do ls rescaled_bams_split/*.$f.bam >> bam_list.txt; done 
# And for those not in all 3
for f in `cat kept_samples/only*gwas`; do ls rescaled_bams_split/gwas.$f.bam >> bam_list_only_some_enrichments.txt; done 
for f in `cat kept_samples/only*neutral`; do ls rescaled_bams_split/neut.$f.bam >> bam_list_only_some_enrichments.txt; done 
for f in `cat kept_samples/only*exon`; do ls rescaled_bams_split/exon.$f.bam >> bam_list_only_some_enrichments.txt; done 


# Get depth statistics for individuals in all 3
samtools depth -Q 20 -q 20 -H -b 2023_rsid_variants.bed -f bam_list.txt  > samtools_depth.tsv
samtools depth -Q 20 -q 20 -H -b 2023_rsid_variants.bed -f bam_list_only_some_enrichments.txt  > samtools_depth_only_some_enrichments.tsv
# mean coverage of 5.4x for all datasets after same coverage 

# in R, get mean coverage per enrichment type and proportion of reads to keep
noms <- read.delim("kept_samples/samples_all_3.txt", header=F); noms$exon <- noms$gwas <- noms$neut <- NA 
d <- read.delim("samtools_depth.tsv", header=T); sites <- read.delim("2023_rsid_variants.bed", header=F); d$site <- paste(d[,1], d[,2], sep="_"); colnames(sites)[4] <- "site"
d <- merge(sites,d, by="site"); colnames(d) <- gsub("rescaled_bams_split.", "", colnames(d)); colnames(d) <- gsub(".bam", "", colnames(d)); colnames(d)[7] <- "type"
for (i in 1:nrow(noms)) {
	tmp_nom <- noms$V1[i]
  td <- d[[paste("neut",tmp_nom, sep=".")]][d$type == "neut"]; noms$neut[i] <- mean(td)
	td <- d[[paste("gwas",tmp_nom, sep=".")]][d$type == "gwas"]; noms$gwas[i] <- mean(td)
	td <- d[[paste("exon",tmp_nom, sep=".")]][d$type == "exon"]; noms$exon[i] <- mean(td)
}
noms$min <- apply(noms[,2:4],1, min)
noms[,2:4] <- noms$min/noms[,2:4]
write.table(noms, "./coverage_per_type.txt", row.names=F, col.names=F, sep="\t", quote=F)

# repeat for only some enrichments, downsampling to the mean coverage of 5.4x when possible 
d <- read.delim("samtools_depth_only_some_enrichments.tsv", header=T); d$site <- paste(d[,1], d[,2], sep="_")
sites <- read.delim("2023_rsid_variants.bed", header=F); colnames(sites)[4] <- "site"
d <- merge(sites,d, by="site"); colnames(d) <- gsub("rescaled_bams_split.", "", colnames(d)); colnames(d) <- gsub(".bam", "", colnames(d)); colnames(d)[7] <- "type"
noms <- as.data.frame(matrix(ncol=1, colnames(d)[-c(1:9)])); library(plyr); ldply(strsplit(noms$V1,"[.]")) -> v1; colnames(v1) <- c("type","name"); noms <- cbind(noms,v1); rm(v1)
noms$cov <- NA; for (i in 1:nrow(noms)) {
	tmp_nom <- noms$V1[i]; tmp_type <- noms$type[i]
	td <- d[[tmp_nom]][d$type == tmp_type]; noms$cov[i] <- mean(td); rm(tmp_nom, tmp_type, td)
}; rm(i)
noms$ratio <- apply(cbind(1, 5.4/noms$cov),1,min)
write.table(noms, "./coverage_per_enrichment_not_all_3.txt", row.names=F, col.names=F, sep="\t", quote=F)


module load samtools; module load htslib
mkdir subsampled

# Get subsampled bams
# If this sample/enrichment type won't be downsampled ("tmp_cov = 1"), then just copy into the new directory
# For samples shared across all enrichments
for nom in `seq 1 167`; do tmp_nom=`cut -f 1 coverage_per_type.txt | head -$nom | tail -1`;
        tmp_cov=`cut -f 2 coverage_per_type.txt | head -$nom | tail -1`; if [ $tmp_cov = "1" ]; then cp rescaled_bams_split/neut.$tmp_nom.bam subsampled/; else samtools view -s $tmp_cov -b -o subsampled/neut.$tmp_nom.bam rescaled_bams_split/neut.$tmp_nom.bam; fi
        tmp_cov=`cut -f 3 coverage_per_type.txt | head -$nom | tail -1`; if [ $tmp_cov = "1" ]; then cp rescaled_bams_split/gwas.$tmp_nom.bam subsampled/; else samtools view -s $tmp_cov -b -o subsampled/gwas.$tmp_nom.bam rescaled_bams_split/gwas.$tmp_nom.bam; fi
        tmp_cov=`cut -f 4 coverage_per_type.txt | head -$nom | tail -1`; if [ $tmp_cov = "1" ]; then cp rescaled_bams_split/exon.$tmp_nom.bam subsampled/; else samtools view -s $tmp_cov -b -o subsampled/exon.$tmp_nom.bam rescaled_bams_split/exon.$tmp_nom.bam; fi ; done

for nom in `seq 1 66`; do tmp_nom=`cut -f 1 coverage_per_enrichment_not_all_3.txt | head -$nom | tail -1`;
        tmp_cov=`cut -f 5 coverage_per_enrichment_not_all_3.txt | head -$nom | tail -1`;
        if [ $tmp_cov = "1" ]; then cp rescaled_bams_split/$tmp_nom.bam subsampled/; else samtools view -s $tmp_cov -b -o subsampled/$tmp_nom.bam rescaled_bams_split/$tmp_nom.bam; fi ; echo $nom $tmp_cov $tmp_nom; done


# Merge all enrichment types for a sample back into a single bam
for tmp_nom in `cat 206_kept_names.txt`; do samtools merge subsampled/$tmp_nom.bam subsampled/*.$tmp_nom.bam; done
# Cleanup single-enrichment bams
rm subsampled/exon.*bam subsampled/gwas.*bam subsampled/neut.*bam

# Add some files we'll need for later
cp 206_kept_names.txt subsampled/samples.txt
cat kept_samples/keep_london_pre* | sort | uniq > subsampled/keep_london_pre
cat kept_samples/keep_london_post* | sort | uniq > subsampled/keep_london_post
cat kept_samples/keep_london_during* | sort | uniq > subsampled/keep_london_during
cat kept_samples/keep_denmark_pre* | sort | uniq > subsampled/keep_denmark_pre
cat kept_samples/keep_denmark_post* | sort | uniq > subsampled/keep_denmark_post

cd subsampled/
mkdir bams_subsampled/ ; mv L*bam bams_subsampled/

# Trim bam files, add RG
mkdir trimmed_bams; for f in `cat samples.txt`; do ~/Programs/bamUtil-master/bin/bam trimBam bams_subsampled/$f.bam trimmed_bams/$f.bam -L 4 -R 4; echo $f; samtools addreplacerg -r ID:$f -r PL:Illumina -r SM:$f -r LB:$f trimmed_bams/$f.bam -o trimmed_bams/rg.$f.bam; samtools index trimmed_bams/rg.$f.bam; rm trimmed_bams/$f.bam ; done

# call genotypes 
module load java; module load python; module load bcftools; module load vcftools
# gVCF files 
mkdir gVCF; for f in `cat samples.txt`; do ~/Programs/gatk-4.1.4.1/gatk HaplotypeCaller -ERC GVCF -I trimmed_bams/rg.$f.bam  -R ~/hg19/hg19.fa  -O gVCF/$f.g.vcf.gz --intervals ../2023_rsid_variants.bed -mbq 20 --sample-name $f ; done
# joint genotype calling 
ls gVCF/*gz > 02_cohortmap; awk 'BEGIN { FS="/t"; OFS="/t" } { print $1 $1}'  02_cohortmap > tmp4 ; sed 's/gzgVCF/gz \t gVCF/g' tmp4 > tmp2; sed 's/ gVCF/gVCF/g' tmp2 > 02_cohortmap; sed -i 's/^gVCF\///g' 02_cohortmap; sed -i 's/.g.vcf.gz //g' 02_cohortmap; rm tmp2; rm tmp4
ls gVCF/*gz > 03_gvcf.list
mkdir vcf
~/Programs/gatk-4.1.4.1/gatk CombineGVCFs -R ~/hg19/hg19.fa -O ./vcf/tmp.g.vcf.gz -V 03_gvcf.list --intervals ../2023_rsid_variants.bed
~/Programs/gatk-4.1.4.1/gatk GenotypeGVCFs -V ./vcf/tmp.g.vcf.gz -R ~/my_genomes/hg19/hg19.fa -O ./vcf/res.vcf.gz
tabix merged.vcf.gz

# Filter, annotate, and get genotype likelihoods
vcftools --gzvcf ./vcf/res.vcf.gz --mac 1 --max-alleles 2 --remove-indels --minQ 30 --freq --out ./vcf/merged
vcftools --gzvcf ./vcf/res.vcf.gz --mac 1 --max-alleles 2 --remove-indels --minQ 30 --singletons --out ./vcf/merged
vcftools --gzvcf ./vcf/res.vcf.gz --mac 1 --max-alleles 2 --remove-indels --minQ 30 --recode --out ./vcf/merged
vcftools --vcf ./vcf/merged.recode.vcf --exclude-positions ./vcf/merged.singletons --recode --out ./vcf/joint
bcftools annotate -x ^INFO/DP,^FORMAT/GT,^FORMAT/AD,^FORMAT/DP,^FORMAT/GQ,^FORMAT/PL ./vcf/joint.recode.vcf > ./vcf/joint.forgenolik.vcf
mkdir GL
vcftools --vcf ./vcf/joint.forgenolik.vcf --keep keep_london_pre --recode --out temp; vcftools --vcf temp.recode.vcf --012 --out ./GL/london_pre
sed '/^#/d' temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ~/Programs/LCLAE/filtbaboon1b 38 > ./GL/genolik.london_pre.genolik
vcftools --vcf ./vcf/joint.forgenolik.vcf --keep keep_london_post --recode --out temp; vcftools --vcf temp.recode.vcf --012 --out ./GL/london_post
sed '/^#/d' temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ~/Programs/LCLAE/filtbaboon1b 63 > ./GL/genolik.london_post.genolik
vcftools --vcf ./vcf/joint.forgenolik.vcf --keep keep_london_during --recode --out temp; vcftools --vcf temp.recode.vcf --012 --out ./GL/london_during
sed '/^#/d' temp.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ~/Programs/LCLAE/filtbaboon1b 42 > ./GL/genolik.london_during.genolik
vcftools --vcf ./vcf/joint.forgenolik.vcf --keep keep_denmark_pre --recode --out temp2; vcftools --vcf temp2.recode.vcf --012 --out ./GL/denmark_pre
sed '/^#/d' temp2.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ~/Programs/LCLAE/filtbaboon1b 29 > ./GL/genolik.denmark_pre.genolik
vcftools --vcf ./vcf/joint.forgenolik.vcf --keep keep_denmark_post --recode --out temp2; vcftools --vcf temp2.recode.vcf --012 --out ./GL/denmark_post
sed '/^#/d' temp2.recode.vcf | cut -f 1,2,10- | sed -e 's/:/ /g' | sed -e 's/\./999/g' | ~/Programs/LCLAE/filtbaboon1b 34 > ./GL/genolik.denmark_post.genolik

# Save output GLs for local transfer
cp -r GL/ ../subsample_GL/
