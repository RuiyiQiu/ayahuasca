# filter 
vcftools --gzvcf merge.vcf.gz \
  --remove-indels \
  --min-alleles 2 --max-alleles 2 \
  --mac 1 \
  --maf 0.05 \
  --max-missing 0.9 \
  --minQ 30 \
  --recode --stdout | bgzip -c > merge_bi_SNP.flt.2.vcf.gz

tabix merge_bi_SNP.flt.2.vcf.gz

###### PCA ######
VCF2PCACluster -InVCF merge_bi_SNP.flt.2.vcf.gz -OutPut $(basename ${VCF} '.vcf.gz')
Rscript PCA_plot.R

##### Admixture #####
#for snmf input:
plink --allow-extra-chr \
  --double-id \
  --out thinned_data \
  --recode \
  --thin 0.01 \
  --vcf merge_bi_SNP.flt.2.vcf.gz

Rscript SNMF.R 

##### Phylogenic network #####
# for splitstree input
~/bin/VCF2Dis-1.55/bin/VCF2Dis \
  -InPut merge_bi_SNP.flt.2.vcf.gz\
  -OutPut $(basename ${VCF} '.vcf.gz').dist