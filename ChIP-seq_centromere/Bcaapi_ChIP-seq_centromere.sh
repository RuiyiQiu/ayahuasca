# ----- 1. Align ChIP-seq reads to reference genome -----

### H3 ###
# rep1
bowtie2 -p 48 --sensitive-local -x $ref_index -1 $dir/7390_O_R1.fastq.gz -2 $dir/7390_O_R2.fastq.gz 2> Bcaapi_TCN_H3_rep1.bowtie2.log | samtools sort -@8 > Bcaapi_TCN_H3_rep1.bam

# rep2
bowtie2 -p 48 --sensitive-local -x $ref_index -1 $dir/7390_P_R1.fastq.gz -2 $dir/7390_P_R2.fastq.gz 2> Bcaapi_TCN_H3_rep2.bowtie2.log | samtools sort -@8 > Bcaapi_TCN_H3_rep2.bam


### CENH3 ###
# rep1
bowtie2 -p 48 --sensitive-local -x $ref_index -1 $dir/7390_M_R1.fastq.gz -2 $dir/7390_M_R2.fastq.gz 2> Bcaapi_TCN_CENH3_rep1.bowtie2.log | samtools sort -@8 > Bcaapi_TCN_CENH3_rep1.bam

# rep2
bowtie2 -p 48 --sensitive-local -x $ref_index -1 $dir/7390_N_R1.fastq.gz -2 $dir/7390_N_R2.fastq.gz 2> Bcaapi_TCN_CENH3_rep2.bowtie2.log | samtools sort -@8 > Bcaapi_TCN_CENH3_rep2.bam

# no filter, merge
samtools merge Bcaapi_TCN_H3.bam Bcaapi_TCN_H3_rep1.bam Bcaapi_TCN_H3_rep2.bam
samtools merge Bcaapi_TCN_CENH3.bam Bcaapi_TCN_CENH3_rep1.bam Bcaapi_TCN_CENH3_rep2.bam

samtools index Bcaapi_TCN_H3.bam
samtools index Bcaapi_TCN_CENH3.bam

# ----- 2. deeptools -----
#mamba create -n deeptools_local python=3.11 deeptools=3.5.5 "numpy<2" -y
mamba activate deeptools_local

bamCompare -b1 Bcaapi_CENH3.bam -b2 Bcaapi_H3.bam \
           -o Bcaapi_CENH3-H3_RPKM_log2.bigwig -of bigwig \
           -p 4 --ignoreDuplicates \
           --normalizeUsing RPKM --scaleFactorsMethod None --operation log2 2>&1 | tee -a Bcaapi_CENH3_H3.bamCompare.log &

mamba deactivate

conda activate bigwigtobedgraph
bigWigToBedGraph Bcaapi_CENH3-H3_RPKM_log2.bigwig Bcaapi_CENH3-H3_RPKM_log2.bedgraph
conda deactivate

# ----- 3. Peak calling using macs3 -----
module load macs3/v3.0.1

nohup macs3 callpeak -t Bcaapi_TCN_CENH3.bam -c Bcaapi_TCN_H3.bam \
               -f BAM -g 1.9e+9 --broad \
               --name Bcaapi_TCN_CENH3 --outdir Bcaapi_TCN_CENH3_peaks 2>&1 | tee -a Bcaapi_TCN_MACS3.log &

cd Bcaapi_TCN_CENH3_peaks
Rscript Bcaapi_TCN_CENH3_model.r 
module unload macs3/v3.0.1

# ----- 4. Centromere boundary refinement -----
python3 merge_macs3_peaks.py Bcaapi_TCN_CENH3_peaks.broadPeak -o Bcaapi_TCN_CENH3_peaks.flt.out -g 100000
grep "chr" Bcaapi_TCN_CENH3_peaks.flt.out > Bcaapi_TCN_CENH3_peaks.chr.flt.out
python3 cen_refine_loc.py -i Bcaapi_TCN_CENH3_peaks.chr.flt.out -b Bcaapi_TCN_CENH3-H3_RPKM_log2.bedgraph -t 0.2 -w 50000 -s 5000 -c 5 -o cen_refine_loc.txt 
python3 make_plot_region.py cen_refine_loc.txt cen_plot_region.txt
