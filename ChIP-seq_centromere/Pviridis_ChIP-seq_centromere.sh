# ----- 1. Align ChIP-seq reads to reference genome -----

### H3 ###
# rep1
bowtie2 -p 16 --sensitive-local -x $ref_index -1 $seqDir/7390_I_R1.fastq.gz -2 $seqDir/7390_I_R2.fastq.gz 2> Pviridis_H3_rep1.bowtie2.log | samtools sort -@8 > Pviridis_H3_rep1.bam

# rep2
bowtie2 -p 16 --sensitive-local -x $ref_index -1 $seqDir/7390_J_R1.fastq.gz -2 $seqDir/7390_J_R2.fastq.gz 2> Pviridis_H3_rep2.bowtie2.log | samtools sort -@8 > Pviridis_H3_rep2.bam


### CENH3 ###
# rep1
bowtie2 -p 16 --sensitive-local -x $ref_index -1 $seqDir/7390_G_R1.fastq.gz -2 $seqDir/7390_G_R2.fastq.gz 2> Pviridis_CENH3_rep1.bowtie2.log | samtools sort -@8 > Pviridis_CENH3_rep1.bam

# rep2
bowtie2 -p 16 --sensitive-local -x $ref_index -1 $seqDir/7390_H_R1.fastq.gz -2 $seqDir/7390_H_R2.fastq.gz 2> Pviridis_CENH3_rep2.bowtie2.log | samtools sort -@8 > Pviridis_CENH3_rep2.bam

# ----- merge -----
samtools merge Pviridis_H3.bam Pviridis_H3_rep1.bam Pviridis_H3_rep2.bam
samtools merge Pviridis_CENH3.bam Pviridis_CENH3_rep1.bam Pviridis_CENH3_rep2.bam

samtools index Pviridis_H3.bam
samtools index Pviridis_CENH3.bam

# ----- 2. deeptools -----
#mamba create -n deeptools_local python=3.11 deeptools=3.5.5 "numpy<2" -y
mamba activate deeptools_local

bamCompare -b1 Pviridis_CENH3.bam -b2 Pviridis_H3.bam \
           -o Pviridis_CENH3-H3_RPKM_log2.bigwig -of bigwig \
           -p 4 --ignoreDuplicates \
           --normalizeUsing RPKM --scaleFactorsMethod None --operation log2 2>&1 | tee -a Pviridis_CENH3_H3.bamCompare.log &

mamba deactivate

conda activate bigwigtobedgraph
bigWigToBedGraph Bcaapi_CENH3-H3_RPKM_log2.bigwig Bcaapi_CENH3-H3_RPKM_log2.bedgraph
conda deactivate

# ----- 3. Peak calling using macs3 -----
module load macs3/v3.0.1

nohup macs3 callpeak -t Pviridis_CENH3.bam -c Pviridis_H3.bam \
               -f BAM -g 7.5e+9 -p 1e-3 --broad --broad-cutoff 5e-2 \
               --name Pviridis_CENH3 --outdir Pviridis_CENH3_peaks 2>&1 | tee -a Pviridis_MACS3.log &


cd Pviridis_CENH3_peaks
Rscript Pviridis_CENH3_model.r 

module unload macs3/v3.0.1

SCRIPT="/home/rqiu/Documents/scripts/cenh3_centromere"
grep "Chr" Pviridis_CENH3_peaks.broadPeak > Pviridis_CENH3_peaks.chr.broadPeak
python3 $SCRIPT/merge_macs3_peaks.py Pviridis_CENH3_peaks.chr.broadPeak -o Pviridis_CENH3_peaks.flt.out
# only 58 chromosomes
# add the rest manually -> Pviridis_CENH3_peaks.flt.r.out
python3 $SCRIPT/cen_refine_loc.py -i Pviridis_CENH3_peaks.flt.r.out -b ../Pviridis_CENH3-H3_RPKM_log2.bedgraph -t 0.1 -w 100000 -s 10000 -c 3 -o cen_refine_loc.out 
python3 $SCRIPT/make_plot_region.py cen_refine_loc.out cen_plot_region.txt


# ----- 4. Centromere boundary refinement -----
grep "Chr" Pviridis_CENH3_peaks.broadPeak > Pviridis_CENH3_peaks.chr.broadPeak
python3 merge_macs3_peaks.py Pviridis_CENH3_peaks.chr.broadPeak -o Pviridis_CENH3_peaks.flt.out
# only 58 chromosomes
# add the rest manually -> Pviridis_CENH3_peaks.flt.r.out
python3 cen_refine_loc.py -i Pviridis_CENH3_peaks.flt.r.out -b ../Pviridis_CENH3-H3_RPKM_log2.bedgraph -t 0.1 -w 100000 -s 10000 -c 3 -o cen_refine_loc.out 
python3 make_plot_region.py cen_refine_loc.out cen_plot_region.txt