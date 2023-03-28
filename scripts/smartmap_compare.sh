#####Bash script to get read depth percent error using Allo and SmartMap. Example for ENCFF631JSV and ENCFF715KYL.
#Requirements: Bowtie2, macs2, samtools, bedtools, cutadapt, SmartMap, Allo, and sm_prep.sh (provided on github).


####Download and alignments####
wget -q https://www.encodeproject.org/files/ENCFF631JSV/@@download/ENCFF631JSV.fastq.gz
wget -q https://www.encodeproject.org/files/ENCFF715KYL/@@download/ENCFF715KYL.fastq.gz
wget -q https://www.encodeproject.org/files/ENCFF873ZZU/@@download/ENCFF873ZZU.fastq.gz
wget -q https://www.encodeproject.org/files/ENCFF611URR/@@download/ENCFF611URR.fastq.gz
gunzip *gz
#Cutting files
cutadapt -j 6 -u -50 -o ENCFF631JSV_50.fastq ENCFF631JSV.fastq
cutadapt -j 6 -u -50 -o ENCFF715KYL_50.fastq ENCFF715KYL.fastq
#Alignments
bowtie2 -x HG38_INDEX -1 ENCFF631JSV.fastq -2 ENCFF715KYL.fastq -S CTCF_MMR.sam -k 25 -p 6 --no-mixed --no-discordant
bowtie2 -x HG38_INDEX -1 ENCFF631JSV_50.fastq -2 ENCFF715KYL_50.fastq -S CTCF_50_MMR.sam -k 25 -p 10 --no-mixed --no-discordant
bowtie2 -x HG38_INDEX -1 ENCFF873ZZU.fastq -2 ENCFF611URR.fastq -S Control_MMR.sam -k 25 -p 6 --no-mixed --no-discordant


#Getting MMRs
samtools collate CTCF_MMR.sam -o CTCF_MMR_sort.sam
python Allo.py CTCF_MMR_sort.sam -seq pe --parser
grep -v "^@" CTCF_MMR_sort.allo.UMR_only.sam | cut -f1 | sort | uniq > temp1
LANG=en_EN join -j 1 <(LANG=en_EN sort -k1 CTCF_50_MMR.sam) <(LANG=en_EN sort -k1 temp1) | gawk 'BEGIN { OFS="\t"}; {$1=$1; print $0}' > temp2
grep "^@" CTCF_50_MMR.sam > head
cat head temp2 > CTCF_map.sam 


##Running smartmap
bash sm_prep.sh -a CTCF_map.sam -o CTCF_SM
SmartMap -g GENOME.INFO -o CTCF_SM CTCF_SM_vf_k_I_X_filt-flag_filt-coord_scores.bed.gz
gunzip CTCF_SM.bedgraph.gz
bedtools sort -i CTCF_SM.bedgraph > CTCF_SM_sort.bedgraph


##Running Allo and making bedgraph
samtools collate CTCF_map.sam -o CTCF_map_sorted.sam
python Allo.py CTCF_map_sorted.sam -seq pe -t 4 -o allo_map.sam
samtools view -Sb allo_map.sam > allo_map.bam
bedtools bamtobed -bedpe -i allo_map.bam | awk -v OFS='\t' '{print $1,$2,$6}' > temp1
bedtools sort -i temp1 > temp2
bedtools genomecov -bga -i temp2 -g GENOME.INFO > CTCF_allo_map.bedgraph
bedtools sort -i CTCF_allo_map.bedgraph > CTCF_allo_map_sort.bedgraph


##Making UMR (ground truth) bedgraph file
grep -v "^@" CTCF_map_sorted.sam | cut -f1 | uniq > temp1
LANG=en_EN join -j 1 <(LANG=en_EN sort -k1 CTCF_MMR_sort.allo.UMR_only.sam) <(LANG=en_EN sort -k1 temp1) | gawk 'BEGIN { OFS="\t"}; {$1=$1; print $0}' > temp2
head -1000 CTCF_map_sorted.sam | grep "^@" > temp3
cat temp3 temp2 > temp4
samtools view -Sb temp4 > temp4.bam
bedtools bamtobed -bedpe -i temp4.bam | awk -v OFS='\t' '{print $1,$2,$6}' > temp5
bedtools sort -i temp5 > temp6
bedtools genomecov -bga -i temp6 -g ~/group/genomes/hg38/hg38.info > CTCF_UMR_map.bedgraph
bedtools sort -i CTCF_UMR_map.bedgraph > CTCF_UMR_sort_map.bedgraph
 

##Getting average read depth across peaks
samtools view -Sb CTCF_MMR_sort.allo.UMR_only.sam > CTCF_MMR_sort.allo.UMR_only.bam
samtools collate Control_MMR.sam -o Control_MMR_sort.sam
python Allo.py Control_MMR_sort.sam -seq pe --parser
samtools view -Sb Control_MMR_sort.allo.UMR_only.sam > Control_MMR_sort.allo.UMR_only.bam
macs2 callpeak -t CTCF_MMR_sort.allo.UMR_only.bam -g hs -f BAMPE -n CTCF_UMR -c Control_MMR_sort.allo.UMR_only.bam
awk -v OFS='\t' '{print $1,$2,$3}' CTCF_UMR_peaks.narrowPeak > temp
bedtools map -a temp -b CTCF_UMR_sort_map.bedgraph -c 4 -o mean -null "0" > temp_actual
bedtools map -a temp -b CTCF_allo_map_sort.bedgraph -c 4 -o mean -null "0" | cut -f4 > temp_allo
bedtools map -a temp -b CTCF_SM_sort.bedgraph -c 4 -o mean -null "0" | cut -f4 > temp_sm
#Getting read depth tables
paste temp_actual temp_allo temp_sm > read_depth_CTCF.txt
rm temp*
samtools view -Sb CTCF_MMR_sort.allo.MMR_only.sam > CTCF_50_MMR_only.bam
bedtools intersect -u -a read_depth_CTCF.txt -b CTCF_50_MMR_only.bam > read_depth_CTCF_MMR.txt