#####Bash script used to create training set. Example for ENCFF843MGJ. 
#Requirements: Bowtie, macs2, samtools, bedtools, cutadapt, Allo, and cnn_training_gen.sh (provided on github)


#Getting files
wget https://www.encodeproject.org/files/ENCFF968ICP/@@download/ENCFF968ICP.fastq.gz
gunzip *gz
bowtie -q HG38_INDEX_LOCATION ENCFF968ICP.fastq -S Control_UMR.sam --best --strata -m 1 -k 1 --chunkmbs 1024 -p 4 

#Cutting reads
cutadapt -j 4 -u -50 -m 50 -o ENCFF843MGJ_50.fastq ENCFF843MGJ.fastq

#Aligning reads with bowtie
bowtie -q HG38_INDEX_LOCATION ENCFF843MGJ.fastq -S ENCFF843MGJ_UMR.sam --best --strata -m 1 -k 1 --chunkmbs 1024 -p 4 > ENCFF843MGJ_bwt.info
bowtie -q HG38_INDEX_LOCATION ENCFF843MGJ_50.fastq -S ENCFF843MGJ_50_MMR.sam --best --strata -m 25 -k 25 --chunkmbs 1024 -p 4 >> ENCFF843MGJ_bwt.info

#Getting MMRs at peaks
macs2 callpeak -t ENCFF843MGJ_UMR.sam -g hs -f SAM -n ENCFF843MGJ_UMR -c Control_UMR.sam
samtools collate ENCFF843MGJ_50_MMR.sam -o ENCFF843MGJ_50_MMR_sort.sam
python Allo.py ENCFF843MGJ_50_MMR_sort.sam -seq se --parser -t 4

#Getting positive and negative set
grep -v "^@" ENCFF843MGJ_50_MMR_sort.allo.MMR_only.sam | awk -v OFS='\t' '{print $3,$4,$4+50}' > ENCFF843MGJ_temp1
bedtools sort -i ENCFF843MGJ_temp1 > ENCFF843MGJ_temp2
bedtools merge -i ENCFF843MGJ_temp2 > ENCFF843MGJ_temp3
bedtools intersect -u -a ENCFF843MGJ_temp3 -b ENCFF843MGJ_UMR_peaks.narrowPeak > ENCFF843MGJ_temp4
bedtools intersect -v -a ENCFF843MGJ_temp3 -b ENCFF843MGJ_UMR_peaks.narrowPeak > ENCFF843MGJ_temp5

#Running training scripts
python cnn_training_gen.py ENCFF843MGJ_50_MMR_sort.allo.UMR_only.sam ENCFF843MGJ_temp4 ENCFF843MGJ_temp5 ENCFF843MGJ.counts

rm ENCFF843MGJ_temp*
