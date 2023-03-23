#####Bash script to get read allocation accuracy. Example for ENCFF873ZZU.
#Requirements: Bowtie, macs2, samtools, bedtools, cutadapt, and Allo


####Control file download, alignment, and allo run####
wget -q https://www.encodeproject.org/files/ENCFF873ZZU/@@download/ENCFF873ZZU.fastq.gz
gunzip ENCFF873ZZU.fastq.gz
bowtie -q HG38_GENOME ENCFF873ZZU.fastq -S ENCFF873ZZU_UMR.sam --best --strata -m 1 -k 1 --chunkmbs 1024 -p 4


####Sample file download, alignment, and allo runs####
wget -q https://www.encodeproject.org/files/ENCFF937ADB/@@download/ENCFF937ADB.fastq.gz
gunzip ENCFF937ADB.fastq.gz
cutadapt -j 4 -u -50 -m 50 -o ENCFF937ADB_50.fastq ENCFF937ADB.fastq
bowtie -q HG38_GENOME ENCFF937ADB.fastq -S ENCFF937ADB_UMR.sam --best --strata -m 1 -k 1 --chunkmbs 1024 -p 4
bowtie -q HG38_GENOME ENCFF937ADB_50.fastq -S ENCFF937ADB_50_MMR.sam --best --strata -m 25 -k 25 --chunkmbs 1024 -p 4
samtools collate ENCFF937ADB_50_MMR.sam -o ENCFF937ADB_50_MMR_sorted.sam
python Allo.py ENCFF937ADB_50_MMR_sorted.sam -seq se -o ENCFF937ADB_50_allo.sam -t 4
python Allo.py ENCFF937ADB_50_MMR_sorted.sam -seq se --readcount -o ENCFF937ADB_50_allo_rc.sam -t 4
python Allo.py ENCFF937ADB_50_MMR_sorted.sam -seq se --random -o ENCFF937ADB_50_allo_rand.sam -t 4


####Running macs2 on UMR sample to get peaks####
macs2 callpeak -t ENCFF937ADB_UMR.sam -g hs -f SAM -n ENCFF937ADB_UMR -c ENCFF873ZZU_UMR.sam


####Getting read accuracy####
#Allo
echo "Allo"
echo "Forward strand reads"
samtools view -F 16 ENCFF937ADB_50_allo.sam | grep "ZA:" | awk -v OFS='\t' '{print $1,$3,$4}' > temp1_for
grep -v "^@" ENCFF937ADB_UMR.sam | awk -v OFS='\t' '{print $1,$3,$4}' | grep -v "*" > temp2
LANG=en_EN join -j 1 -o 1.1,1.2,1.3,2.2,2.3 <(LANG=en_EN sort -k1 temp1_for) <(LANG=en_EN sort -k1 temp2) > temp3
awk -v OFS='\t' '{print $4,$5,$5+100,$2,$3}' temp3 > temp4
bedtools intersect -u -a temp4 -b ENCFF937ADB_UMR_peaks.narrowPeak > temp5
echo "Total MMRs in peaks"
wc -l temp5
echo "Correctly assigned MMRs in peaks"
awk '{if($1==$4 && $2==$5){print}}' temp5 | wc -l
echo "Reverse strand reads"
samtools view -f 16 ENCFF937ADB_50_allo.sam | grep "ZA:" | awk -v OFS='\t' '{print $1,$3,$4}'> temp1_rev
LANG=en_EN join -j 1 -o 1.1,1.2,1.3,2.2,2.3 <(LANG=en_EN sort -k1 temp1_rev) <(LANG=en_EN sort -k1 temp2) > temp3
awk -v OFS='\t' '{print $4,$5,$5+100,$2,$3}' temp3 > temp4
bedtools intersect -u -a temp4 -b ENCFF937ADB_UMR_peaks.narrowPeak > temp5
echo "Total MMRs in peaks"
wc -l temp5
echo "Correctly assigned MMRs in peaks"
awk '{if($1==$4 && $2+50==$5){print}}' temp5 | wc -l

###Read count only
echo "Read count only"
echo "Forward strand reads"
samtools view -F 16 ENCFF937ADB_50_allo_rc.sam | grep "ZA:" | awk -v OFS='\t' '{print $1,$3,$4}' > temp1_for
grep -v "^@" ENCFF937ADB_UMR.sam | awk -v OFS='\t' '{print $1,$3,$4}' | grep -v "*" > temp2
LANG=en_EN join -j 1 -o 1.1,1.2,1.3,2.2,2.3 <(LANG=en_EN sort -k1 temp1_for) <(LANG=en_EN sort -k1 temp2) > temp3
awk -v OFS='\t' '{print $4,$5,$5+100,$2,$3}' temp3 > temp4
bedtools intersect -u -a temp4 -b ENCFF937ADB_UMR_peaks.narrowPeak > temp5
echo "Total MMRs in peaks"
wc -l temp5
echo "Correctly assigned MMRs in peaks"
awk '{if($1==$4 && $2==$5){print}}' temp5 | wc -l
echo "Reverse strand reads"
samtools view -f 16 ENCFF937ADB_50_allo_rc.sam | grep "ZA:" | awk -v OFS='\t' '{print $1,$3,$4}'> temp1_rev
LANG=en_EN join -j 1 -o 1.1,1.2,1.3,2.2,2.3 <(LANG=en_EN sort -k1 temp1_rev) <(LANG=en_EN sort -k1 temp2) > temp3
awk -v OFS='\t' '{print $4,$5,$5+100,$2,$3}' temp3 > temp4
bedtools intersect -u -a temp4 -b ENCFF937ADB_UMR_peaks.narrowPeak > temp5
echo "Total MMRs in peaks"
wc -l temp5
echo "Correctly assigned MMRs in peaks"
awk '{if($1==$4 && $2+50==$5){print}}' temp5 | wc -l

###Random
echo "Read count only"
echo "Forward strand reads"
samtools view -F 16 ENCFF937ADB_50_allo_rand.sam | grep "ZA:" | awk -v OFS='\t' '{print $1,$3,$4}' > temp1_for
grep -v "^@" ENCFF937ADB_UMR.sam | awk -v OFS='\t' '{print $1,$3,$4}' | grep -v "*" > temp2
LANG=en_EN join -j 1 -o 1.1,1.2,1.3,2.2,2.3 <(LANG=en_EN sort -k1 temp1_for) <(LANG=en_EN sort -k1 temp2) > temp3
awk -v OFS='\t' '{print $4,$5,$5+100,$2,$3}' temp3 > temp4
bedtools intersect -u -a temp4 -b ENCFF937ADB_UMR_peaks.narrowPeak > temp5
echo "Total MMRs in peaks"
wc -l temp5
echo "Correctly assigned MMRs in peaks"
awk '{if($1==$4 && $2==$5){print}}' temp5 | wc -l
echo "Reverse strand reads"
samtools view -f 16 ENCFF937ADB_50_allo_rand.sam | grep "ZA:" | awk -v OFS='\t' '{print $1,$3,$4}'> temp1_rev
LANG=en_EN join -j 1 -o 1.1,1.2,1.3,2.2,2.3 <(LANG=en_EN sort -k1 temp1_rev) <(LANG=en_EN sort -k1 temp2) > temp3
awk -v OFS='\t' '{print $4,$5,$5+100,$2,$3}' temp3 > temp4
bedtools intersect -u -a temp4 -b ENCFF937ADB_UMR_peaks.narrowPeak > temp5
echo "Total MMRs in peaks"
wc -l temp5
echo "Correctly assigned MMRs in peaks"
awk '{if($1==$4 && $2+50==$5){print}}' temp5 | wc -l