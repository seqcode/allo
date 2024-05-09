# Allo

A multi-mapped read rescue strategy for gene regulatory analyses.

[Check out our pre-print to learn more!](https://www.biorxiv.org/content/10.1101/2023.09.12.556916v1)

### Releases

As of **v1.1.0**, Allo has neural networks trained for DNase-seq and ATAC-seq under the MACS2 parameters "--nomodel --shift -100 --extsize 200" for ATAC-seq and MACS2 default parameters for DNase-seq. Additionally, Allo now has the option to remove introns as identified by splice junction information in the CIGAR string of an aligned read. This affects the window used to sum uniquely mapped reads. Information below regarding the use of Allo for RNA-seq data processing.

## Installation
### Package managers

*  Bioconda: [![Anaconda-Server Badge](https://anaconda.org/bioconda/allo/badges/version.svg)](https://anaconda.org/bioconda/allo)
[![Anaconda-Server Badge](https://anaconda.org/bioconda/allo/badges/downloads.svg)](https://anaconda.org/bioconda/allo)


*  PyPI: [![PyPI version](https://badge.fury.io/py/bio-allo.svg)](https://badge.fury.io/py/bio-allo)


### Git clone and pip

```
git clone https://github.com/seqcode/allo.git
cd allo
pip install -e .
```

## Usage
### Peak-based applications (ChIP-seq, ATAC-seq, DNase-seq, etc)
#### Pre-processing and alignment
Using Allo requires a few pre-processing steps. In most ChIP-seq, ATAC-seq, and DNase-seq pipelines, the default behavior of aligners is to assign multi-mapped reads to random locations within their mappings without retaining information on the other locations. Both Bowtie1/2 and BWA can be used for single-end. Unfortunately, BWA cannot be used for paired-end reads prior to Allo due to constraints in how it outputs multi-mapped reads. The following arguments should be used:

*Bowtie1*

```
#Single-end
bowtie -x INDEX -q FASTQ -S SAMOUT --best --strata -m 25 -k 25 -p THREADS
#Paired-end
bowtie -x INDEX -1 READ1 -2 READ2 -S SAMOUT --best --strata -m 25 -k 25 -p THREADS
```
*Bowtie2*
```
#Single-end
bowtie2 -x INDEX -q FASTQ -S SAMOUT -k 25 -p THREADS
#Paired-end
bowtie2 -x INDEX -1 READ1 -2 READ2 -S SAMOUT -k 25 --no-mixed --no-discordant -p THREADS
```
*BWA*
```
#BWA is not supported for paired-end processing due to constraints with multi-mapped read reporting
#Single-end
bwa mem -a INDEX FASTQ > SAMOUT
```


Finally, the output of the aligners must be sorted by read name in order to use Allo. If using Bowtie1 in paired-end mode, you may want to use samtools to filter unproperly mapped pairs. It is possible to use Allo with discordant or mixed pairs but Allo will treat them equally to proper pairs. Proceed accordingly. 
```
samtools collate -o ALIGNEROUTPUT_SORT.SAM ALIGNEROUTPUT_FILTER.SAM
```

#### Running Allo
The basic command for Allo:
```
allo ALIGNEROUTPUT_SORT.SAM -seq PAIRED_OR_SINGLE -o OUTPUTNAME
```
Allo also accepts BAM files as input. See other options below..

#### Additional tips
It is recommended to run Allo on both the control and target sequencing files in order to balance out background in the samples. We recommend running Allo using the --random argument on the control file. This generally results in higher confidence peaks.

During each run, Allo will create temporary files as it allocates the data. UM files are reads designated as uniquely mapped (has to be parsed in Bowtie2 or BWA). MM files are unallocated multi-mapped reads. AL files are allocated reads. Checking the size of the AL files during the run will give you an estimate of how many reads have already been allocated at that time.

Very short test files are supplied to make sure Allo runs to completion on your machine. Imports can take a minute so be patient. Using the paired-end example:

```
allo testRunPE.sam -seq pe
```

### RNA-seq application
#### Pre-processing and alignment
Allo is compatible with STAR alignments. We recommend using the "--outFilterType BySJout" argument if you choose to use the "--splice" function in Allo in order to only consider high-quality junctions. An example of a paired-end STAR alignment, keeping up to 25 locations per read, is shown below:
```
STAR --genomeDir GENOMEDIR --readFilesIn fASTQ_1 FASTQ_2 --outSAMtype BAM Unsorted --outSAMmultNmax 25 --outFilterType BySJout --outFileNamePrefix ALIGNEROUTPUT
```

Next, sort your file based on read name:
```
samtools collate -o ALIGNEROUTPUT_SORT.BAM ALIGNEROUTPUT_FILTER.BAM
```

#### Running Allo
Following this, we recommend running Allo on read count only mode as the neural networks available are not trained on RNA-seq profiles. Additionally, the --splice argument can be used if the user would like Allo to splice introns out when summing uniquely mapped reads.
```
allo ALIGNEROUTPUT_SORT.BAM -seq PAIRED_OR_SINGLE -o OUTPUTNAME --readcount --splice
```

#### Downstream analysis
Following the use of Allo, users can utilize FeatureCounts with the argument "-M" which retains multi-mapped reads.
```
featureCounts -a GTF_FILE -o COUNTS.out *.bam -M
```


## Output information
Allo adds a ZA tag to every MMR that is allocated. For reads that are allocated to regions that all contain 0 UMRs (random assignment), a ZZ tag is used instead. This allows users to remove reads that only map to zero UMR regions if they wish. The value within either tag corresponds to the number of places a read/pair mapped to. In order to get only uniquely mapped reads, grep could be used with the -v option to exclude lines with ZA or ZZ tags. On the same note, awk can used to filter reads with a specific number of mapping locations (can also be done with the -max option within Allo). Outside of adding these tags, Allo does not change anything within the read alignment columns for allocated reads.

### Options
| Argument  | Options | Explanation |
| ------------- | ------------- | ------------- |
| -o  | any string | Output file name  |
| -seq | "se" "pe" | Single-end or paired-end sequencing mode, REQUIRED | 
| --mixed | | Use CNN trained on histone ChIP-seq datasets with mixed peaks, narrow by default |
| --dnase | | Use CNN trained on histone DNase-seq datasets, narrow by default |
| --atac | | Use CNN trained on histone ATAC-seq datasets, narrow by default |
| --splice | | Remove introns as identified by splice junctions when summing the uniquely-mapped read counts |
| -p  | any int | Number of processes, 1 by default |
| --keep-unmap |  | Keep unmapped reads and reads that include N in their sequence | 
| --remove-zeros |  | Do not report multi-mapped reads that map to regions with 0 uniquely mapped reads (random assignment) |
| -max | any int | Maximum value for number of locations a read can map to |
| --r2 |  | Use read 2 for allocation procedure instead of read 1 (only for paired-end sequencing) |
| --readcount |  | CNN will not be used in allocation, only read counts |
| --random |  | Reads will be randomly assigned (similar to Bowtie and BWA by default) |
| -s  | any int | Seed for random number generator to keep consistency, 7 by default |
| --ignore |  | Ignore warning about read sorting |
| --parser |  | Allo utility that produces separate files for unique and multi-mapped reads from a SAM or BAM file. Bowtie2 and BWA will output all alignments that meet the given threshold even if one alignment has the highest score. |


## Help and contact information
Please contact Alexis Morrissey (anm5579@psu.edu) with any questions concerning Allo or open an issue here on GitHub. 
