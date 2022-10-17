# Usage
## Required imports
Allo was written with the following imported packages. All other packages used are built into Python. Especially concerning Tensorflow, it is recommended to create a conda environment specific to this version as Allo may not function without it:

joblib	1.2.1

numpy		1.23.3

tensorflow	2.8.1

## Pre-processing
Using Allo requires a few pre-processing steps. In most ChIP pipelines, the default behavior of aligners is to assign multi-mapped reads to random locations within their mappings without retaining information on the other locations. Both Bowtie1/2 and BWA can be used for single-end. Unfortunately, BWA cannot be used for paired-end reads prior to Allo due to constraints in how it outputs multi-mapped reads. The following arguments should be used:

*Bowtie1*

```
#Single-end
bowtie -x INDEX -q FASTQ -S SAMOUT --best --strata -m 50 -k 50 -p THREADS
#Paired-end
bowtie -x INDEX -1 READ1 -2 READ2 -S SAMOUT --best --strata -m 50 -k 50 -p THREADS
```
*Bowtie2*
```
#Single-end
bowtie2 -x INDEX -q FASTQ -S SAMOUT -k 50 -p THREADS
#Paired-end
bowtie2 -x INDEX -1 READ1 -2 READ2 -S SAMOUT -k 50 --no-mixed --no-discordant -p THREADS
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

## Running Allo
The basic command for Allo:
```
python PATH/Allo/Allo.py ALIGNEROUTPUT_SORT.SAM -seq PAIRED_OR_SINGLE -o OUTPUTNAME -m MIXED_OR_NARROW_PEAKS
```
See other options below..

During each run, Allo will create temporary files as it allocates the data. UM files are reads designated as uniquely mapped (has to be parsed in Bowtie2 or BWA). MM files are unallocated multi-mapped reads. AL files are allocated reads. Checking the size of the AL files during the run will give you an estimate of how many reads have already been allocated at that time.

## Post-processing and tips
Allo adds a ZA optional tag to every read that is allocated. The value within the tag corresponds to the number of places a read/pair mapped to. In order to get only uniquely mapped reads, grep could be used with the -v option. On the same note, awk can used to filter reads with a specific number of mapping locations (can also be done with the -max option within Allo). Outside of adding this optional tag, Allo does not change anything within the read alignment columns for allocated reads.

Tip: It is recommended to run Allo on both the control and target sequencing files in order to balance out background in the samples. This generally results in higher confidence peaks.


## Options
| Argument  | Options | Explanation |
| ------------- | ------------- | ------------- |
| -o  | any string | Output file name  |
| -seq | "se" "pe" | Single-end or paired-end sequencing mode, REQUIRED | 
| -m  | "mixed" "narrow" | Use CNN trained on either a narrow peak dataset or a dataset with mixed peaks, narrow by default |
| -t  | any int | Number of threads, 1 by default |
| --keep-unmap |  | Keep unmapped reads and reads that include N in their sequence | 
| --remove-zeros |  | Disregard multi-mapped reads that map to regions with 0 uniquely mapped reads (random assignment) |
| -max | any int | Maximum value for number of locations a read can map to |
| --r2 |  | Use read 2 for allocation procedure instead of read 1 (only for paired-end sequencing) |
| --readcount |  | CNN will not be used in allocation, only read counts |
| --random |  | Reads will be randomly assigned (similar to Bowtie and BWA by default) |
| -s  | any int | Seed for random number generator to keep consistency, 7 by default |
| --ignore |  | Ignore warning about read sorting |


## Test files
Very short test files are supplied to make sure Allo runs to completion on your machine. Imports can take a minute so be patient.

## Contact information
Please contact Alexis Morrissey (anm5579@psu.edu) with any questions or issues concerning Allo.
