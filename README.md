# Optimized_GATK4_pipeline
Part1.sh performes fastq and bam processing and quality check. 
Part2.sh performes variant calling, filtering and annotation
The current version of these scripts should be submitted as slurm jobs. 
You need to ensure that you have access to enough cores and memories on your server based on the number of samples and the size of your genome of interest before you run the packages. Downsload the files and follow directions below.

## Running optimized GATK4 pipeline
* Download Part1.sh and Part2.sh packages and execute them with arguments
* sh Part1.sh -i [OPTION] -f [OPTION] -r [OPTION] -b [OPTION] -v [OPTION] -u [OPTION] -g [OPTION] -s [OPTION]
* Argument definitions:

-i: input file (sample ID list) in tab separated format
-f: fastq file directory
-r: reference genome and bed file directory
-b: bam file directory
-v : vcf file (or output) directory
-u: directory where unpaired fastq files are kept
-g: gvcf file directory
-s: directory to save qc output files

* sh Part2.sh 
## Creating sample ID list 
This step is required for both packages prior to submitting jobs.

Create a sample ID file and call it ID_list.tsv in your bam file directory.

## Submitting slurm jobs
sbatch Part1.sh

sbatch Part2.sh
