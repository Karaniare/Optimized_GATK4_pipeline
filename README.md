# Optimized_GATK4_pipeline
Part1.sh performes fastq and bam processing and quality check. 
Part2.sh performes variant calling, filtering and annotation
The current version of these scripts should be submitted as slurm jobs.
# Specifying necessary directories
This step is required for both packages prior to submitting your job
ref_dir="/full/path/to/reference/genome/directory"
fastq_dir="/full/path/to/fastq/directory"
bam_dir="/full/path/to/bam/file/directory"
vcf_dir="/full/path/to/vcfs/directory"
unpaired_dir="/users/path/to/keep/unpaired/fastq/directory"
gvcf_dir="/full/path/to/gvcf/directory"
stat_dir="/full/path/to/qc/directory"
