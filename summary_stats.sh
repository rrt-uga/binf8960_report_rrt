#!/bin/bash
#SBATCH --job-name=summary_stats_generation     #Job name
#SBATCH --partition=batch       #Partition name
#SBATCH --ntasks=1      #Task
#SBATCH --cpus-per-task=4      #CPU core number
#SBATCH --mem=16GB     #RAM
#SBATCH --time=48:00:00 #Runtime
#SBATCH --output=%x_%j.out      #Standard output log
#SBATCH --mail-user=rt36111@uga.edu     #Mail
#SBATCH --mail-type=ALL #Mail events (BEGIN, END, FAIL, ALL)

# Location of directories
data="/scratch/rt36111/binf8960_report_rrt/data"
results="/scratch/rt36111/binf8960_report_rrt/results"

ml SAMtools/1.18-GCC-12.3.0	# Load SAMtools

mkdir $results/stats/	# Create Directory for Stats

# Retrieve Sample Name
# First column of the file contains FASTQ file names
# grep with _1 returns one name for each sample
# awk retrieves only the first column
# Removing _1 using sed from the fastq file name to create the sample name
echo "Sample" > $results/stats/Sample.tsv
grep -r "_1" $data/raw_fastqc_results/multiqc_data/multiqc_fastqc.txt |\
 awk 'BEGIN{FS=OFS="\t"}{print $1}' |\
 sed "s/_1//g" \
 >> $results/stats/Sample.tsv

# Retrieve Raw Read Count
# Fifth column of the file contains read count
# grep with _1 returns one value for each sample
# awk retrieves only the fifth column
# Multiplying the value by 2 to represent both forward and reverse reads
# Removing .0 using sed from the value
echo "RawReadCount" > $results/stats/RawReadCount.tsv
grep -r "_1" $data/raw_fastqc_results/multiqc_data/multiqc_fastqc.txt |\
 awk 'BEGIN{FS=OFS="\t"}{print $5*2}' |\
 sed "s/\.0//g"\
 >> $results/stats/RawReadCount.tsv

# Retrieve Trimmed Read Count
# Fifth column of the file contains read count
# grep with _1 returns one value for each sample
# awk retrieves only the fifth column
# Multiplying the value by 2 to represent both forward and reverse reads
# Removing .0 using sed from the value
echo "TrimmedReadCount" > $results/stats/TrimmedReadCount.tsv
grep -r "_1" $data/trimmed_fastqc_results/multiqc_data/multiqc_fastqc.txt |\
 grep -r "\.paired" - |\
 awk 'BEGIN{FS=OFS="\t"}{print $5*2}' |\
 sed "s/\.0//g"\
 >> $results/stats/TrimmedReadCount.tsv

# Retrieve Aligned Read Count
# Running samtools stats on bam file generates statistics
# Total number of reads mapped and paired is for both forward and reverse files
# grep used to retrieve the line with the value
# The value is in the third column of the line
# awk retrieves only the third column
echo "AlignedReadCount" > $results/stats/AlignedReadCount.tsv
for bam in $results/bam/*.sorted.bam
do samtools stats\
 -@ 4\
 $bam |\
 grep "reads mapped and paired:" |\
 awk 'BEGIN{FS=OFS="\t"}{print $3}' \
 >> $results/stats/AlignedReadCount.tsv
done

# Retrieve Variant Count
# Line count of all lines that has no "#" in front
# grep -c counts the lines
# grep -v is inverse match, so only lines that do not have "#" are returned
echo "VariantCount" > $results/stats/VariantCount.tsv
for vcf in results/vcf/*.vcf
do
 grep -c -v -r "#" $vcf \
 >> $results/stats/VariantCount.tsv
done

# Combine all tsv files
# paste adds files together with tab as delimited
paste -d '\t'\
 $results/stats/Sample.tsv\
 $results/stats/RawReadCount.tsv\
 $results/stats/TrimmedReadCount.tsv\
 $results/stats/AlignedReadCount.tsv\
 $results/stats/VariantCount.tsv\
 > $results/stats/reads_stat.tsv

ml purge
