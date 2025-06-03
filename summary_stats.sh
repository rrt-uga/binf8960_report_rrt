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

ml SAMtools/1.18-GCC-12.3.0     #Load SAMtools

mkdir $results/stats/

echo "Sample" > $results/stats/Sample.tsv
grep -r "_1" $data/raw_fastqc_results/multiqc_data/multiqc_fastqc.txt |\
 awk 'BEGIN{FS=OFS="\t"}{print $1}' |\
 sed "s/_1//g" \
 >> $results/stats/Sample.tsv

echo "RawReadCount" > $results/stats/RawReadCount.tsv
grep -r "_1" $data/raw_fastqc_results/multiqc_data/multiqc_fastqc.txt |\
 awk 'BEGIN{FS=OFS="\t"}{print $5}' |\
 sed "s/\.0//g"\
 >> $results/stats/RawReadCount.tsv

echo "TrimmedReadCount" > $results/stats/TrimmedReadCount.tsv
grep -r "_1" $data/trimmed_fastqc_results/multiqc_data/multiqc_fastqc.txt |\
 grep -r "\.paired" - |\
 awk 'BEGIN{FS=OFS="\t"}{print $5}' |\
 sed "s/\.0//g"\
 >> $results/stats/TrimmedReadCount.tsv

echo "AlignedReadCount" > $results/stats/AlignedReadCount.tsv
for bam in $results/bam/*.sorted.bam
do samtools stats\
 -@ 4\
 $bam |\
 grep "reads mapped and paired:" |\
 awk 'BEGIN{FS=OFS="\t"}{print $3/2}' \
 >> $results/stats/AlignedReadCount.tsv
done

echo "VariantCount" > $results/stats/VariantCount.tsv
for vcf in results/vcf/*.vcf
do
 grep -c -v -r "#" $vcf \
 >> $results/stats/VariantCount.tsv
done

paste -d '\t'\
 $results/stats/sample.tsv\
 $results/stats/RawReadCount.tsv\
 $results/stats/TrimmedReadCount.tsv\
 $results/stats/AlignedReadCount.tsv\
 $results/stats/VariantCount.tsv\
 > $results/stats/reads_stat.tsv

ml purge
