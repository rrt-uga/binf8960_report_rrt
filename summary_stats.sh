grep -r "_1" $data/raw_fastqc_results/multiqc_data/multiqc_fastqc.txt |\
 awk 'BEGIN{FS=OFS="\t"}{print $1}' |\
 sed "s/_1//g"

grep -r "_1" $data/raw_fastqc_results/multiqc_data/multiqc_fastqc.txt |\
 awk 'BEGIN{FS=OFS="\t"}{print $5}'

grep -r "_1" $data/trimmed_fastqc_results/multiqc_data/multiqc_fastqc.txt |\
 grep -r "\.paired" - |\
 awk 'BEGIN{FS=OFS="\t"}{print $5}'
