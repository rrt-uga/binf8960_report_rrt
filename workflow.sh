#!/bin/bash
#SBATCH --job-name=binf8960_run_rrt	#Job name
#SBATCH --partition=batch       #Partition name
#SBATCH --ntasks=1      #Task
#SBATCH --cpus-per-task=16      #CPU core number
#SBATCH --mem=64GB     #RAM
#SBATCH --time=48:00:00 #Runtime
#SBATCH --output=%x_%j.out      #Standard output log
#SBATCH --mail-user=rt36111@uga.edu     #Mail
#SBATCH --mail-type=ALL #Mail events (BEGIN, END, FAIL, ALL)

# Genome download
genome_url="ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/017/985/GCA_000017985.1_ASM1798v1/GCA_000017985.1_ASM1798v1_genomic.fna.gz"

# Location of directories
raw_fastq="/work/binf8960/instructor_data/raw_fastq"
data="/scratch/rt36111/binf8960_report_rrt/data"
results="/scratch/rt36111/binf8960_report_rrt/results"

# Copying raw fastq data
cp -r $raw_fastq $data/

# Load FastQC
ml FastQC/0.11.9-Java-11

# FastQC of Raw Data
mkdir $data/raw_fastqc_results/	# Directory for FastQC data of Raw Reads

fastqc $data/raw_fastq/*.fastq.gz\
 -t 16\
 -o $data/raw_fastqc_results/

ml purge

# Load MultiQC
ml MultiQC/1.14-foss-2022a

# MultiQC of Raw Data
multiqc $data/raw_fastqc_results/\
 -o $data/raw_fastqc_results/

ml purge

# Load Trimmomatic
ml Trimmomatic/0.39-Java-13

Trimmomatic="java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar"

# Trim Reads
mkdir $data/trimmed_fastq
for fwd in $data/raw_fastq/*_1.fastq.gz
do
        sample=$(basename $fwd _1.fastq.gz)
        $Trimmomatic PE $data/raw_fastq/${sample}_1.fastq.gz $data/raw_fastq/${sample}_2.fastq.gz\
         $data/trimmed_fastq/${sample}_1.paired.fastq.gz $data/trimmed_fastq/${sample}_1.unpaired.fastq.gz\
         $data/trimmed_fastq/${sample}_2.paired.fastq.gz $data/trimmed_fastq/${sample}_2.unpaired.fastq.gz\
         ILLUMINACLIP:$EBROOTTRIMMOMATIC/adapters/NexteraPE-PE.fa:2:30:10:5:True SLIDINGWINDOW:4:20
done

ml purge

# Load FastQC
ml FastQC/0.11.9-Java-11

# FastQC of Trimmed Data
mkdir $data/trimmed_fastqc_results/	# Directory for FastQC data of Trimmed Reads

fastqc $data/trimmed_fastq/*.fastq.gz\
 -t 16\
 -o $data/trimmed_fastqc_results/

ml purge

# Load MultiQC
ml MultiQC/1.14-foss-2022a

# MultiQC of Trimmed Data
multiqc $data/trimmed_fastqc_results/\
 -o $data/trimmed_fastqc_results/

ml purge

# Download E. coli genome

mkdir $data/genomes	# Directory for storing genome

wget -O $data/genomes/ecoli_rel606.fna.gz\
 $genome_url
gzip -d $data/genomes/ecoli_rel606.fna.gz	# extract genome FASTA

# Index the Genome

ml BWA/0.7.18-GCCcore-13.3.0	# Load BWA
bwa index\
 -p $data/genomes/ecoli_rel606\
 $data/genomes/ecoli_rel606.fna

# Align reads to the genome

mkdir $results/sam	# Directory of SAM files

for fastq in $data/trimmed_fastq/*.paired.fastq.gz
do
	gzip -d $fastq	# extract fastq.gz
done

for fwd in $data/trimmed_fastq/*_1.paired.fastq
do
	sample=$(basename $fwd _1.paired.fastq)
	echo "Aligning $sample"
	rev=$data/trimmed_fastq/${sample}_2.paired.fastq
	bwa mem\
	 -t 16\
	 $data/genomes/ecoli_rel606\
	 $fwd\
	 $rev\
	 > $results/sam/$sample.sam
done

ml purge

# Call variants in the aligned reads

ml SAMtools/1.18-GCC-12.3.0	# Load SAMtools

mkdir $results/bam	# Directory for BAM files
mkdir $results/bcf	# Directory for BCF files
mkdir $results/vcf	# Directory for VCF files

for fwd in $data/trimmed_fastq/*_1.paired.fastq
do
        sample=$(basename $fwd _1.paired.fastq)
	echo "Generate and Sort $sample BAM"
	samtools view\
	 -@ 16\
	 -S -b $results/sam/$sample.sam\
	 > $results/bam/$sample.bam	# SAM to BAM
	samtools sort\
	 -@ 16\
	 -o $results/bam/$sample.sorted.bam\
	 $results/bam/$sample.bam	# Sort BAM
done

ml purge

ml BCFtools/1.18-GCC-12.3.0	# Load BCFtools

for fwd in $data/trimmed_fastq/*_1.paired.fastq
do
        sample=$(basename $fwd _1.paired.fastq)
	echo "Calling variants for $sample"
	bcftools mpileup\
	 --threads 16\
	 -O b\
	 -o $results/bcf/$sample.bcf\
	 -f $data/genomes/ecoli_rel606.fna\
	 $results/bam/$sample.sorted.bam	# Generate BCF
	bcftools call\
	 --threads 16\
	 --ploidy 1\
	 -m -v\
	 -o $results/vcf/$sample.vcf\
	 $results/bcf/$sample.bcf	# Generate VCF to Call Variant
done
