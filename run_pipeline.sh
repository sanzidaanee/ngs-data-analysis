# Part 1: Data Retrieval


# create a new directory 
mkdir ngs_internship
cd ngs_internship
# create sub-directory
mkdir ref_genome
cd ref_genome

# get ref genome
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr20.fa.gz

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz

# unzip file
gunzip chr20.fa.gz
Gunzip chr21.fa.gz

# download SRA TOOL KIT and configure it
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -xzvf sratoolkit.current-ubuntu64.tar.gz
# Rename file
mv sratoolkit.3.1.0-ubuntu64 sra_toolkit # version can be changed

# enter sra_toolkit/bin
cd sra_toolkit/bin
./vdb-config -i #set permission and default folders, go to TOOLS > set to current directory
cd ../../

# download raw data
# using the three samples from bioproject PRJEB62494 ( SRA Accessions: ERR11468775, ERR11468776 and ERR11468777 
# LINK: https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJEB62494
##download the raw data for the accessions ERR11468775, ERR11468776 and ERR11468777 each with 1000000 reads (hint: use -X)
## Use a new directory (raw_data_mapping) to store the above files (Important!!)
## convert into paired end FASTQ files

mkdir raw_data
cd raw_data
../sra_toolkit/bin/fastq-dump --split-files --gzip  ERR11468775 -X 1000000

../sra_toolkit/bin/fastq-dump --split-files --gzip  ERR11468776 -X 1000000

../sra_toolkit/bin/fastq-dump --split-files --gzip  ERR11468777 -X 1000000
# Part 2:  Quality Control and Sequence Alignment

# install fastqc
sudo apt-get -y install fastqc  

##may give errors while installation
##If get any error on fixing issue after giving above command then use this one

sudo apt-get update
sudo apt --fix-broken install
sudo apt-get -y install fastqc

fastqc ERR11468775_* 

fastqc ERR11468776_* 

fastqc ERR11468777_* 

## view QC results
# trimming
# install fastp 
wget http://opengene.org/fastp/fastp
chmod a+x fastp

# run fastp

./fastp -i ERR11468775_1.fastq.gz -o ERR11468775_trimmed_1.fastq.gz -I ERR11468775_2.fastq.gz -O ERR11468775_trimmed_2.fastq.gz --detect_adapter_for_pe -f 10 -g -l 50 -c -h ERR11468775_fastp.html -w 10


./fastp -i ERR11468776_1.fastq.gz -o ERR11468776_trimmed_1.fastq.gz -I ERR11468776_2.fastq.gz -O ERR11468776_trimmed_2.fastq.gz --detect_adapter_for_pe -f 10 -g -l 50 -c -h ERR11468776_fastp.html -w 10

./fastp -i ERR11468777_1.fastq.gz -o ERR11468777_trimmed_1.fastq.gz -I ERR11468777_2.fastq.gz -O ERR11468777_trimmed_2.fastq.gz --detect_adapter_for_pe -f 10 -g -l 50 -c -h ERR11468777_fastp.html -w 10



# install BWA
git clone https://github.com/lh3/bwa.git
cd bwa

## make (use this command) to compile all c file in bwa folder 
./bwa (to check whether it’s running good or not)


# index ref genome with bwa (sam file)
## reference genome downloaded in part 1
## install bwa
sudo apt-get install bwa

## run bwa

bwa index -a bwtsw -p chr20_ref chr20.fa
bwa index -a bwtsw -p chr21_ref chr21.fa

## for each sample
bwa mem chr20_ref raw_data/ERR11468775_trimmed_1.fastq.gz raw_data/ERR11468775_trimmed_2.fastq.gz  -t 10 -o sample1.sam

bwa mem chr20_ref ../raw_data/ERR11468776_trimmed_1.fastq.gz ../raw_data/ERR11468776_trimmed_2.fastq.gz  -t 10 -o sample2.sam

bwa mem chr20_ref ../raw_data/ERR11468777_trimmed_1.fastq.gz ../raw_data/ERR11468777_trimmed_2.fastq.gz  -t 10 -o sample3.sam

bwa mem chr21_ref ../raw_data/ERR11468775_trimmed_1.fastq.gz ../raw_data/ERR11468776_trimmed_2.fastq.gz  -t 10 -o sample4.sam

bwa mem chr21_ref ../raw_data/ERR11468776_trimmed_1.fastq.gz ../raw_data/ERR11468776_trimmed_2.fastq.gz  -t 10 -o sample5.sam

bwa mem chr21_ref ../raw_data/ERR11468777_trimmed_1.fastq.gz ../raw_data/ERR11468776_trimmed_2.fastq.gz  -t 10 -o sample6.sam






# Part 3: Mark Duplication and Recalibration

# create new file 
mkdir bam_file
cd bam_file
## copy the sam and chr20 and chr20 file into this present directory
cp ref_genome/sample1.sam bam_file
## data conversion

##convert sam to bam format

# install samtools

sudo apt-get -y install samtools
samtools view -bo sample1.bam sample1.sam (do it for different samples)
samtools view -bo sample2.bam sample2.sam
samtools view -bo sample1.bam sample3.sam 
samtools view -bo sample2.bam sample4.sam
samtools view -bo sample1.bam sample5.sam 
samtools view -bo sample2.bam sample6.sam

# using GATK
## install GATK latest docker image
## to get known variants file 
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx


## check if u install docker, if not follow previous command
docker images

## run docker in interactive mode
docker run -it -v $PWD:/data broadinstitute/gatk:latest
cd /data/


### data cleaning using GATK
### adding read group info

gatk AddOrReplaceReadGroups -I sample1.bam -O sample1_withRG.bam -ID 1 -LB lib1 -PL ILLUMINA -PU unit1 -SM sample1

gatk AddOrReplaceReadGroups -I sample2.bam -O sample2_withRG.bam -ID 1 -LB lib1 -PL ILLUMINA -PU unit1 -SM sample2


gatk AddOrReplaceReadGroups -I sample3.bam -O sample3_withRG.bam -ID 1 -LB lib1 -PL ILLUMINA -PU unit1 -SM sample3

gatk AddOrReplaceReadGroups -I sample4.bam -O sample4_withRG.bam -ID 1 -LB lib1 -PL ILLUMINA -PU unit1 -SM sample4

gatk AddOrReplaceReadGroups -I sample5.bam -O sample5_withRG.bam -ID 1 -LB lib1 -PL ILLUMINA -PU unit1 -SM sample5

gatk AddOrReplaceReadGroups -I sample6.bam -O sample6_withRG.bam -ID 1 -LB lib1 -PL ILLUMINA -PU unit1 -SM sample6















# sorting reads (sort read in sam and bam file according to genomics coordination)

gatk SortSam -I sample1_withRG.bam -O sorted_sample1.bam -SO coordinate
gatk SortSam -I sample2_withRG.bam -O sorted_sample2.bam -SO coordinate
gatk SortSam -I sample3_withRG.bam -O sorted_sample3.bam -SO coordinate
gatk SortSam -I sample4_withRG.bam -O sorted_sample4.bam -SO coordinate
gatk SortSam -I sample5_withRG.bam -O sorted_sample5.bam -SO coordinate
gatk SortSam -I sample6_withRG.bam -O sorted_sample6.bam -SO coordinate


# overview of the mapping results

samtools flagstat sorted_sample1.bam

samtools flagstat sorted_sample2.bam

samtools flagstat sorted_sample3.bam

samtools flagstat sorted_sample4.bam

samtools flagstat sorted_sample5.bam

samtools flagstat sorted_sample6.bam



# marking duplicates (duplicate reads  provide a matrices)

gatk MarkDuplicates -I sorted_sample1.bam -O sample1_markedDups.bam -M metrics_duplicates1
gatk MarkDuplicates -I sorted_sample2.bam -O sample2_markedDups.bam -M metrics_duplicates2
gatk MarkDuplicates -I sorted_sample3.bam -O sample3_markedDups.bam -M metrics_duplicates3
gatk MarkDuplicates -I sorted_sample4.bam -O sample4_markedDups.bam -M metrics_duplicates4
gatk MarkDuplicates -I sorted_sample5.bam -O sample5_markedDups.bam -M metrics_duplicates5
gatk MarkDuplicates -I sorted_sample6.bam -O sample6_markedDups.bam -M metrics_duplicates6



#BQSR (Base Quality Score Recalibration)
## creating reference genome files 
gatk CreateSequenceDictionary -R chr20.fa
samtools faidx chr20.fa
gatk CreateSequenceDictionary -R chr21.fa
samtools faidx chr21.fa

## run docker in interactive mode
docker run -it -v $PWD:/data broadinstitute/gatk:latest
cd /data/


## running BQSR
## move assembly files into under current working directory (bam_file)

## run docker in interactive mode
docker run -it -v $PWD:/data broadinstitute/gatk:latest
cd /data/
### Check in docker if assembly file have or not


gatk BaseRecalibrator -I sample1_markedDups.bam -R chr20.fa --known-sites Homo_sapiens_assembly38.dbsnp138.vcf -O sample1_recal_data.table

gatk BaseRecalibrator -I sample2_markedDups.bam -R chr20.fa --known-sites Homo_sapiens_assembly38.dbsnp138.vcf -O sample2_recal_data.table

gatk BaseRecalibrator -I sample3_markedDups.bam -R chr20.fa --known-sites Homo_sapiens_assembly38.dbsnp138.vcf -O sample3_recal_data.table

gatk BaseRecalibrator -I sample4_markedDups.bam -R chr21.fa --known-sites Homo_sapiens_assembly38.dbsnp138.vcf -O sample4_recal_data.table

gatk BaseRecalibrator -I sample5_markedDups.bam -R chr21.fa --known-sites Homo_sapiens_assembly38.dbsnp138.vcf -O sample5_recal_data.table

gatk BaseRecalibrator -I sample6_markedDups.bam -R chr21.fa --known-sites Homo_sapiens_assembly38.dbsnp138.vcf -O sample6_recal_data.table

## link to interpret the table
https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR



# Part 4: Variant Calling

## apply BQSR to bam file

gatk ApplyBQSR -R chr20.fa -I sample1_markedDups.bam --bqsr-recal-file sample1_recal_data.table -O sample1_recal.bam

gatk ApplyBQSR -R chr20.fa -I sample2_markedDups.bam --bqsr-recal-file sample2_recal_data.table -O sample2_recal.bam

gatk ApplyBQSR -R chr20.fa -I sample3_markedDups.bam --bqsr-recal-file sample3_recal_data.table -O sample3_recal.bam

gatk ApplyBQSR -R chr21.fa -I sample4_markedDups.bam --bqsr-recal-file sample4_recal_data.table -O sample4_recal.bam

gatk ApplyBQSR -R chr21.fa -I sample5_markedDups.bam --bqsr-recal-file sample5_recal_data.table -O sample5_recal.bam

gatk ApplyBQSR -R chr21.fa -I sample6_markedDups.bam --bqsr-recal-file sample6_recal_data.table -O sample6_recal.bam


#sam tools view (convert /filter)

samtools index sample1_recal.bam

samtools index sample2_recal.bam

samtools index sample3_recal.bam

samtools index sample4_recal.bam

samtools index sample5_recal.bam

samtools index sample6_recal.bam



## calling somatic variants using GATK 
# running tumor only pipeline

gatk Mutect2 -I sample1_recal.bam -R chr20.fa -O sample1.vcf.gz

gatk Mutect2 -I sample2_recal.bam -R chr20.fa -O sample2.vcf.gz

gatk Mutect2 -I sample3_recal.bam -R chr20.fa -O sample3.vcf.gz

gatk Mutect2 -I sample4_recal.bam -R chr21.fa -O sample4.vcf.gz

gatk Mutect2 -I sample5_recal.bam -R chr21.fa -O sample5.vcf.gz

gatk Mutect2 -I sample6_recal.bam -R chr21.fa -O sample6.vcf.gz

## filtering variants

gatk FilterMutectCalls -R chr20.fa -V sample1.vcf.gz -O filtered_sample1.vcf.gz
gatk FilterMutectCalls -R chr20.fa -V sample2.vcf.gz -O filtered_sample2.vcf.gz
gatk FilterMutectCalls -R chr20.fa -V sample3.vcf.gz -O filtered_sample3.vcf.gz
gatk FilterMutectCalls -R chr21.fa -V sample4.vcf.gz -O filtered_sample4.vcf.gz
gatk FilterMutectCalls -R chr21.fa -V sample5.vcf.gz -O filtered_sample5.vcf.gz
gatk FilterMutectCalls -R chr21.fa -V sample6.vcf.gz -O filtered_sample6.vcf.gz


## getting indels

bcftools view --types indels sample1.vcf.gz >> sample1_indels.vcf
bcftools view --types indels sample2.vcf.gz >> sample2_indels.vcf
bcftools view --types indels sample3.vcf.gz >> sample3_indels.vcf
bcftools view --types indels sample4.vcf.gz >> sample4_indels.vcf
bcftools view --types indels sample5.vcf.gz >> sample5_indels.vcf
bcftools view --types indels sample6.vcf.gz >> sample6_indels.vcf

## getting snps

bcftools view --types snps sample1.vcf.gz >> sample1_snps.vcf
bcftools view --types snps sample2.vcf.gz >> sample2_snps.vcf
bcftools view --types snps sample3.vcf.gz >> sample3_snps.vcf
bcftools view --types snps sample4.vcf.gz >> sample4_snps.vcf
bcftools view --types snps sample5.vcf.gz >> sample5_snps.vcf
bcftools view --types snps sample6.vcf.gz >> sample6_snps.vcf
## filtering VCF file

bcftools filter -i '%QUAL>50' sample1.vcf.gz

bcftools filter -i '%QUAL>50' sample2.vcf.gz

bcftools filter -i '%QUAL>50' sample3.vcf.gz

bcftools filter -i '%QUAL>50' sample4.vcf.gz

bcftools filter -i '%QUAL>50' sample5.vcf.gz

bcftools filter -i '%QUAL>50' sample6.vcf.gz


## germline variant calling

### create a folder 
mkdir germline_call
cd germline_call

### copy sample_recal.bam and sample_recal.bam.bai file into current working directory
cp sample1_recal.bam ../germline_cell

## The go to germline_cell folder and run docker

### run docker in interactive mode
docker run -it -v $PWD:/data broadinstitute/gatk:latest
cd /data/
# Then download the following link data 
## germline variant calling
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/hapmap_3.3.hg38.vcf.gz.tbi



gatk HaplotypeCaller -R chr20.fa -I sample1_recal.bam -O sample1_germline.g.vcf.gz -ERC GVCF
gatk HaplotypeCaller -R chr20.fa -I sample2_recal.bam -O sample2_germline.g.vcf.gz -ERC GVCF
gatk HaplotypeCaller -R chr20.fa -I sample3_recal.bam -O sample3_germline.g.vcf.gz -ERC GVCF
gatk HaplotypeCaller -R chr21.fa -I sample4_recal.bam -O sample4_germline.g.vcf.gz -ERC GVCF
gatk HaplotypeCaller -R chr21.fa -I sample5_recal.bam -O sample5_germline.g.vcf.gz -ERC GVCF
gatk HaplotypeCaller -R chr21.fa -I sample6_recal.bam -O sample6_germline.g.vcf.gz -ERC GVCF
## demonstrating with data as sample1: father; sample2: mother; sample3: child

#combine all vcf file into genomic data base

# importing vcf to genomicDB (germline_call folder)
gatk GenomicsDBImport -V sample1_germline.g.vcf.gz -V sample2_germline.g.vcf.gz -V sample3_germline.g.vcf.gz --genomicsdb-workspace-path test_db1 --intervals chr20

gatk GenomicsDBImport -V sample3_germline.g.vcf.gz -V sample4_germline.g.vcf.gz -V sample5_germline.g.vcf.gz --genomicsdb-workspace-path test_db2 --intervals chr21

# genotyping
gatk GenotypeGVCFs -R chr20.fa -V gendb://test_db1 -O try_new.vcf.gz1 

gatk GenotypeGVCFs -R chr21.fa -V gendb://test_db2 -O try_new.vcf.gz2

# check snps

sudo apt update
sudo apt install bcftools -y
bcftools view --types snps try_new.vcf.gz1 
bcftools view --types snps try_new.vcf.gz2


## annotation of variants

### Use VEP of ENSEMBL
### use for GRCH38 version of human genome only
### https://grch37.ensembl.org/Multi/Tools/VEP


