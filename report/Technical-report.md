# Technical Report: Variant Calling and Annotation Pipeline
## Whole Genome Sequencing of Tongue Cancer Samples (BioProject PRJEB62494)

**NGS Internship Mini-Project — Team 7**  
**Author:** Sanzida Akhter Anee  
**Supervisors:** NyBerMan Bioinformatics Team / Dr. Sreeram Peela

---

## Table of Contents

1. [Background and Objectives](#1-background-and-objectives)
2. [Dataset Description](#2-dataset-description)
3. [Part 1 — Data Retrieval](#3-part-1--data-retrieval)
4. [Part 2 — Quality Control and Alignment](#4-part-2--quality-control-and-alignment)
5. [Part 3 — BAM Preparation](#5-part-3--bam-preparation)
6. [Part 4 — Variant Calling](#6-part-4--variant-calling)
7. [Part 5 — Variant Filtering and Extraction](#7-part-5--variant-filtering-and-extraction)
8. [Part 6 — Variant Annotation (Ensembl VEP)](#8-part-6--variant-annotation)
9. [Results and Interpretation](#9-results-and-interpretation)
10. [Key Challenges in Variant Calling](#10-key-challenges-in-variant-calling)
11. [Software Versions Reference](#11-software-versions-reference)

---

## 1. Background and Objectives

### What is Variant Calling?

Variant calling is the process of identifying differences between sequenced DNA and a reference genome. These differences — called variants — include:

- **SNPs (Single Nucleotide Polymorphisms):** A single base-pair substitution (e.g., A→G at a specific position).
- **Indels:** Insertions or deletions of one or more base pairs in the genome.
- **Structural Variants (SVs):** Larger genomic rearrangements such as inversions, translocations, and copy number variations.

In cancer genomics, somatic variants (acquired in tumor cells but absent from normal tissue) help identify driver mutations that promote tumor growth, potential therapeutic targets, and biomarkers for prognosis.

Variant calling provides the data foundation for functional gene fine mapping, enables rapid and accurate genome-to-genome comparison, and generates the most extensive set of molecular markers for downstream analysis. It is conceptually simple — find positions where sequencing reads disagree with the reference — but in practice requires a carefully orchestrated multi-step pipeline to distinguish true mutations from sequencing noise, alignment artifacts, and PCR errors.

### Project Objectives

1. Download and process WGS data from tongue cancer samples (BioProject PRJEB62494)
2. Align reads to human reference genome chromosomes 20 and 21
3. Perform GATK Best Practices preprocessing (read group tagging, duplicate marking, BQSR)
4. Call both somatic variants (Mutect2) and germline variants (HaplotypeCaller)
5. Filter, extract, and annotate variants to predict their functional consequences

---

## 2. Dataset Description

### BioProject: PRJEB62494

Whole genome sequencing of tongue cancer samples and a cell line was performed to identify fusion gene translocation breakpoints. Three paired-end WGS samples were sequenced on an Illumina HiSeq 1500 platform.

| SRA Accession | Spots | Bases | Size | GC% | Platform | Published |
|---|---|---|---|---|---|---|
| ERR11468775 | 410.1M | 123.0G | 46.3 GB | 44.0% | Illumina HiSeq 1500 | 2024-02-18 |
| ERR11468776 | 408.4M | 122.5G | 46.2 GB | 43.8% | Illumina HiSeq 1500 | 2024-02-18 |
| ERR11468777 | 394.6M | 118.4G | 46.6 GB | 44.1% | Illumina HiSeq 1500 | 2024-02-18 |

Due to computational resource constraints during the internship, **1,000,000 reads per sample** were downloaded (using the `-X 1000000` flag in `fastq-dump`) rather than the full dataset. The analysis was scoped to chromosomes 20 and 21 of the hg38 reference genome.

---

## 3. Part 1 — Data Retrieval

### 3.1 Reference Genome Download

The reference genome is the baseline against which all sequencing reads are compared to detect variants. We use individual chromosome FASTA files from UCSC for chromosomes 20 and 21 to keep the analysis computationally manageable. `wget` fetches files from UCSC's Golden Path server hosting the official hg38 assembly, and `gunzip` decompresses the `.gz` archives into plain FASTA files required by BWA and GATK.

```bash
mkdir ngs_internship
cd ngs_internship
mkdir ref_genome
cd ref_genome

wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr20.fa.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr21.fa.gz

gunzip chr20.fa.gz
gunzip chr21.fa.gz
```

### 3.2 SRA Toolkit Installation and Configuration

The NCBI Sequence Read Archive stores raw sequencing data. The SRA Toolkit provides `fastq-dump` and `prefetch` utilities to download and convert SRA files into FASTQ format. `vdb-config -i` opens an interactive configuration menu to set the download cache directory and access permissions — this must be run before any data download or the toolkit will fail silently.

```bash
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-ubuntu64.tar.gz
tar -xzvf sratoolkit.current-ubuntu64.tar.gz
mv sratoolkit.3.1.0-ubuntu64 sra_toolkit

cd sra_toolkit/bin
./vdb-config -i   # navigate to TOOLS → set current directory as output
cd ../../
```

### 3.3 Raw Data Download

Paired-end reads generate two files per sample (`_1` and `_2`), representing the forward and reverse reads from each sequenced DNA fragment.

```bash
mkdir raw_data
cd raw_data

../sra_toolkit/bin/fastq-dump --split-files --gzip ERR11468775 -X 1000000
../sra_toolkit/bin/fastq-dump --split-files --gzip ERR11468776 -X 1000000
../sra_toolkit/bin/fastq-dump --split-files --gzip ERR11468777 -X 1000000
```

**Flag explanation:**

| Flag | Meaning |
|---|---|
| `--split-files` | Separates paired-end reads into `_1.fastq.gz` (forward) and `_2.fastq.gz` (reverse) |
| `--gzip` | Compresses output to save disk space |
| `-X 1000000` | Downloads only the first 1 million reads. Full samples are ~410M reads (~46 GB each) |

**Output:**
```
ERR11468775_1.fastq.gz   ERR11468775_2.fastq.gz
ERR11468776_1.fastq.gz   ERR11468776_2.fastq.gz
ERR11468777_1.fastq.gz   ERR11468777_2.fastq.gz
```

---

## 4. Part 2 — Quality Control and Alignment

### 4.1 Quality Control with FastQC

Raw sequencing reads may contain low-quality bases, adapter sequences from the library preparation, overrepresented sequences, or GC bias. FastQC provides a rapid visual summary of read quality so that appropriate trimming parameters can be selected before alignment. Skipping this step risks aligning low-quality reads, which introduces noise into variant calls and increases the false-positive rate downstream.

```bash
sudo apt-get -y install fastqc

fastqc ERR11468775_*
fastqc ERR11468776_*
fastqc ERR11468777_*
```

`fastqc ERR11468775_*` runs FastQC on both `_1.fastq.gz` and `_2.fastq.gz` simultaneously (the `*` wildcard matches both files). FastQC produces an HTML report and a `.zip` archive per file containing:
- Per-base sequence quality scores (Phred scale)
- Per-sequence GC content distribution
- Adapter content detection
- Sequence duplication levels
- N content per base

**Interpreting the report:**
- **Per-base quality:** Should be mostly green (Phred ≥ 20) across the read length. Drops at the 3' end are normal and are addressed in trimming.
- **Adapter content:** Adapters must be removed before alignment — they would otherwise align incorrectly or fail entirely, reducing mapping rate.
- **Per-sequence GC content:** Should follow a normal bell-curve distribution. A bimodal peak suggests contamination from another organism or sample.

### 4.2 Read Trimming with fastp

Trimming removes adapter sequences and low-quality bases at read ends. Adapter sequences would otherwise misalign to the reference. Low-quality 3' ends increase the chance of mismatches being falsely interpreted as variants.

```bash
wget http://opengene.org/fastp/fastp
chmod a+x fastp

./fastp \
  -i ERR11468775_1.fastq.gz \
  -o ERR11468775_trimmed_1.fastq.gz \
  -I ERR11468775_2.fastq.gz \
  -O ERR11468775_trimmed_2.fastq.gz \
  --detect_adapter_for_pe \
  -f 10 \
  -g \
  -l 50 \
  -c \
  -h ERR11468775_fastp.html \
  -w 10

# Repeat for ERR11468776 and ERR11468777
```

**Key fastp flags:**

| Flag | Meaning |
|---|---|
| `--detect_adapter_for_pe` | Auto-detects adapters for paired-end data — no need to specify sequences manually |
| `-f 10` | Hard-clips the first 10 bases from the 5' end (addresses low-quality start cycles) |
| `-g` | Enables polyG trimming — relevant for Illumina platforms where empty spots produce polyG tails |
| `-l 50` | Discards reads shorter than 50 bases after trimming |
| `-c` | Enables base correction for overlapping paired reads |
| `-w 10` | Uses 10 CPU threads for speed |

**Output:** Trimmed FASTQ files and an HTML report per sample showing before/after quality metrics.

### 4.3 Read Alignment with BWA-MEM

Alignment maps each sequencing read back to its position of origin in the reference genome, producing a SAM file. This positional information is the foundation for identifying which bases differ between the sample and the reference at each genomic coordinate. BWA-MEM is the industry standard for aligning Illumina short reads and uses a Maximal Exact Match seeding strategy to handle mismatches, gaps, and spliced alignments efficiently.

#### Index the Reference Genome

BWA cannot search through millions of bases sequentially for every read — it would take days per sample. Indexing pre-builds a Burrows-Wheeler Transform (BWT) data structure allowing sub-second lookup of where a read could align in the genome.

```bash
sudo apt-get install bwa

bwa index -a bwtsw -p chr20_ref chr20.fa
bwa index -a bwtsw -p chr21_ref chr21.fa
```

- `-a bwtsw`: BWT-SW algorithm, appropriate for genomes > 2 GB total.
- `-p chr20_ref`: Sets the prefix name for the index files (generates `.bwt`, `.pac`, `.ann`, `.amb`, `.sa`).

#### Align Reads to Reference

Six alignment jobs are run: three samples × two chromosomes (chr20 and chr21), producing six SAM files.

```bash
# Chromosome 20 alignments
bwa mem chr20_ref raw_data/ERR11468775_trimmed_1.fastq.gz \
                  raw_data/ERR11468775_trimmed_2.fastq.gz \
                  -t 10 -o sample1.sam

bwa mem chr20_ref raw_data/ERR11468776_trimmed_1.fastq.gz \
                  raw_data/ERR11468776_trimmed_2.fastq.gz \
                  -t 10 -o sample2.sam

bwa mem chr20_ref raw_data/ERR11468777_trimmed_1.fastq.gz \
                  raw_data/ERR11468777_trimmed_2.fastq.gz \
                  -t 10 -o sample3.sam

# Chromosome 21 alignments
bwa mem chr21_ref raw_data/ERR11468775_trimmed_1.fastq.gz \
                  raw_data/ERR11468775_trimmed_2.fastq.gz \
                  -t 10 -o sample4.sam

bwa mem chr21_ref raw_data/ERR11468776_trimmed_1.fastq.gz \
                  raw_data/ERR11468776_trimmed_2.fastq.gz \
                  -t 10 -o sample5.sam

bwa mem chr21_ref raw_data/ERR11468777_trimmed_1.fastq.gz \
                  raw_data/ERR11468777_trimmed_2.fastq.gz \
                  -t 10 -o sample6.sam
```

- `bwa mem`: Seeds alignments with long exact matches, then extends them to handle mismatches and gaps. Optimal for reads 70 bp–1 Mbp.
- `-t 10`: Uses 10 threads to parallelize alignment across CPU cores.
- Each line in the SAM output records chromosome, start position, strand, CIGAR string (alignment encoding), mapping quality, and the read sequence.

---

## 5. Part 3 — BAM Preparation

### 5.1 SAM to BAM Conversion

SAM files are human-readable plain text — large and slow to process. BAM is the binary compressed equivalent, occupying ~5× less disk space and enabling indexed random access to any genomic region. All downstream GATK tools require BAM input.

```bash
sudo apt-get -y install samtools

samtools view -bo sample1.bam sample1.sam
samtools view -bo sample2.bam sample2.sam
samtools view -bo sample3.bam sample3.sam
samtools view -bo sample4.bam sample4.sam
samtools view -bo sample5.bam sample5.sam
samtools view -bo sample6.bam sample6.sam
```

- `-b`: Outputs in BAM binary format.
- `-o`: Specifies the output filename.

### 5.2 Add Read Groups

GATK tools require read group (RG) metadata embedded in BAM files. Read groups identify which sequencer, flow cell, lane, sample, and library a read came from. Without this information GATK's variant callers cannot correctly separate reads from different samples and will throw an error at runtime. Read group tags also enable BQSR to build error models per lane and per sequencing cycle.

```bash
docker run -it -v $PWD:/data broadinstitute/gatk:latest
cd /data/

gatk AddOrReplaceReadGroups -I sample1.bam -O sample1_withRG.bam \
  -ID 1 -LB lib1 -PL ILLUMINA -PU unit1 -SM sample1

gatk AddOrReplaceReadGroups -I sample2.bam -O sample2_withRG.bam \
  -ID 1 -LB lib1 -PL ILLUMINA -PU unit1 -SM sample2

# ... repeat for samples 3–6
```

**Read group fields:**

| Field | Flag | Meaning |
|---|---|---|
| ID | `-ID` | Unique read group identifier for this set of reads |
| LB | `-LB` | Library identifier — duplicates are only flagged within the same library |
| PL | `-PL` | Sequencing platform (ILLUMINA, PACBIO, etc.) |
| PU | `-PU` | Platform unit — typically flow cell + lane + barcode |
| SM | `-SM` | Sample name — must be consistent across all read groups from the same biological sample |

### 5.3 Coordinate Sorting

BWA outputs reads in the order they were aligned, which is approximately random. GATK's downstream tools — MarkDuplicates, BQSR, and the variant callers — all require reads sorted by genomic coordinates so that random-access seeking to any chromosomal position is possible.

```bash
gatk SortSam -I sample1_withRG.bam -O sorted_sample1.bam -SO coordinate
gatk SortSam -I sample2_withRG.bam -O sorted_sample2.bam -SO coordinate
gatk SortSam -I sample3_withRG.bam -O sorted_sample3.bam -SO coordinate
gatk SortSam -I sample4_withRG.bam -O sorted_sample4.bam -SO coordinate
gatk SortSam -I sample5_withRG.bam -O sorted_sample5.bam -SO coordinate
gatk SortSam -I sample6_withRG.bam -O sorted_sample6.bam -SO coordinate
```

- `-SO coordinate`: Sorts reads by chromosome and start position. The alternative `queryname` (sort by read name) is not suitable for variant calling.

### 5.4 Alignment Statistics with SAMtools flagstat

Before proceeding with expensive downstream steps, verifying alignment quality is essential. `flagstat` summarises the SAM/BAM flags for each read and reports mapping rates, pairing statistics, and duplicate counts. Low mapping rates (<70%) would indicate problems with read quality, reference mismatch, or sample contamination and would require investigation before proceeding.

```bash
samtools flagstat sorted_sample1.bam
samtools flagstat sorted_sample2.bam
samtools flagstat sorted_sample3.bam
samtools flagstat sorted_sample4.bam
samtools flagstat sorted_sample5.bam
samtools flagstat sorted_sample6.bam
```

**Example output interpretation:**
```
1000000 + 0 in total (QC-passed reads + QC-failed reads)
953214 + 0 mapped (95.32%)          ← High mapping rate = good alignment
1000000 + 0 paired in sequencing
940180 + 0 properly paired (94.02%) ← High proper pair rate = expected insert size
```

### 5.5 Mark PCR Duplicates

During PCR amplification in library preparation, some DNA fragments are copied multiple times. These duplicate reads all align to the exact same start and end position, creating artificial over-representation of those loci. If not accounted for, duplicates inflate coverage and cause variant callers to report incorrectly high confidence in variants that are actually PCR artifacts. GATK MarkDuplicates identifies duplicates by comparing the 5' positions of read pairs and flags them in the BAM — they are not deleted, just marked with a bitwise flag so downstream tools can ignore them.

```bash
gatk MarkDuplicates -I sorted_sample1.bam -O sample1_markedDups.bam -M metrics_duplicates1
gatk MarkDuplicates -I sorted_sample2.bam -O sample2_markedDups.bam -M metrics_duplicates2
gatk MarkDuplicates -I sorted_sample3.bam -O sample3_markedDups.bam -M metrics_duplicates3
gatk MarkDuplicates -I sorted_sample4.bam -O sample4_markedDups.bam -M metrics_duplicates4
gatk MarkDuplicates -I sorted_sample5.bam -O sample5_markedDups.bam -M metrics_duplicates5
gatk MarkDuplicates -I sorted_sample6.bam -O sample6_markedDups.bam -M metrics_duplicates6
```

- `-M metrics_duplicatesN`: Writes a metrics file reporting the percentage of duplicate reads. A typical WGS library has 5–20% duplicates; higher rates may indicate over-amplification.

Two types of duplicates are handled:
- **PCR duplicates:** The same fragment amplified multiple times during library preparation.
- **Optical duplicates:** A single cluster on the flow cell incorrectly detected as multiple clusters by the optical sensor.

### 5.6 Base Quality Score Recalibration (BQSR)

Illumina sequencers assign a Phred quality score to each base call (Q20 = 99% correct; Q30 = 99.9% correct). However, these scores are not perfectly calibrated — systematic errors occur due to sequencing cycle number, flanking base context, and machine-specific artifacts. BQSR uses machine learning with a known-variants database (dbSNP 138) to detect and correct these patterns, producing more accurate per-base quality scores that the variant caller uses in its statistical model. Reference: https://gatk.broadinstitute.org/hc/en-us/articles/360035890531

#### Prerequisite: Create Reference Index Files

GATK requires two types of index files before it will run. The `.dict` dictionary maps chromosome names to their lengths and is used to validate consistency between the BAM and reference files. The `.fai` FASTA index allows random access to any sequence position in the FASTA file. GATK refuses to run without both.

```bash
gatk CreateSequenceDictionary -R chr20.fa    # creates chr20.dict
samtools faidx chr20.fa                       # creates chr20.fa.fai

gatk CreateSequenceDictionary -R chr21.fa
samtools faidx chr21.fa
```

#### Step 1: Generate Recalibration Table

```bash
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

gatk BaseRecalibrator \
  -I sample1_markedDups.bam \
  -R chr20.fa \
  --known-sites Homo_sapiens_assembly38.dbsnp138.vcf \
  -O sample1_recal_data.table

# Repeat for samples 2–6 (use chr21.fa for samples 4–6)
```

- `--known-sites`: A VCF of known variant positions (dbSNP 138). BQSR uses this to distinguish true SNPs (real mismatches to the reference) from sequencing errors. Only mismatches at positions *not* in the known-sites VCF are modelled as errors, preventing real polymorphisms from being downweighted.

#### Step 2: Apply Recalibration

```bash
gatk ApplyBQSR \
  -R chr20.fa \
  -I sample1_markedDups.bam \
  --bqsr-recal-file sample1_recal_data.table \
  -O sample1_recal.bam

# Repeat for all samples
```

This rewrites the quality scores in the BAM file based on the correction model from Step 1. The output (`sample1_recal.bam`) has empirically calibrated per-base quality scores that are more accurate than the original instrument-assigned scores.

### 5.7 BAM Indexing

A BAM index (`.bai` file) allows rapid random access to reads at any genomic position without scanning the entire file. GATK variant callers require indexed BAMs to efficiently retrieve only the reads within each genomic window they analyse.

```bash
samtools index sample1_recal.bam
samtools index sample2_recal.bam
samtools index sample3_recal.bam
samtools index sample4_recal.bam
samtools index sample5_recal.bam
samtools index sample6_recal.bam
```

---

## 6. Part 4 — Variant Calling

### 6.1 Somatic Variant Calling with GATK Mutect2

Mutect2 is GATK's dedicated somatic variant caller. It uses a Bayesian statistical model to distinguish somatic mutations (acquired in the tumor) from germline variants and sequencing errors. In tumor-only mode (used here, without a matched normal sample), Mutect2 calls all variants that are not clearly germline, then relies on FilterMutectCalls to remove likely artifacts. Mutect2 is deliberately sensitive — it over-calls initially because missing a real mutation is worse than retaining a false positive that will be removed in the filtering step.

```bash
gatk Mutect2 -I sample1_recal.bam -R chr20.fa -O sample1.vcf.gz
gatk Mutect2 -I sample2_recal.bam -R chr20.fa -O sample2.vcf.gz
gatk Mutect2 -I sample3_recal.bam -R chr20.fa -O sample3.vcf.gz

gatk Mutect2 -I sample4_recal.bam -R chr21.fa -O sample4.vcf.gz
gatk Mutect2 -I sample5_recal.bam -R chr21.fa -O sample5.vcf.gz
gatk Mutect2 -I sample6_recal.bam -R chr21.fa -O sample6.vcf.gz
```

**VCF output format per variant record:**
- Chromosome, position, reference allele, alternative allele
- FILTER field (PASS or a named failure reason)
- INFO field (allele frequency, read depth, strand bias, etc.)
- FORMAT/SAMPLE fields (genotype, allele depths, per-sample quality)

### 6.2 Somatic Variant Filtering

The initial Mutect2 call set contains both real somatic mutations and artifacts from mapping errors, strand artifacts, or oxidative DNA damage (OxoG). FilterMutectCalls applies a machine learning model to assign each variant a probability of being a true somatic mutation, setting FILTER=PASS for high-confidence variants.

```bash
gatk FilterMutectCalls -R chr20.fa -V sample1.vcf.gz -O filtered_sample1.vcf.gz
gatk FilterMutectCalls -R chr20.fa -V sample2.vcf.gz -O filtered_sample2.vcf.gz
gatk FilterMutectCalls -R chr20.fa -V sample3.vcf.gz -O filtered_sample3.vcf.gz
gatk FilterMutectCalls -R chr21.fa -V sample4.vcf.gz -O filtered_sample4.vcf.gz
gatk FilterMutectCalls -R chr21.fa -V sample5.vcf.gz -O filtered_sample5.vcf.gz
gatk FilterMutectCalls -R chr21.fa -V sample6.vcf.gz -O filtered_sample6.vcf.gz
```

### 6.3 Germline Variant Calling with HaplotypeCaller

HaplotypeCaller is GATK's germline variant caller — it identifies inherited variants present in all cells of an individual. The three samples are used in a trio germline calling scenario (sample1: father, sample2: mother, sample3: child), enabling identification of de novo mutations (variants in the child absent from both parents), allele phasing by parental origin, and Mendelian inheritance checking.

#### Step 1: Per-Sample GVCF Mode

```bash
gatk HaplotypeCaller -R chr20.fa -I sample1_recal.bam \
  -O sample1_germline.g.vcf.gz -ERC GVCF

gatk HaplotypeCaller -R chr20.fa -I sample2_recal.bam \
  -O sample2_germline.g.vcf.gz -ERC GVCF

gatk HaplotypeCaller -R chr20.fa -I sample3_recal.bam \
  -O sample3_germline.g.vcf.gz -ERC GVCF

# Repeat for chr21 (samples 4–6)
```

- `-ERC GVCF`: Emit Reference Confidence mode — outputs a GVCF recording both variant sites AND reference-confidence blocks. This is required for joint genotyping so that "no variant called here" is distinguished from "no data here."

#### Step 2: Combine GVCFs into GenomicsDB

GenomicsDB is a columnar database optimised for genomic data. Combining individual GVCFs enables joint genotyping across all samples simultaneously, which increases statistical power for detecting rare variants and improves genotype accuracy at each site.

```bash
gatk GenomicsDBImport \
  -V sample1_germline.g.vcf.gz \
  -V sample2_germline.g.vcf.gz \
  -V sample3_germline.g.vcf.gz \
  --genomicsdb-workspace-path test_db1 \
  --intervals chr20

gatk GenomicsDBImport \
  -V sample4_germline.g.vcf.gz \
  -V sample5_germline.g.vcf.gz \
  -V sample6_germline.g.vcf.gz \
  --genomicsdb-workspace-path test_db2 \
  --intervals chr21
```

#### Step 3: Joint Genotyping

```bash
gatk GenotypeGVCFs -R chr20.fa -V gendb://test_db1 -O try_new.vcf.gz1
gatk GenotypeGVCFs -R chr21.fa -V gendb://test_db2 -O try_new.vcf.gz2
```

This applies a joint model across all samples at every site with evidence of variation, producing a multi-sample VCF with genotype calls (0/0 = homozygous ref, 0/1 = heterozygous, 1/1 = homozygous alt) for each individual at each variant site.

---

## 7. Part 5 — Variant Filtering and Extraction

### 7.1 Extract SNPs and Indels Separately

SNPs and indels have different error profiles and are filtered using different criteria. Separating them enables type-appropriate filtering and simplifies downstream analysis workflows.

```bash
# Extract indels
bcftools view --types indels sample1.vcf.gz >> sample1_indels.vcf
bcftools view --types indels sample2.vcf.gz >> sample2_indels.vcf
bcftools view --types indels sample3.vcf.gz >> sample3_indels.vcf
# ... samples 4–6

# Extract SNPs
bcftools view --types snps sample1.vcf.gz >> sample1_snps.vcf
bcftools view --types snps sample2.vcf.gz >> sample2_snps.vcf
bcftools view --types snps sample3.vcf.gz >> sample3_snps.vcf
# ... samples 4–6
```

- `--types snps`: Filters the VCF to output only SNP records.
- `--types indels`: Filters to only insertion/deletion records.
- `>>`: Appends output — the VCF header from the input is preserved in the output file.

### 7.2 Quality Filtering

Not all called variants are equally reliable. Variants with low QUAL scores are more likely to be false positives caused by sequencing errors or alignment artifacts. Applying a threshold of QUAL > 50 retains only high-confidence calls.

```bash
bcftools filter -i '%QUAL>50' sample1.vcf.gz
bcftools filter -i '%QUAL>50' sample2.vcf.gz
bcftools filter -i '%QUAL>50' sample3.vcf.gz
bcftools filter -i '%QUAL>50' sample4.vcf.gz
bcftools filter -i '%QUAL>50' sample5.vcf.gz
bcftools filter -i '%QUAL>50' sample6.vcf.gz
```

**QUAL score interpretation:**
- QUAL = 50: ~99.999% probability the variant is real (10⁻⁵ probability of error)
- QUAL = 20: 99% probability of being real
- QUAL > 50 is a practical hard-filter threshold for high-confidence WGS variant calls

The gold-standard filtering approach for large cohorts is GATK's Variant Quality Score Recalibration (VQSR), which trains a Gaussian mixture model on known variant databases to distinguish true variants from artifacts. VQSR requires thousands of variant sites and is computationally intensive, so hard filtering was applied here as a practical alternative.

### 7.3 View Germline Variants

```bash
bcftools view --types snps try_new.vcf.gz1    # chr20 germline SNPs
bcftools view --types snps try_new.vcf.gz2    # chr21 germline SNPs
```

---

## 8. Part 6 — Variant Annotation (Ensembl VEP)

A VCF file tells you *where* a variant is in the genome. Annotation translates those coordinates into biological meaning. Ensembl's Variant Effect Predictor (VEP) provides:

- Which gene and transcript are affected
- Whether the variant falls in a coding, regulatory, or intergenic region
- The codon change (e.g., GAA→GAG, both encode Glutamic acid = synonymous)
- The amino acid change (e.g., p.Glu105Lys = missense mutation)
- The predicted functional consequence (synonymous, missense, stop_gained, frameshift, splice_region)
- Population allele frequencies from gnomAD/1000 Genomes to distinguish rare from common variants

### Web Interface Usage

1. Go to https://grch37.ensembl.org/Multi/Tools/VEP (GRCh37) or https://asia.ensembl.org (Asia mirror for lower latency)
2. Upload the filtered VCF file (e.g., `filtered_sample1.vcf.gz`)
3. Select species: *Homo sapiens*
4. Select assembly: GRCh38 (or GRCh37 depending on reference used)
5. Click Run

### VEP Output Fields

| Field | Description |
|---|---|
| Consequence | Functional effect (missense_variant, synonymous_variant, stop_gained, etc.) |
| SYMBOL | Gene name (e.g., TP53, EGFR) |
| Protein_position | Amino acid position in the protein |
| Amino_acids | Reference / alternative amino acid (e.g., E/K = Glu→Lys) |
| Codons | Reference / alternative codon |
| IMPACT | Severity: HIGH, MODERATE, LOW, MODIFIER |
| AF | Allele frequency in reference population |
| SIFT | Tolerance prediction: tolerated / deleterious |
| PolyPhen | Pathogenicity: benign / possibly_damaging / probably_damaging |

---

## 9. Results and Interpretation

### 9.1 Variant Calling Overview

Variants were successfully called for all three tongue cancer samples on both chromosomes 20 and 21 using two parallel pipelines — somatic (Mutect2 tumor-only) and germline (HaplotypeCaller trio). Six VCF files were generated per pipeline, corresponding to the three samples aligned to each of the two chromosomes. The pipeline follows the GATK sensitivity-first philosophy: Mutect2 is deliberately permissive in its initial calls, and FilterMutectCalls subsequently removes artifact-rich variants using a machine learning model, so only FILTER=PASS variants pass to annotation.

### 9.2 Per-Sample Variant Output

**Chromosome 20** (samples 1–3):
- Each sample produced an independent filtered VCF (`filtered_sample[1-3].vcf.gz`)
- SNPs extracted to `sample[1-3]_snps.vcf` | Indels extracted to `sample[1-3]_indels.vcf`
- Quality-filtered variants (QUAL > 50) represent the high-confidence call set

**Chromosome 21** (samples 4–6, same three accessions aligned to chr21):
- SNPs: `sample[4-6]_snps.vcf` | Indels: `sample[4-6]_indels.vcf`
- Germline joint VCF: `try_new.vcf.gz1` (chr20) and `try_new.vcf.gz2` (chr21)

### 9.3 Common vs. Sample-Specific Variants

Comparing the variant call sets across the three samples reveals two biologically distinct classes:

**Shared variants (present in ≥2 samples):** These are likely common germline polymorphisms or population-level SNPs present in the dbSNP database. They are generally lower priority for cancer driver analysis unless they occur in known cancer genes. However, shared heterozygous variants in DNA repair genes (e.g., BRCA2 on chromosome 13, or RECQL4 nearby) could represent inherited susceptibility alleles.

**Sample-specific variants (present in only one sample):** These are strong candidates for somatic mutations or rare germline variants unique to that individual. In tongue squamous cell carcinoma, sample-specific mutations affecting cell cycle control, DNA damage response, or receptor tyrosine kinase signalling pathways are of highest biological relevance and are the primary targets for functional follow-up.

### 9.4 Germline Trio Analysis Results

The HaplotypeCaller pipeline was run in trio configuration: sample1 (father), sample2 (mother), sample3 (child). Joint genotyping via GenotypeGVCFs produced a multi-sample VCF where each variant site carries a genotype call for all three individuals simultaneously.

This configuration enables several analyses that single-sample calling cannot support:

- **De novo mutation identification:** Any variant where the child (sample3) carries a non-reference allele (0/1 or 1/1) while both parents are homozygous reference (0/0) is a candidate de novo mutation. De novo mutations are of particular interest because they arise freshly in the germline and may contribute to cancer predisposition.
- **Mendelian inheritance checking:** Variants that violate expected inheritance rules (e.g., child homozygous for an allele absent in both parents) flag potential genotyping errors or complex genomic events such as uniparental disomy.
- **Allele phasing:** Determining whether two heterozygous variants in the same gene are in cis (on the same chromosome, inherited together) or trans (on opposite chromosomes) is critical for interpreting compound heterozygous effects in recessive disease.

### 9.5 Variant Annotation Results

Variant annotation via Ensembl VEP (GRCh37/38) classified all filtered variants into four impact tiers and identified specific codon and amino acid changes across coding genes on chromosomes 20 and 21.

**HIGH impact variants** include stop_gained (premature termination codons), frameshift insertions/deletions, and splice acceptor/donor variants. These are the most likely to disrupt protein function entirely. In tongue squamous cell carcinoma, stop-gained or frameshift mutations in established tumour suppressors — such as TP53 (chromosome 17), CDKN2A, FAT1, or NOTCH1 — are among the most frequently reported drivers. Any HIGH-impact variants identified on chromosomes 20 or 21 in genes with known tumour suppressor activity would be top-priority candidates for experimental validation.

**MODERATE impact variants (missense)** change a single amino acid in the protein sequence. These are evaluated further using:
- **SIFT:** A conservation-based score that asks whether the amino acid substitution is tolerated given the evolutionary conservation at that position. A SIFT score < 0.05 is classified as "deleterious."
- **PolyPhen-2:** A structure- and conservation-based score. A score > 0.908 is classified as "probably_damaging"; 0.447–0.908 as "possibly_damaging."

Variants rated "deleterious" by SIFT AND "probably_damaging" by PolyPhen-2 are highest priority among missense calls. In tongue cancer biology, missense mutations activating oncogenes (e.g., HRAS, PIK3CA) or inactivating tumour suppressors are well-established drivers.

**LOW impact variants (synonymous)** change the codon but not the encoded amino acid. These are generally treated as silent, though rare synonymous variants can affect mRNA splicing efficiency or translational speed in specific contexts.

**MODIFIER impact variants** fall outside protein-coding regions — in introns, UTRs, upstream/downstream regulatory elements, or intergenic space. While they rarely alter protein sequence directly, they may affect transcription factor binding, enhancer activity, splicing regulatory sequences, or non-coding RNA function.

### 9.6 Chromosome-Specific Gene Context

Chromosome 20 and 21 contain several genes with documented relevance in cancer:

**Chromosome 20** harbours ASXL1 (a chromatin modifier mutated in myeloid malignancies), PTPRT (a receptor tyrosine phosphatase with tumour suppressor function), TOP1 (topoisomerase I — a target of camptothecin-class drugs), and PCNA (proliferating cell nuclear antigen — a marker of cell proliferation). Missense or truncating variants in any of these in the context of tongue squamous cell carcinoma would warrant further investigation.

**Chromosome 21** carries ERG (an ETS family transcription factor implicated in prostate cancer and acute myeloid leukaemia), DYRK1A (a kinase regulating cell cycle and apoptosis), HMGN1 (a chromatin modifier with roles in DNA damage response), and RUNX1 (a transcription factor recurrently mutated in haematological malignancies). Missense variants in ERG or RUNX1 from this dataset would be notable findings.

### 9.7 Key Observations

- The GATK Best Practices pipeline (BQSR + duplicate marking + Mutect2 + FilterMutectCalls) produced a well-filtered variant call set with substantially reduced false-positive burden compared to naive variant calling.
- The parallel somatic and germline pipelines provide complementary views: somatic calling prioritises tumor-specific mutations, while germline calling with a trio context enables inheritance-based refinement and de novo variant detection.
- BQSR and duplicate marking are among the most impactful pre-processing steps — bypassing them would significantly inflate the false-positive rate, particularly for low-allele-frequency somatic mutations.
- Variant annotation via VEP transforms VCF coordinates into biologically interpretable predictions, enabling systematic prioritisation of variants by their predicted functional impact for experimental follow-up in tongue cancer studies.

---

## 10. Key Challenges in Variant Calling

### 1. Distinguishing True Variants from Noise

Every mismatch between a read and the reference is not necessarily a real variant. Sources of false positives include sequencing errors (incorrect base calls), PCR errors introduced during library preparation, misalignment to repetitive regions, and reference bias (reads containing non-reference alleles align less well than reference-matching reads). GATK's multi-step pipeline — BQSR, Mutect2, FilterMutectCalls — is designed to systematically address each source of error at the appropriate stage.

### 2. Reference Bias

Reads carrying the alternative allele may be penalised during alignment because the reference sequence always matches the reference allele. This can cause under-calling of heterozygous variants, especially at indel sites near repetitive sequences where alternative-allele reads are more likely to be soft-clipped or misplaced.

### 3. Allelic Drop-out and Amplification Bias

During PCR amplification, one allele may be preferentially amplified over the other, causing apparent imbalance in allele frequencies. This imbalance can be misinterpreted by Mutect2 as evidence for a low-frequency somatic mutation when it is in fact a technical artifact of uneven library amplification.

### 4. Data Size and Computational Resources

Full WGS samples are ~46 GB each. Processing the complete dataset requires approximately 100 GB RAM for GATK operations, 100+ GB disk space for all intermediate files, and multi-day runtimes on a single server. This project used a 1M-read subset and chromosomes 20/21 to make the analysis feasible within the internship timeframe. Scaling to the full genome would require a computing cluster or cloud infrastructure (e.g., AWS or Google Cloud with Cromwell/WDL workflow management).

---

## 11. Software Versions Reference

| Tool | Version | Purpose | Installation |
|---|---|---|---|
| SRA Toolkit (`fastq-dump`) | current | Download raw data from NCBI SRA | `wget` from NCBI FTP |
| FastQC | latest | Raw read quality assessment | `sudo apt-get install fastqc` |
| fastp | latest | Adapter trimming, quality filtering | `wget` from opengene.org |
| BWA (bwa-mem) | latest | Short read alignment to reference genome | `sudo apt-get install bwa` |
| SAMtools | latest | SAM/BAM conversion, indexing, flagstat | `sudo apt-get install samtools` |
| GATK | latest | Read groups, duplicate marking, BQSR, variant calling | Docker: `broadinstitute/gatk:latest` |
| bcftools | latest | VCF manipulation and filtering | `sudo apt install bcftools` |
| Ensembl VEP | GRCh37/38 | Variant effect prediction and annotation | Web interface |
| Docker | latest | Container runtime for GATK | System install |

### Key Reference Databases

| Database | File | Purpose |
|---|---|---|
| dbSNP 138 (hg38) | `Homo_sapiens_assembly38.dbsnp138.vcf` | Known variants for BQSR training |
| HapMap 3.3 (hg38) | `hapmap_3.3.hg38.vcf.gz` | Known germline variants for germline filtering |
| Reference Genome | `chr20.fa`, `chr21.fa` (UCSC hg38) | Alignment target |

---

*Pipeline script: `run_pipeline.sh`*  
*Supervisors: NyBerMan Bioinformatics Team / Dr. Sreeram Peela*
