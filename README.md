# Variant Calling and Annotation Pipeline for Tongue Cancer WGS Data

**NGS Internship Mini-Project**  
**BioProject:** [PRJEB62494](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJEB62494)

---

## Project Overview

This project implements a comprehensive variant calling and annotation pipeline applied to whole genome sequencing (WGS) data from tongue cancer samples and a cell line. The goal is to identify somatic and germline variants — including SNPs and indels — on chromosomes 20 and 21, and predict their functional effects using variant annotation tools.

Tongue cancer (oral squamous cell carcinoma of the tongue) is an aggressive malignancy where understanding the mutational landscape is critical for identifying potential driver mutations, therapeutic targets, and disease mechanisms. This pipeline follows GATK Best Practices for variant discovery.

### Dataset

| SRA Accession | Spots | Bases | Size | GC% | Published |
|---|---|---|---|---|---|
| ERR11468775 | 410.1M | 123.0G | 46.3 GB | 44.0% | 2024-02-18 |
| ERR11468776 | 408.4M | 122.5G | 46.2 GB | 43.8% | 2024-02-18 |
| ERR11468777 | 394.6M | 118.4G | 46.6 GB | 44.1% | 2024-02-18 |

- **Organism:** *Homo sapiens*
- **Sequencing platform:** Illumina HiSeq 1500
- **Data type:** Paired-end whole genome sequencing
- **Target chromosomes:** chr20 and chr21 (hg38/GRCh38)

---

## Tools Used

| Tool | Version | Purpose |
|---|---|---|
| SRA Toolkit (`fastq-dump`) | current | Download raw FASTQ from NCBI SRA |
| FastQC | latest | Raw read quality assessment |
| fastp | latest | Adapter trimming and quality filtering |
| BWA-MEM | latest | Read alignment to reference genome |
| SAMtools | latest | SAM/BAM conversion, sorting, indexing, flagstat |
| GATK | latest (Docker) | Read group tagging, duplicate marking, BQSR, variant calling |
| bcftools | latest | VCF filtering and variant type extraction |
| Ensembl VEP | GRCh37/38 | Variant effect prediction and annotation |

---

## Pipeline Workflow


<img width="1400" height="2284" alt="workflow (1)" src="https://github.com/user-attachments/assets/213fa789-10f0-47b2-9696-f47ceb7a5115" />



---

## Key Results Summary

- Variants (SNPs and indels) were successfully called on chromosomes 20 and 21 for all three samples.
- Somatic variants were called using GATK Mutect2 (tumor-only pipeline).
- Germline variants were called using GATK HaplotypeCaller with GVCF joint genotyping, demonstrating a trio-like setup (sample1: father, sample2: mother, sample3: child).
- SNPs and indels were extracted separately per sample for downstream filtering.
- Variants passing quality filter (QUAL > 50) were retained for annotation.
- Functional annotation via Ensembl VEP predicted effects including codon changes, amino acid changes, genomic region classification, and functional consequences (silent, missense).

---

## Repository Structure

```
ngs_internship/
├── ref_genome/
│   ├── chr20.fa            # Reference chromosome 20 (hg38)
│   ├── chr21.fa            # Reference chromosome 21 (hg38)
│   ├── chr20_ref.*         # BWA index files for chr20
│   └── chr21_ref.*         # BWA index files for chr21
├── raw_data/
│   ├── ERR11468775_1/2.fastq.gz      # Raw reads (1M reads subset)
│   ├── ERR11468776_1/2.fastq.gz
│   ├── ERR11468777_1/2.fastq.gz
│   └── *_trimmed_*.fastq.gz          # fastp-trimmed reads
├── bam_file/
│   ├── sample[1-6].sam               # Aligned reads (SAM)
│   ├── sample[1-6].bam               # Converted BAM files
│   ├── sample[1-6]_withRG.bam        # With read groups added
│   ├── sorted_sample[1-6].bam        # Coordinate-sorted BAM
│   ├── sample[1-6]_markedDups.bam    # Duplicate-marked BAM
│   ├── sample[1-6]_recal.bam         # BQSR-applied BAM
│   └── metrics_duplicates[1-6]       # Duplicate metrics
├── germline_call/
│   ├── sample[1-6]_germline.g.vcf.gz # Per-sample GVCFs
│   ├── test_db1/                      # GenomicsDB for chr20
│   ├── test_db2/                      # GenomicsDB for chr21
│   └── try_new.vcf.gz[1/2]           # Joint-genotyped VCF
├── *.vcf.gz                           # Mutect2 somatic VCF output
├── filtered_sample[1-6].vcf.gz       # FilterMutectCalls output
├── sample[1-6]_snps.vcf              # Extracted SNPs
├── sample[1-6]_indels.vcf            # Extracted indels
├── Homo_sapiens_assembly38.dbsnp138.vcf  # Known variants for BQSR
└── run_pipeline.sh                    # Complete reproducible pipeline
```

---

## How to Run

```bash
# Clone or copy the pipeline script
chmod +x run_pipeline.sh

# Run all parts in sequence (requires ~50GB disk space and 16+ CPU cores)
bash run_pipeline.sh

# Or run each part individually:
# Part 1: Data retrieval
# Part 2: QC and alignment
# Part 3: BAM preparation (mark duplicates, BQSR)
# Part 4: Variant calling (somatic + germline)
# Part 5: Variant filtering and extraction
# Part 6: Annotation via Ensembl VEP (web interface)
```

> **Note:** Parts 1–5 run on the command line. Part 6 (VEP annotation) uses the Ensembl web interface at https://grch37.ensembl.org/Multi/Tools/VEP

---

## Contributors

| Name | Role |
|---|---|
| Sanzida Akhter Anee | Pipeline design and implementation |

**Supervised by:** NyBerMan Bioinformatics Team / Dr. Sreeram Peela

---

## References

1. GATK Best Practices — https://gatk.broadinstitute.org
2. Ensembl VEP — https://www.ensembl.org/Tools/VEP
3. BioProject PRJEB62494 — https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJEB62494
4. BWA Manual — http://bio-bwa.sourceforge.net/bwa.shtml
5. SAMtools documentation — http://www.htslib.org
