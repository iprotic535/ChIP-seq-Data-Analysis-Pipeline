# ChIP-seq Data Analysis Pipeline

This repository documents the complete workflow for analyzing H3K27ac ChIP-seq data using a combination of Bash and R scripts. The dataset involves honeybee genome with a genome size of ~225Mb and includes preprocessing, alignment, peak calling, visualization, and annotation.

---

## 🧪 Experimental Setup

- **Target histone mark**: H3K27ac
- **Control**: Input DNA (matched sample)
- **Organism**: Apis mellifera (Honeybee)
- **Genome Size (effective)**: 2.25e8
- **Sequencing**: Paired-end
- **Environment**: Conda virtual environment (`ChIPseq`, `chipseqTools`)
- **Tools**: `fastp`, `bwa`, `samtools`, `macs3`, `deepTools`, `ChIPseeker`, `TxDb`, `ggplot2`

---

## 📁 Directory Structure

```
ChIPseq_project/
├── input/                    # Raw and trimmed FASTQ files
├── aligned/                 # BAM and sorted BAM files
├── peaks/                   # MACS3 peak output
├── bigwig/                  # bigWig files for visualization
├── annotation/              # Annotated peaks
└── scripts/                 # Bash and R scripts
```

---

## 🧬 Pipeline Steps

### 1. Preprocessing with `fastp`
```bash
fastp   -i input/SRR28865860_1.fastq   -I input/SRR28865860_2.fastq   -o cleaned/SRR28865860_1.trim.fastq   -O cleaned/SRR28865860_2.trim.fastq   --detect_adapter_for_pe   --thread 8   --html SRR28865860_fastp.html   --json SRR28865860_fastp.json
```

### 2. Alignment with BWA and BAM cleanup
```bash
# Index genome
bwa index genome.fa

# Align
bwa mem -t 8 genome.fa SRR28865860_1.trim.fastq SRR28865860_2.trim.fastq > SRR28865860.sam

# Convert SAM to BAM
samtools view -Sb SRR28865860.sam | samtools sort -o SRR28865860.sorted.bam

# Add mate info and mark duplicates
samtools fixmate -m SRR28865860.sorted.bam SRR28865860.fixmate.bam
samtools sort -o SRR28865860.positionsorted.bam SRR28865860.fixmate.bam
samtools markdup -r SRR28865860.positionsorted.bam SRR28865860.dedup.bam
samtools index SRR28865860.dedup.bam
```

### 3. Peak Calling with `MACS3`
```bash
macs3 callpeak   -t aligned/SRR28865737.dedup.bam   -c input/aligned/SRR28865860.dedup.bam   -f BAMPE   --broad   --broad-cutoff 0.1   -g 2.25e8   -n SRR28865737_H3K27ac   --outdir peaks/
```

### 4. Signal Track Generation (bigWig)
```bash
bamCoverage   -b aligned/SRR28865737.dedup.bam   -o bigwig/SRR28865737_H3K27ac.bw   --binSize 50   --normalizeUsing RPKM   -p 8 

#RPKM = Reads Per Kilobase per Million mapped reads (adapted for bins instead of genes):
#Adjusts for library size (so samples with more total reads don’t look artificially higher).
#Adjusts for region length (in this case, bin length).
```

### 5. Visualization with IGV
- Load the `bigwig/SRR28865737_H3K27ac.bw` file in IGV along with `genome.fa` and `genome.fa.fai`.

---

## 📊 Peak Annotation (R Script)

```r
library(ChIPseeker)
library(GenomicFeatures)

txdb <- makeTxDbFromGFF("Amel_4.5.gff3", format="gff3")
peakfile <- "peaks/SRR28865737_H3K27ac_peaks.broadPeak"

peak_anno <- annotatePeak(
  peakfile,
  TxDb = txdb,
  tssRegion = c(-1000, 200),
  annoDb = NA
)

plotAnnoBar(peak_anno)
plotDistToTSS(peak_anno)
```

---

## 📑 Methods (for Manuscript)

Raw H3K27ac ChIP-seq paired-end FASTQ files were preprocessed using `fastp` to remove adapters and low-quality reads. The reads were aligned to the Apis mellifera reference genome using `bwa`, followed by sorting, duplicate removal, and indexing using `samtools`. Peak calling was performed using `MACS3` with the `--broad` option for H3K27ac histone marks, using a genome size of 2.25e8. Signal tracks were generated using `deepTools` and visualized in IGV. Peaks were annotated using `ChIPseeker` in R with a custom `TxDb` object generated from the honeybee GFF3 file. Figures were generated using `ggplot2`.

---

## 📂 Scripts Directory

- `scripts/ChIPseq_pipeline.sh` - Full bash pipeline
- `scripts/peak_annotation.R` - R script for peak annotation and plots

---

## 🧵 Versioning and Reproducibility

- Python: 3.11
- MACS3: v3.0.3
- deepTools: latest via conda
- R: v4.3.0
- ChIPseeker: Bioconductor v1.34+

---

## ✍️ Author

Ismam Ahmed Protic  
PhD Student, University of Nevada, Reno

---

## 📜 License

This project is for academic research use only. Contact the author for permissions.
