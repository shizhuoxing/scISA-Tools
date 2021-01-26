# scIsoSeq-Analysis-Tool

**Authors: Zhuo-Xing Shi, SKLO, Sun Yat-Sen University; Ying-Dong He, BGI-Shenzhen**

*(If you are looking for a single cell matrix interactive transformation tool you can visit this link directly: https://github.com/shizhuoxing/scISA-Tools/wiki/scMatrixInteraction)*

Single-cell RNA sequencing is a powerful technique that advances gene expression regulation research to a higher resolution level for single cells. Recent advancements in scRNA-Seq methods (e.g., Droplet-based) allow thousands of cells to be captured and sequenced in a relatively short time frame at a fraction of the cost. They rely on capture of polyadenylated mRNA transcripts, followed by cDNA synthesis, pooling, amplification, library construction, and cDNA sequencing.

A major challenge for scRNA research is that the most widely used single cell RNA sequencing only generates short reads from one end of a cDNA template, thus limiting the reconstruction of highly diverse isoform such as alternative splicing and structural variation. Recent advances in long-read sequencing technologies present a potential solution to the shortcomings of short-read sequencing, as full-length cDNA reads can now encompass the entire sequence of transcripts. However, these long-read technologies typically suffer from higher error rates (Oxford Nanopore) and/or lower sequencing depth (PacBio) than short-read instruments.

Based on our previous development works, and by combining 10x genomics platform and HIT-scISOseq library construction technology, we have significantly improved the throughput of single-cell RNA sequencing using PacBio Sequel II platform. This has made it possible for us to obtain sufficient high-accuracy reads, and subsequently parse cell barcode, UMI and full-length transcript information in one dataset for thousands of single cells. It's really amazing!

In this repo, I will provide the tool set for analyzing HIT-scISOseq data starting from raw data. This is a comprehensive tool including basic data processing, cellBC extraction, expression matrix generation, and so on. It has been highly modularized so that you can combine various modules to address your need.

For citation and more detailed information, please see our preprint paper: [HIT-scISOseq: High-throughput and High-accuracy Single-cell Full-length Isoform Sequencing for Corneal Epithelium](https://www.biorxiv.org/content/10.1101/2020.07.27.222349v1.full)

# Dependencies
* SMRTlink 8.0 or later, you can install it in light way: `smrtlink_*.run --rootdir smrtlink --smrttools-only`
* ncbi-blast-2.10.0+ or later
* R-3.4.1 or later with `ggplot2 | gridExtra | grid | Seurat | DropletUtils`
* CellRanger 3.1.0 or later
* cDNA_Cupcake

# Usage
Export `smrtlink` `blast` `R` to your path first.
```
export PATH=$PATH:/smrtlink/smrtcmds/bin
export PATH=$PATH:/ncbi-blast-2.10.0+/bin
export PATH=$PATH:/R-3.4.1/bin
```

## Step1 Run CCS
```
ccs *.subreads.bam ccs.bam --min-passes 0 --min-length 50 --max-length 21000 --min-rq 0.75 -j 20
```
Start from SMRTlink8.0, CCS4.0 significantly speeds up the analysis and can be easily parallelized by using `--chunk`.

**Note:** the example_data folder in this repo have an `example.ccs.bam` for test use, which contain 5000 CCS reads produced from BGI patented HIT-scISOseq library construction protocol.

## Step2 Classify CCS by primer blast

### 2.1) convert ccs.bam to fasta and fastq
```
samtools view ccs.bam | awk '{print ">"$1"\n"$10}' > ccs.fa
samtools view ccs.bam | awk '{print "@"$1"\n"$10"\n+\n"$11}' > ccs.fq
```

### 2.2) make primer blast to CCS
```
makeblastdb -in primer.fa -dbtype nucl
blastn -query ccs.fa -db primer.fa -outfmt 7 -word_size 5 > mapped.m7
```

The following primer sequence is commonly used for PacBio original scIsoSeq library construction protocol with 10x genomics and BGI patented HIT-scISOseq library construction protocol with 10x genomics.

```
$ cat primer.fa
>primer_F
AAGCAGTGGTATCAACGCAGAGTACATGGG
>primer_S
CTACACGACGCTCTTCCGATCT
```

Have to note that, here primer_F means 5’ primer, primer_S means 3’ primer, you must set the name in this format for `classify_by_primer` utility.

### 2.3) classify CCS by primer

Here is an example for classifying CCS generate from PacBio original scIsoSeq library construction protocol with 10x genomics and BGI patented HIT-scISOseq library construction protocol with 10x genomics.

```
perl classify_by_primer.pl -blastm7 mapped.m7 -ccsfq ccs.fq -min_primerlen 16 -min_isolen 50 -outdir ./
```

`classify_by_primer` wraps a tool to detect full-length transcript from CCS base on PacBio original scIsoSeq library construction protocol with 10x genomics and BGI patented HIT-scISOseq library construction protocol with 10x genomics.

```
$ perl classify_by_primer.pl -h

Usage: perl classify_by_primer.pl -blastm7 mapped.m7 -ccsfq ccs.fq -min_primerlen 16 -min_isolen 50 -outdir ./

Options:
        -blastm7*:              result of primer blast to ccs.fa in blast -outfmt 7 format
        -ccsfq*:                the ccs.fq you want to classify to get full-length transcript
        -min_primerlen*:        the minimum primer alignment length in ccs.fa
        -min_isolen*:           the minimum output's transcript length whithout polyA tail
        -outdir*:               output directory
        -help:			print this help
```

After running finished, will find three file in your specified output path: `isoseq_flnc.BarcodeUMI.fastq`, `isoseq_flnc.Transcript.fastq`, `isoseq.PrimerStat.csv`

## Step3 Mapping FLNC reads to reference genome

After FLNC detection and primer, cellBC, UMI and polyA tail trimming, the remaining fraction of each isoforms FLNC reads were aligned to the reference genome with `minimap2`.

```
minimap2 -ax splice -t 6 -uf --secondary=no -C5 -t 5 ref.genome.fa isoseq_flnc.Transcript.fastq > isoseq_flnc.Transcript.sam
```

## Step4 Using gffcompare mapping FLNC to reference annotation

After FLNC mapping to reference genome, using `gffcompare`(https://ccb.jhu.edu/software/stringtie/gffcompare.shtml) to mapping FLNC to reference annotation.

### 4.1) transform FLNC mapped SAM file to GFF format
```
perl sam2gff.pl isoseq_flnc.Transcript.sam > isoseq_flnc.Transcript.gff
```

The script `sam2gff.pl` were original from: https://github.com/gpertea/gscripts/blob/master/sam2gff.pl

### 4.2) using gffcompare to mapping GFF to reference annotation.

```
gffcompare -r ref.genes.gtf isoseq_flnc.Transcript.gff -o flnc
```

### 4.3) filtering tmap for keep uniq mapping result only.

```
perl filter.tmap.pl flnc.isoseq_flnc.Transcript.gff.tmap > flnc.filtered.tmap
```

## Step5 Cell Barcode and UMI correction

We adopted a similar strategy of 10x Genomics `CellRanger` for cellBC and UMI correction, we have warps `CellRanger` cell barcode and UMI correction function as a module in our pipeline, named cellBC_UMI_corrector, these methods using long reads data independently, no need additional short read data for guiding.

```
/path/cellranger-3.1.0/miniconda-cr-cs/4.3.21-miniconda-cr-cs-c10/bin/python cellBC_UMI_corrector.py -w 3M-february-2018.txt -t 0.95 --input isoseq_flnc.BarcodeUMI.fastq --tmap flnc.filtered.tmap --output ./flnc.cellBC_UMI_correction.xls --cellranger_path /path/cellranger-3.1.0/
```

**Note:** you can find `3M-february-2018.txt.gz` file in example_data folder in this repo or in your cellranger path: `/path/cellranger-3.1.0/cellranger-cs/3.1.0/lib/python/cellranger/barcodes/3M-february-2018.txt.gz`

## Step6 Generation of Single Cell Gene Count matrix

Based on gffcompare output and corrected cellBC and UMI on each FLNC, we are using `scGene_matrix` utility to get the single cell gene expression quantity in each sample.

If NGS top rank cellBC list is provided, you can refer to the command line below:

```
perl scGene_matrix.pl -ngsbc NGScellBC.list -tgsbcumi flnc.cellBC_UMI_correction.xls -tmap flnc.filtered.tmap -outdir ./ -sample TEST
```

If you want using TGS top rank cellBC only, you can refer to the command line below:

```
perl scGene_matrix.pl -tgsbcumi flnc.cellBC_UMI_correction.xls -tmap flnc.filtered.tmap -topbc 2500 -outdir ./ -sample TEST
```

## Step7 Collapse Redundant Isoforms, Classification and Quality Control of Non-Redundant Isoforms

After mapping of FLNC to genome, we then using `cDNA_Cupcake`(https://github.com/Magdoll/cDNA_Cupcake) python script `collapse_isoforms_by_sam.py` to collapse all samples isoforms, and then we using `SQANTI3`(https://github.com/ConesaLab/SQANTI3) to characterization non-redundant isoforms.

### 7.1) example for collapsing.

```
sort -k 3,3 -k 4,4n isoseq_flnc.Transcript.sam > isoseq_flnc.Transcript.sorted.sam
/path/anaCogent/bin/python /path/anaCogent/bin/collapse_isoforms_by_sam.py --input isoseq_flnc.BarcodeUMI.fastq --fq -s isoseq_flnc.Transcript.sorted.sam -o TEST
```

### 7.2) example using SQANTI3.

```
python /path/SQANTI3-master/sqanti3_qc.py ./all.collapsed.rep.fa ref.genes.gtf ref.genome.fa -d ./TEST --isoAnnotLite -t 30
```

## Step8 Generation of Single Cell Isoform Count matrix

Based on collapsed output and corrected cellBC and UMI on each FLNC, we are using `scIsoform_matrix` utility to get the single cell isoform expression quantity in each sample.

If NGS top rank cellBC list provided, you can refer to the command line below:

```
perl scIsoform_matrix.pl -ngsbc NGScellBC.list -tgsbcumi flnc.cellBC_UMI_correction.xl -group all.collapsed.group.txt -minUMIcount 3 -outdir ./ -sample TEST
```

If you want using TGS top rank cellBC only, you can refer to the command line below:

```
perl scIsoform_matrix.pl -tgsbcumi flnc.cellBC_UMI_correction.xl -topbc 2500 -group all.collapsed.group.txt -minUMIcount 3 -outdir ./ -sample TEST
```

## Step9 Expression Matrix’s Quality Control

Each single cell gene and isoform expression matrix was using the `Seurat` R package (version 3.1.5) for further quality filtering analysis.

### 9.1) example using Seurat for expression matrix’s quality control and using DropletUtils for text matrix convert to CellRanger output format.

```
rm(list=ls())
library(dplyr)
library(Seurat)
library(Matrix)
library(DropletUtils)
metadata <- read.csv(TEST. isoform.matrix, row.names =1, sep=",")
test <- CreateSeuratObject(counts = metadata, min.cells = 5, min.features = 100, project = sample_name)
test[["percent.mt"]] <- PercentageFeatureSet(test, pattern = "^MT-")
test <- subset(test, subset = nFeature_RNA > 100 & nFeature_RNA < 6000 & percent.mt < 25)
write10xCounts(x =test@assays$RNA@counts, path = ./TEST)
```

### 9.2) example for fine-tuning.

```
cd ./TEST && awk '{print $1"\t"$2"\tGene Expression"}' genes.tsv > features.tsv && rm genes.tsv && sed -i 's/\./\-/' barcodes.tsv && gzip *
```

## Step10 Cell Clustering and Cell Type Annotation

After quality filtering of each sample’s single cell gene and isoform expression matrix, we using `scMatrix2CellRangerH5` utility to convert matrix to CellRanger h5 format, than we using CellRanger reanalyze pipeline for PCA and cell clustering with default parameters, the resulting cloupe file was loaded to `Loupe Browser` for fine manual annotation cell type and tune adjustments.

### 10.1) example for convert matrix to CellRanger h5 format.

```
/youpath/cellranger-3.1.0/miniconda-cr-cs/4.3.21-miniconda-cr-cs-c10/bin/python scMatrix2CellRangerH5.py –m ./TEST -c  example.filtered_feature_bc_matrix.h5 -o ./test.h5
```

**Note:** the example_data folder in this repo have an `example.filtered_feature_bc_matrix.h5`, you can copy and use this file for any `scMatrix2CellRangerH5` task.

### 10.2) example for CellRanger reanalyze.

```
cellranger reanalyze --id=TEST --matrix=./test.h5 --localcores=5
```

# Contact
If you have any questions, encounter problems or potential bugs, don’t hesitate to contact us. Report issues on github are welcome.
