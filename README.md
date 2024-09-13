## Cut&Tag Data Processing Guide

### Intro:

Cleavage Under Targets & Tagmentation (CUT&Tag)is an ChIP-like profiling strategy in which primary antibodies are bound to chromatin proteins in situ after nuclei have been permeabilized. Secondary antibodies covalently bound to the cut-and-paste transposase Tn5 are then bound to the primary antibody. Activation of the transposase simultaneously digests DNA and adds high throughput sequencing adapters (referred to as ‘tagmentation’) for paired-end DNA sequencing. Cut&Tag libraries have a much higher signal-to-noise ratio than traditional ChIP libraries and need far fewer cells as well.

Below is documentation for set-up and a single script for parsing Cut&Tag data based off of the data analysis protocol presented [here](https://yezhengstat.github.io/CUTTag_tutorial/#I_Introduction) from the Cut&Tag [protocol](https://www.nature.com/articles/s41596-020-0373-x) and [paper](https://www.nature.com/articles/s41467-019-09982-5). Additionally, there is an [nf-core pipeline](https://nf-co.re/cutandrun/3.2.2/) for Cut&Tag. Because I want  more control over the processing of the data, I will write my own wrapper based on the Henikoff Lab's data analysis protocol.

Here are related papers and resources for Cut&Tag:
- [A 2023 review on Cut&Tag](https://www.tandfonline.com/doi/full/10.1080/15592294.2023.2293411#abstract)
- [Integrative Analysis of CUT&Tag and RNA-Seq Data Through Bioinformatics: A Unified Workflow for Enhanced Insights](https://link.springer.com/protocol/10.1007/978-1-0716-4071-5_13)
- [CUT&RUNTools 2.0: a pipeline for single-cell and bulk-level CUT&RUN and CUT&Tag data analysis](https://academic.oup.com/bioinformatics/article/38/1/252/6318389)
- [Spatial-CUT&Tag: Spatially resolved chromatin modification profiling at the cellular level](https://www.science.org/doi/full/10.1126/science.abg7216?casa_token=xLnWrQG_IG0AAAAA%3AdHWxN4ylifpH74stpHLjaS4jEISVEni-8akaksmJ-kuijNiCDKb515VaHKQWkmHSk3O_kX4ZZTTh4I4)
- [Active Motif's Spike in Strategy](https://www.activemotif.com/documents/AACR-2024-Spike-In-Poster.pdf)
- [GoPeaks: histone modification peak calling for CUT&Tag](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-022-02707-w)
- [CUT&Tag recovers up to half of ENCODE ChIP-seq peaks in modifications of H3K27](https://www.biorxiv.org/content/10.1101/2022.03.30.486382v2.full)
- [ChIPseqSpikeInFree: a package for Spike-in Free analysis](https://github.com/stjude/ChIPseqSpikeInFree)
- [Peak calling by Sparse Enrichment Analysis for CUT&RUN chromatin profiling](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-019-0287-4)
- [Another Henikoff tutorial on Cut&Tag Analysis](https://www.protocols.io/view/cut-amp-tag-data-processing-and-analysis-tutorial-5jyl8py98g2w/v2) 
- [Really good resource on Bash operators](https://tldp.org/LDP/abs/html/comparison-ops.html)



I submitted three groups of Cut&Tag biological replicates on 6.03.2024,  07.01.2024, and 7.19.2024. The later two submissions have technical replicates in them that should be combined. The sample identifiers for each are below:

**6.03.2024 - Submission 1:**

|Cells          | AB	 | Index     | NuSeq ID |
|---------------|--------|-----------|----------|
| WT	        | Pol1AB | i71-i51	 | BM101    |
| Pol1 degraded | Pol1AB | i71-i52   | BM102    |
| Pol1 degraded | K4Me3AB| i71-i53	 | BM104    |
| WT	        |K4Me3AB |	i71-i54	 |BM103     |
| WT	        |Pol2AB	 |i72-i51	 |BM105     |
| Pol 1 degraded|Pol2AB	 |i72-i52	 |BM106     |
| WT	        |K9Me3AB |i72-i53	 |BM107     |
| Pol 1 degraded|K9Me3AB |i72-i54	 |BM108     |

**07.01.2024 - Submission 2:**

|Cells  Rep1    | AB	 | Index     | NuSeq ID |
|---------------|--------|-----------|----------|
|WT	            |Pol1AB  |i71-i51	|BM109 |
|Pol 1 degraded	|Pol1AB	 |i71-i52	|BM110 |
|Pol 2 degraded	|Pol1AB	 |i71-i53	|BM111 |
|WT	            |Pol2AB	 |i71-i54	|BM112 |
|Pol 1 degraded	|Pol2AB	 |i72-i51	|BM113 |
|Pol 2 degraded	|Pol2AB	 |i72-i52	|BM114 |
|Pol 1 degraded	|K27AcAB |i72-i53	|BM115 |
|Pol 2 degraded	|K27AcAB |i72-i54	|BM116 |

|Cells  Rep2    | AB	 | Index     | NuSeq ID |
|---------------|--------|-----------|----------|
|WT	| Pol1AB	| i73-i51	| BM117 |
|Pol 1 degraded	| Pol1AB |	i73-i52	| BM118 |
|Pol 2 degraded	| Pol1AB | i73-i53	| BM119 |
|WT	| Pol2AB	|	i73-i54	| BM120 | 
|Pol 1 degraded	| Pol2AB	| i74-i51	| BM121 |
|Pol 2 degraded	| Pol2AB	| i74-i52| BM122 |
|Pol 1 degraded	| K27AcAB	|	i74-i53	| BM123 |
|Pol 2 degraded	| K27AcAB	|	i74-i54	| BM124 |

**07.19.2024 - Submission 3:**

|Cells  Rep1    | AB	 | Index     | NuSeq ID |
|---------------|--------|-----------|----------|
|WT	|Pol1AB	|i71-i51	| BM125 |
|Pol 1 degraded	| Pol1AB	| i71-i52 | 	BM126 |
|Pol 2 degraded	| Pol1AB	| i71-i53 | BM127 | 
|WT	| CTCFAB	| i71-i54	| BM128 |
|WT | CTCFAB	| i72-i51	| BM129 |
|WT	| CTCFAB	| i72-i52	| BM130 |
|Pol 1 degraded	| CTCFAB	| i72-i53	| BM131 |
|Pol 2 degraded	| CTCFAB	| i72-i54	| BM132 |

|Cells  Rep2    | AB	 | Index     | NuSeq ID |
|---------------|--------|-----------|----------|
| WT	| K9AB	| i73-i51	| BM133 |
| Pol 1 degraded	| K9AB	| i73-i52	| BM134 |
| Pol 2 degraded	| K9AB	| i73-i53	| BM135 |
| WT	| K9AB	| i73-i54	| BM136 |
| Pol 1 degraded	| K9AB	| i74-i51	| BM137 |
| Pol 2 degraded	| K9AB	| i74-i52	| BM138 |
| Pol 1 degraded	| CTCFAB	| i74-i53	| BM139 |
| Pol 2 degraded	| CTCFAB	| i74-i54	| BM140 |



Cut&Tag uses [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml) for alignment.

Build Bowtie indices for reference genome prior to alignment with Tophat using the following command <kbd> bowtie2-build path_to_reference_genome.fa prefix_to_name_indexes </kbd> using default bowtie2 on quest, bowtie2/2.4.5

For mapping to rDNA repeats + human genome, use the [following reference](https://github.com/vikramparalkar/rDNA-Mapping-Genomes/blob/main/Human_hg38-rDNA_genome_v1.0_annotation.tar.gz). For general mapping, use Ensembl 

```shell
VERSION=108
wget -L ftp://ftp.ensembl.org/pub/release-$VERSION/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
wget -L ftp://ftp.ensembl.org/pub/release-$VERSION/gtf/homo_sapiens/Homo_sapiens.GRCh38.$VERSION.gtf.gz
```

OR UCSC: [UCSC Genome](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/) and [UCSC GTF](https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/) (preferred). There are 4 different GTF references to use.

```shell
wget -L https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.fa.gz
wget -L https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/genes/hg38.ncbiRefSeq.gtf.gz

```

example SLURM header for Quest HPC

```shell
#!/bin/bash
#SBATCH -A p32171 ##--## SLURM start
#SBATCH -p genhimem
#SBATCH -N 1
#SBATCH --ntasks-per-node=12
#SBATCH --mem-per-cpu=10G
#SBATCH --mail-user=lucascarter2025@u.northwestern.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH -t 24:00:00
#SBATCH --job-name=build_genome
#SBATCH --output=outlog
#SBATCH --error=errlog
```


```shell
##-- header here --##

cd /projects/b1042/BackmanLab/Lucas/090124_CnT/genome/ref

module load bowtie2/2.4.5

bowtie2-build /projects/b1042/BackmanLab/Lucas/090124_CnT/genome/ref/hg38.fa hg38
```

GTF file can be found at ``` /projects/b1042/BackmanLab/Lucas/090124_CnT/genome/gtf ``` and genome index can be found at ``` /projects/b1042/BackmanLab/Lucas/090124_CnT/genome/ref ```

### Spike in normalization

I am using [ChIPseqSpikeInFree](https://github.com/stjude/ChIPseqSpikeInFree) to normalize the data in lieu of a spike-in. I installed the command line version at:
```  /home/lmc0633/executables/ChIPseqSpikeInFree ``` using the command version of the installation detailed at the github and **R/4.3.0**

Need a metadata file with the following columns:

ANTIBODY	GROUPANTIBODY	GROUP
|ID | ANTIBODY | GROUP |
|---------------|--------|-----------|
|WT.bam	   |Pol1AB  | control	|
|Pol1.bam	|Pol1AB	 | pol1ko	|
|Pol2.bam |Pol1AB	 | pol2ko	|


Code is written in such a way that I can run all BAMs for a particular group at once. I deposited metadata at ``` /projects/b1042/BackmanLab/Lucas/090124_CnT/metadata ``` as individual .txt files for each submission. each submission is labelled S1, S2, or S3, with corresponding fastqs and alignment files in each respective directory. This can be run on SLURM or it can be run on an interaction session. Script should look something like below:

Shell submission script (if submitting via SLURM):

```shell

module load R/4.3.0

Rscript --vanilla --verbose spike.R

```

R script for Spike In normalization factor. Must have more than 1 core enabled for this to work:

```R

## HPC Script for normalization factor generation
library("ChIPseqSpikeInFree")

## Must have more than 1 core enabled for this to work
library(BiocParallel)
numCores <- parallel::detectCores()
register(MulticoreParam(workers = numCores-1), default = TRUE)

n=2 ## 1, 2, or 3
## Get metadata path (Must be tab delimited)
dir.meta <- "/projects/b1042/BackmanLab/Lucas/090124_CnT/metadata/"
dir.s <- c("06032024.txt", "07012024.txt", "07192024.txt")
meta <- paste0(dir.meta,dir.s[n])

## Get bam paths
dir.bam <- "/projects/b1042/BackmanLab/Lucas/090124_CnT/fastq/"
subdirs <- c("S1", "S2", "S3")
bam.path <- paste0(dir.bam, subdirs[n],"/BAM/")

## Get bam list
tbl <- read.table(meta, sep= "\t", header= T)
bams <- tbl$ID

## Write out txt files
#out <- c("06032024.txt", "07012024.txt", "07192024.txt")
#write.table(tbl, file = paste0(dir.meta, out[n]), quote=F, sep = "\t",row.names = F, col.names = TRUE)

library(stringr)
input <- str_trim(paste0(bam.path,bams), side = c("both", "left", "right"))

## Call function
ChIPseqSpikeInFree(bamFiles = input, chromFile = "hg38", metaFile = meta, prefix = paste0ls("/projects/b1042/BackmanLab/Lucas/090124_CnT/spikein/",subdirs[n]), ncores=4)

```

We need the size of each chromosome as a text file. this can be retrieved at UCSC:
```shell
wget -L https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/latest/hg38.chrom.sizes
```

Once we have the scaling factors from , we want to integrate them into generation of our coverage files like so:

```shell
libSize=`cat sample1.bed|wc -l`
scale=15000000/($libSize*$SF)
genomeCoverageBed -bg -scale $scale -i sample1.bed  -g hg38.chromSizes > sample1.bedGraph
bedGraphToBigWig sample1.bedGraph hg38.chromSizes sample1.bw
```
OR

```shell
scale_factor=`echo "10000 / $seqDepth" | bc -l`
    echo "Scaling factor for $histName is: $scale_factor!"
    bedtools genomecov -bg -scale $scale_factor -i $projPath/alignment/bed/${histName}_bowtie2.fragments.bed -g $chromSize > $projPath/alignment/bedgraph/${histName}_bowtie2.fragments.normalized.bedgraph

```

### Peak calling

There are a few options for calling peaks on Cut&Tag data. The more leniant and traditional option is [MACS2](https://github.com/macs3-project/MACS). Using MACS2 will increase the false positive rate but will capture more peaks overall. [SEACR](https://github.com/FredHutch/SEACR/tree/master) was designed with Cut&Run/TAG data (high SN) in mind. More recent and less used, [GoPeaks](https://github.com/maxsonBraunLab/gopeaks) is also specifically designed for Cut&Run/TAG data. There are several other peak callers, but for my purposes, I will use MACS2 and SEACR. [Here](https://help.pluto.bio/en/articles/choosing-between-macs2-and-seacr-for-peak-calling-with-dna-sequencing-datasets) is a short explanation on the difference between SEACR and MACS2. This [HBC workshop](https://github.com/hbctraining/Intro-to-ChIPseq-flipped/tree/main) covers calling peaks on Cut&Run/TAG data. 

For SEACR, we need to retrieve two scripts and put them in a directory for use later:

```shell
mkdir SEACR ## in bin or scripts
wget https://github.com/FredHutch/SEACR/raw/master/SEACR_1.3.sh
wget https://github.com/FredHutch/SEACR/raw/master/SEACR_1.3.R

```
SEACR usage and example: 
```shell 
## usage
bash SEACR_1.3.sh experimental bedgraph [control bedgraph | numeric threshold] ["norm" | "non"] ["relaxed" | "stringent"] output prefix 
## example
bash SEACR_1.3.sh target.bedgraph 0.01 non stringent output
```

For the script, we call the following:
```shell
seacr=/home/lmc0633/executables/SEACR/SEACR_1.3.sh
bash $seacr $projPath/alignment/bedgraph/${histName}_bowtie2.fragments.normalized.bedgraph 0.01 non stringent $projPath/peakCalling/SEACR/${histName}_seacr_top0.01.peaks


```


```shell
macs2 callpeak -t $projPath/alignment/bam/${org}/${time}/${hist}/bowtie2.mapped.bam \
      -c $projPath/alignment/bam/${org}/${timeControl}/${hist}/bowtie2.mapped.bam \
      -g ${genome} -f BAMPE -n macs2_peak_q0.1 --outdir $projPath/peakCalling/MACS2/${org}/${timeControl}_control/${time}/${hist}/
       -q 0.1 --keep-dup all 2>${projPath}/peakCalling/MACS2/${org}/${timeControl}_control/${time}/${hist}/macs2Peak_summary.txt

```

```shell
mkdir -p $projPath/peakCalling
macs2 callpeak -t ${projPath}/alignment/bam/${histName}_rep1_bowtie2.mapped.bam \
      -c ${projPath}/alignment/bam/${controlName}_rep1_bowtie2.mapped.bam \
      -g hs -f BAMPE -n macs2_peak_q0.1 --outdir $projPath/peakCalling/MACS2 -q 0.1 --keep-dup all 2>${projPath}/peakCalling/MACS2/macs2Peak_summary.txt

```
double array in linux

```shell
#!/bin/bash

array=( "Vietnam" "Germany" "Argentina" )
array2=( "Asia" "Europe" "America" )

for i in "${!array[@]}"; do
    printf "%s is in %s\n" "${array[i]}" "${array2[i]}"
done

```

Script to call the CutTag alignment script  is below

```shell
##-- header here --##

echo -e "Usage: sh $0 -f <Forward Read> -w <Path to Working Directory> -g <Path to Bowtie Indices> -r <Path to Annotation File> -t <threads> -m <min MAPQ score> -s <Spike in score (optional)> -p < path to SEACR (optional)> -b < sets MAC2 to --broad (optional)> \n"

## with MAC2
echo "Processing data"
for file in *R1_001.fastq.gz;
do echo "${file}";
   bash CutTag.sh -f ${file} -w /projects/b1042/BackmanLab/Lucas/090124_CnT/fastq/S3/Backman26_7.19.2024/ -g /projects/b1042/BackmanLab/Lucas/090124_CnT/genome/ref/hg38 -r /projects/b1042/BackmanLab/Lucas/090124_CnT/genome/gtf/hg38.ncbiRefSeq.gtf -c hg38.chrom.sizes -t 24 -m 2;
done

## with SEACR
echo "Processing data"
for file in *R1_001.fastq.gz;
do echo "${file}";
   bash CutTag.sh -f ${file} -w /projects/b1042/BackmanLab/Lucas/090124_CnT/fastq/test/ -g /projects/b1042/BackmanLab/Lucas/090124_CnT/genome/ref/hg38 -r /projects/b1042/BackmanLab/Lucas/090124_CnT/genome/gtf/hg38.ncbiRefSeq.gtf -p /home/lmc0633/executables/SEACR/SEACR_1.3.sh -c hg38.chrom.sizes -t 16 -m 2;
done

```
