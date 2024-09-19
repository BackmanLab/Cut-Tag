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

Script to call the EU-seq alignment script (Tophat) is below

```shell
#!/bin/bash
#SBATCH -A p32171 ##--## SLURM alt #SBATCH -A b1042
#SBATCH -p genhimem

#!/bin/bash
#SBATCH -A b1042
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=13G
#SBATCH --mail-user=lucascarter2025@u.northwestern.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH -t 48:00:00
#SBATCH --job-name=07.19.2024
#SBATCH --output=outlog
#SBATCH --error=errlog
##-- header here --##

echo -e "Usage: sh $0 -f <Forward Read> -w <Path to Working Directory> -g <Path to Bowtie Indices> -r <Path to Annotation File> -c <chromosome size file name > -t <threads> -m <min MAPQ score> -s <Spike in score (optional)> -p < path to SEACR (optional)> -b < sets MAC2 to --broad (optional)> -d <remove duplicates with Picard (optional)> \n"

## with MAC2
echo "Processing data"
for file in *R1_001.fastq.gz;
do echo "${file}";
   bash CutTag.sh -f ${file} -w /projects/b1042/BackmanLab/Lucas/090124_CnT/fastq/S3/Backman26_7.19.2024/ -g /projects/b1042/BackmanLab/Lucas/090124_CnT/genome/ref/hg38 -r /projects/b1042/BackmanLab/Lucas/090124_CnT/genome/gtf/hg38.ncbiRefSeq.gtf -c hg38.chrom.sizes -t 16 -m 2 -d ;
done

## with SEACR
echo "Processing data"
for file in *R1_001.fastq.gz;
do echo "${file}";
   bash CutTag.sh -f ${file} -w /projects/b1042/BackmanLab/Lucas/090124_CnT/fastq/test/ -g /projects/b1042/BackmanLab/Lucas/090124_CnT/genome/ref/hg38 -r /projects/b1042/BackmanLab/Lucas/090124_CnT/genome/gtf/hg38.ncbiRefSeq.gtf -p /home/lmc0633/executables/SEACR/SEACR_1.3.sh -c hg38.chrom.sizes -t 16 -m 2;
done

```

For the correlation plot generation script, i adapted this code for use as follows:

```R
##== R command ==##

library(corrplot)

list.files('/projects/b1042/BackmanLab/Lucas/090124_CnT/fastq/S1/BED/', pattern = ".bin.bed")

n=1 # c(1,2,3)
path <- '/projects/b1042/BackmanLab/Lucas/090124_CnT/fastq/'
subdir <- c("S1", "S2", "S3")
dir <- paste0(path, subdir[n], "/BED/")
hists <- list.files(dir, pattern = ".bin.bed")

reprod = c()
fragCount = NULL
for(hist in hists){
  
  if(is.null(fragCount)){
    
    fragCount = read.table(paste0(dir, hist), header = FALSE) 
    colnames(fragCount) = c("chrom", "bin", hist)
    
  }else{
    
    fragCountTmp = read.table(paste0(dir, hist), header = FALSE)
    colnames(fragCountTmp) = c("chrom", "bin", hist)
    fragCount = full_join(fragCount, fragCountTmp, by = c("chrom", "bin"))
    
  }
}

M = cor(fragCount %>% select(-c("chrom", "bin")) %>% log2(), use = "complete.obs") 

corrplot(M, method = "color", outline = T, addgrid.col = "darkgray", order="hclust", addrect = 3, rect.col = "black", rect.lwd = 3,cl.pos = "b", tl.col = "indianred4", tl.cex = 0.5, cl.cex = 0.5, addCoef.col = "black", number.digits = 2, number.cex = 0.5, col = colorRampPalette(c("midnightblue","white","darkred"))(100))


```
Script for calling just SEACR peaks. Can't use NORM setting unless I have an input control

```shell

#!/bin/bash
#SBATCH -A p32171 ##--## SLURM start
#SBATCH -p genhimem
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=200G
#SBATCH --mail-user=lucascarter2025@u.northwestern.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH -t 24:00:00
#SBATCH --job-name=build_genome
#SBATCH --output=outlog
#SBATCH --error=errlog
##-- header here --##

seacr=/home/lmc0633/executables/SEACR/SEACR_1.3.sh
bgfold=/projects/b1042/BackmanLab/Lucas/090124_CnT/fastq/S3/BEDGRAPH
pkfold=/projects/b1042/BackmanLab/Lucas/090124_CnT/fastq/S3/peaks

module load R/4.3.0
module load bedtools

bash /home/lmc0633/executables/SEACR/SEACR_1.3.sh /projects/b1042/BackmanLab/Lucas/090124_CnT/fastq/S3/BEDGRAPH/BM125_WT_S3.bedgraph 0.01 non stringent /projects/b1042/BackmanLab/Lucas/090124_CnT/fastq/S3/peaks/test.peaks
cd ${bgfold}

echo "Processing data"
for file in *.bedgraph;
do echo "${file}";


P1=$(echo $file | sed 's/.bedgraph/.seacr.peaks/g')

echo "8. Calling SEACR Peaks on ${file}"
echo -e "\t bash $seacr ${bgfold}/${file} 0.01 norm relaxed ${pkfold}/${P1}"

bash $seacr ${bgfold}/${file} 0.01 non stringent ${pkfold}/${P1}

done

echo -e "Peak calling complete \n"


```
#### Later modification

I looked into merging replicates using [Coverage()](https://ro-che.info/articles/2018-07-11-chip-seq-consensus) in R and unfortunately, you lose all the metadate (qValue, pValue, signal) when doing so. I decided to merge my technical replicates at the FASTQ level and then call the peaks again on those. Following calling peaks, I plan to use Bedtools intersect to find the conensus peaks since it preserves the metadata and other information. Merging all technical replicates in submission 2 and in submission 3

```shell

out=/projects/b1042/BackmanLab/Lucas/090124_CnT/fastq/merged/fastq/ 
in=/projects/b1042/BackmanLab/Lucas/090124_CnT/fastq/S3/Backman26_7.19.2024/

read1=BM135_Pol_2_degraded_S13_R2_001.fastq.gz  ; read2=BM138_Pol_2_degraded_S16_R2_001.fastq.gz ;  file=Merged.K9AB.P2KO.S2_R2_001.fastq.gz

#Merge input BAMS
cat ${in}${read1} ${in}${read2} > ${out}${file}

```

```shell
#!/bin/bash
#SBATCH -A b1042 ##--## SLURM start
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH --ntasks-per-node=10
#SBATCH --mem-per-cpu=10G
#SBATCH --mail-user=lucascarter2025@u.northwestern.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH -t 6:00:00
#SBATCH --job-name=compute matrix
#SBATCH --output=outlog
#SBATCH --error=errlog

module load deeptools

##== linux command ==##

inpath=/projects/b1042/BackmanLab/Lucas/090124_CnT/fastq/S2/peaks/
outpath=/projects/b1042/BackmanLab/Lucas/090124_CnT/fastq/S2/peaks/test/
bgpath=/projects/b1042/BackmanLab/Lucas/090124_CnT/fastq/S2/coverage/

bgfile=BM109_S1.rpgc.bigWig
pkfile=BM109_S1.seacr.peaks.relaxed.bed
otfile=$(echo $pkfile | sed 's/.seacr.peaks.relaxed.bed/.seacr.peaks.summitRegion.bed/g')
matfile=$(echo $pkfile | sed 's/.seacr.peaks.relaxed.bed/_SEACR.mat.gz/g')

label=P1AB_WT
rep=S2
cores=10

awk '{split($6, summit, ":"); split(summit[2], region, "-"); print summit[1]"\t"region[1]"\t"region[2]}' $inpath$pkfile > $outpath$otfile

computeMatrix reference-point -S $bgpath$bgfile \
              -R $outpath$otfile \
              --skipZeros -o $outpath$matfile -p $cores -a 3000 -b 3000 --referencePoint center

plotHeatmap -m $outpath$matfile -out ${outpath}${matfile}_SEACR_heatmap.png --sortUsing sum --startLabel "Peak Start" -\
-endLabel "Peak End" --xAxisLabel "" --regionsLabel "Peaks" --samplesLabel "${label} ${rep}"


inpath=/projects/b1042/BackmanLab/Lucas/090124_CnT/fastq/S2/peaks/
outpath=/projects/b1042/BackmanLab/Lucas/090124_CnT/fastq/S2/peaks/test/
bgpath=/projects/b1042/BackmanLab/Lucas/090124_CnT/fastq/S2/coverage/

bgfile=BM110_S2.rpgc.bigWig
pkfile=BM110_S2.seacr.peaks.relaxed.bed
otfile=$(echo $pkfile | sed 's/.seacr.peaks.relaxed.bed/.seacr.peaks.summitRegion.bed/g')
matfile=$(echo $pkfile | sed 's/.seacr.peaks.relaxed.bed/_SEACR.mat.gz/g')

label=P1AB_P1KO
rep=S2

awk '{split($6, summit, ":"); split(summit[2], region, "-"); print summit[1]"\t"region[1]"\t"region[2]}' $inpath$pkfile > $outpath$otfile

computeMatrix reference-point -S $bgpath$bgfile \
              -R $outpath$otfile \
              --skipZeros -o $outpath$matfile -p $cores -a 3000 -b 3000 --referencePoint center

plotHeatmap -m $outpath$matfile -out ${outpath}${matfile}_SEACR_heatmap.png --sortUsing sum --startLabel "Peak Start" -\
-endLabel "Peak End" --xAxisLabel "" --regionsLabel "Peaks" --samplesLabel "${label} ${rep}"

```

### downsampling reads 

I'm downsampling my reads because some are over-sequenced and some are undersequenced. I worry this will contribute to variability between conditions, since some samples did experience overamplification. For this task, I am using [seqtk](https://github.com/lh3/seqtk)

usage: seqtk sample -s <seed> <path/in/fastq.gz> <#ofsamples> > <path/out/fastq.gz>

First thing to do is count the total number of reads per file
```shell
#!/bin/bash
#SBATCH -A b1042 ##--## SLURM start
#SBATCH -p genomics
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=50G
#SBATCH --mail-user=lucascarter2025@u.northwestern.edu
#SBATCH -t 2:00:00
#SBATCH --job-name=count_reads

for file in *R1_001.fastq.gz;
do echo "${file}";
READS=$(expr $(zcat ${file} | wc -l) / 4)
echo -e "\t file: $file has $READS reads  "
echo -e "\t file: $file reads: $READS " >> reads.txt
done
```

S1:

|File                                |      Reads    	   | Name      |
|------------------------------------|-------------------|-----------|
| file: K4AB-Pol1_S4_R1_001.fastq.gz | reads: 50,054,872 | K4AB_P1KO |
| file: K4AB-WT_S3_R1_001.fastq.gz   | reads: 68,741,393 | K4AB_WT   |
| file: K9AB-Pol1_S8_R1_001.fastq.gz | reads: 5,858,924  | K9AB_P1KO |
| file: K9AB-WT_S7_R1_001.fastq.gz   | reads: 2,249,616  | K9AB_WT   |
| file: P1AB-Pol1_S2_R1_001.fastq.gz | reads: 21,542,758 | P1AB_P1KO |
| file: P1AB-WT_S1_R1_001.fastq.gz   | reads: 4,246,499  | P1AB_WT   |
| file: P2AB-Pol1_S6_R1_001.fastq.gz | reads: 56,878,563 | P2AB_P1KO |
| file: P2AB-WT_S5_R1_001.fastq.gz   | reads: 21,233,193 | P2AB_WT   |

S2:

|File                              |      Reads      | Name      |
|----------------------------------|-----------------|-----------|
|	 file: BM109_S1_R1_001.fastq.gz  |reads: 2,831,267 | P1AB_WT   |
|	 file: BM110_S2_R1_001.fastq.gz  |reads: 1,997,534 | P1AB_P1KO |
|	 file: BM111_S3_R1_001.fastq.gz  |reads: 2,761,436 | P1AB_P2KO |
|	 file: BM112_S4_R1_001.fastq.gz  |reads: 2,904,866 | P2AB_WT   |
|	 file: BM113_S5_R1_001.fastq.gz  |reads: 8,911,631 | P2AB_P1KO |
|	 file: BM114_S6_R1_001.fastq.gz  |reads: 3,923,244 | P2AB_P2KO |
|	 file: BM115_S7_R1_001.fastq.gz  |reads: 11,335,882| K27AB_P1KO|
|	 file: BM116_S8_R1_001.fastq.gz  |reads: 4,809,199 | K27AB_P2KO|
|	 file: BM117_S9_R1_001.fastq.gz  |reads: 408,680   | P1AB_WT   |
|	 file: BM118_S10_R1_001.fastq.gz |reads: 3,241,433 | P1AB_P1KO |
|	 file: BM119_S11_R1_001.fastq.gz |reads: 3,001,093 | P1AB_P2KO |
|	 file: BM120_S12_R1_001.fastq.gz |reads: 1,265,939 | P2AB_WT   |
|	 file: BM121_S13_R1_001.fastq.gz |reads: 4,573,860 | P2AB_P1KO |
|	 file: BM122_S14_R1_001.fastq.gz |reads: 6,784,824 | P2AB_P2KO |
|	 file: BM123_S15_R1_001.fastq.gz |reads: 2,668,701 | K27AB_P1KO|
|	 file: BM124_S16_R1_001.fastq.gz |reads: 2,413,773 | K27AB_P2KO|

Changing names to: 

BM109_S1_R1_001.fastq.gz P1AB_WT_S1_R1_001.fastq.gz
BM110_S2_R1_001.fastq.gz P1AB_P1KO_S2_R1_001.fastq.gz
BM111_S3_R1_001.fastq.gz P1AB_P2KO_S3_R1_001.fastq.gz
BM112_S4_R1_001.fastq.gz P2AB_WT_S4_R1_001.fastq.gz
BM113_S5_R1_001.fastq.gz P2AB_P1KO_S5_R1_001.fastq.gz
BM114_S6_R1_001.fastq.gz P2AB_P2KO_S6_R1_001.fastq.gz
BM115_S7_R1_001.fastq.gz K27AB_P1KO_S7_R1_001.fastq.gz
BM116_S8_R1_001.fastq.gz K27AB_P2KO_S8_R1_001.fastq.gz
BM117_S9_R1_001.fastq.gz P1AB_WT_S9_R1_001.fastq.gz
BM118_S10_R1_001.fastq.gz P1AB_P1KO_S10_R1_001.fastq.gz
BM119_S11_R1_001.fastq.gz P1AB_P2KO_S11_R1_001.fastq.gz
BM120_S12_R1_001.fastq.gz P2AB_WT_S12_R1_001.fastq.gz
BM121_S13_R1_001.fastq.gz P2AB_P1KO_S13_R1_001.fastq.gz
BM122_S14_R1_001.fastq.gz P2AB_P2KO_S14_R1_001.fastq.gz
BM123_S15_R1_001.fastq.gz K27AB_P1KO_S15_R1_001.fastq.gz
BM124_S16_R1_001.fastq.gz K27AB_P2KO_S16_R1_001.fastq.gz

S3:

|File                                                |  Reads           | Name      |
|----------------------------------------------------|------------------|-----------|
file: BM125_WT_S3_R1_001.fastq.gz                    | reads: 45,823,405  | P1AB_WT 
file: BM126_Pol_1_degraded_S4_R1_001.fastq.gz        | reads: 34,609,056  | P1AB_P1KO
file: BM127_Pol_2_degraded_S5_R1_001.fastq.gz        | reads: 28,356,451  | P1AB_P2KO
file: BM128_WT_S6_R1_001.fastq.gz                    | reads: 34,639,433  | CTCFAB_WT
file: BM129_WT_S7_R1_001.fastq.gz                    | reads: 49,503,419  | CTCFAB_WT
file: BM130_WT_S8_R1_001.fastq.gz                    | reads: 47,559,479  | CTCFAB_WT
file: BM131_Pol_1_degraded_S9_R1_001.fastq.gz        | reads: 38,030,748  | CTCFAB_P1KO
file: BM132_Pol_2_degraded_S10_R1_001.fastq.gz       | reads: 28,916,144  | CTCFAB_P2KO
file: BM133_WT_S11_R1_001.fastq.gz                   | reads: 63,835,066  | K9AB_WT
file: BM134_Pol_1_degraded_S12_R1_001.fastq.gz       | reads: 54,217,202  | K9AB_P1KO
file: BM135_Pol_2_degraded_S13_R1_001.fastq.gz       | reads: 53,024,658  | K9AB_P2KO
file: BM136_WT_S14_R1_001.fastq.gz                   | reads: 49,808,186  | K9AB_WT
file: BM137_Pol_1_degraded_S15_R1_001.fastq.gz       | reads: 53,189,723  | K9AB_P1KO
file: BM138_Pol_2_degraded_S16_R1_001.fastq.gz       | reads: 72,659,514  | K9AB_P2KO
file: BM139_Pol_1_degraded_S17_R1_001.fastq.gz       | reads: 21,882,985  | CTCFAB_P1KO
file: BM140_Pol_2_degraded_S18_R1_001.fastq.gz       | reads: 23,661,905  | CTCFAB_P2KO

Changing names to:

mv BM125_WT_S3_R2_001.fastq.gz P1AB_WT_S3_R2_001.fastq.gz
mv BM126_Pol_1_degraded_S4_R2_001.fastq.gz P1AB_P1KO_S4_R2_001.fastq.gz
mv BM127_Pol_2_degraded_S5_R2_001.fastq.gz P1AB_P2KO_S5_R2_001.fastq.gz
mv BM128_WT_S6_R2_001.fastq.gz CTCFAB_WT_S6_R2_001.fastq.gz
mv BM129_WT_S7_R2_001.fastq.gz CTCFAB_WT_S7_R2_001.fastq.gz
mv BM130_WT_S8_R2_001.fastq.gz CTCFAB_WT_S8_R2_001.fastq.gz
mv BM131_Pol_1_degraded_S9_R2_001.fastq.gz CTCFAB_P1KO_S9_R2_001.fastq.gz
mv BM132_Pol_2_degraded_S10_R2_001.fastq.gz CTCFAB_P2KO_S10_R2_001.fastq.gz
mv BM133_WT_S11_R2_001.fastq.gz K9AB_WT_S11_R2_001.fastq.gz
mv BM134_Pol_1_degraded_S12_R2_001.fastq.gz K9AB_P1KO_S12_R2_001.fastq.gz
mv BM135_Pol_2_degraded_S13_R2_001.fastq.gz K9AB_P2KO_S13_R2_001.fastq.gz
mv BM136_WT_S14_R2_001.fastq.gz K9AB_WT_S14_R2_001.fastq.gz
mv BM137_Pol_1_degraded_S15_R2_001.fastq.gz K9AB_P1KO_S15_R2_001.fastq.gz
mv BM138_Pol_2_degraded_S16_R2_001.fastq.gz K9AB_P2KO_S16_R2_001.fastq.gz
mv BM139_Pol_1_degraded_S17_R2_001.fastq.gz CTCFAB_P1KO_S17_R2_001.fastq.gz
mv BM140_Pol_2_degraded_S18_R2_001.fastq.gz CTCFAB_P2KO_S18_R2_001.fastq.gz

Then we can normalize based on read #

```shell

#!/bin/bash
#SBATCH -A p32171 ##--## SLURM start
#SBATCH -p short
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=150G
#SBATCH --mail-user=lucascarter2025@u.northwestern.edu
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH -t 3:00:00
#SBATCH --job-name=downsample
#SBATCH --output=outlog
#SBATCH --error=errlog

cutoff=1000000
norm=2000000
highnorm=10000000
seed=100

## Paths
inpath=$1 ## usage sbatch <sample.sh> </inpath/>
outpath=/projects/b1042/BackmanLab/Lucas/090124_CnT/fastq/sampled/fastq/

module load seqtk/Dec17

for file in *R1_001.fastq.gz;
do echo "${file}";

## File names
reads=$(expr $(zcat ${file} | wc -l) / 4)
R2=$(echo $file | sed 's/R1_001.fastq.gz/R2_001.fastq.gz/g')
R1out=$(echo $file | sed 's/R1_001.fastq.gz/sampled_R1_001.fastq/g')
R2out=$(echo $R2 | sed 's/R2_001.fastq.gz/sampled_R2_001.fastq/g')

if [[ $reads -lt $cutoff ]] ## If there is no path to seacr software, use MACS2
then

    echo -e "\t Sample ${file} has too few reads: ${reads} \n"

elif [[ $reads -lt $norm ]] && [[ $reads -gt $cutoff ]]
then

echo -e "\t Sample ${file} with ${reads} reads: sampled to 1,500,000 "
echo -e "\t seqtk sample -s $seed ${inpath}${file} 1500000 > ${outpath}${R1out} \n"

seqtk sample -s $seed ${inpath}${file} 1500000 > ${outpath}${R1out}
seqtk sample -s $seed ${inpath}${R2} 1500000 > ${outpath}${R2out}

elif [[ $reads -gt $highnorm ]] 
then

echo -e "\t Sample ${file} with ${reads} reads: sampled to $highnorm "
echo -e "\t seqtk sample -s $seed ${inpath}${file} $highnorm > ${outpath}${R1out} \n"

seqtk sample -s $seed ${inpath}${file} $highnorm > ${outpath}${R1out}
seqtk sample -s $seed ${inpath}${R2} $highnorm > ${outpath}${R2out}

else

echo -e "\t Sample ${file} with ${reads} reads: sampled to $norm "
echo -e "\t seqtk sample -s $seed ${inpath}${file} $norm > ${outpath}${R1out} \n"

seqtk sample -s $seed ${inpath}${file} $norm > ${outpath}${R1out}
seqtk sample -s $seed ${inpath}${R2} $norm > ${outpath}${R2out}

fi
done

gzip *.fastq


```

I checked the duplication rate using the following code and found that the rates for each sample were quite high (consistent with overamplication by PCR (my fault) and oversequencing by the core (NuSeq's fault)). [This publication](https://www.biorxiv.org/content/10.1101/2022.03.30.486382v1.full.pdf) addresses these issues and found similar duplication issues. Here is the code I used:

```R
## Summarize the duplication information from the picard summary outputs.

n=3
dir.sam <- "/projects/b1042/BackmanLab/Lucas/090124_CnT/fastq/"
subdirs <- c("S1", "S2", "S3")
sam.path <- paste0(dir.sam, subdirs[n],"/SAM/")
hists <- list.files(sam.path, pattern = ".dupMark.txt")

dupResult = c()
for(hist in hists){
  dupRes = read.table(paste0(sam.path,hist), header = TRUE, fill = TRUE)
  
  histInfo = strsplit(hist, "_")[[1]]
  dupResult = data.frame(Histone = histInfo[1], Replicate = histInfo[2], MappedFragNum_hg38 = dupRes$READ_PAIRS_EXAMINED[1] %>% as.character %>% as.numeric, DuplicationRate = dupRes$PERCENT_DUPLICATION[1] %>% as.character %>% as.numeric * 100, EstimatedLibrarySize = dupRes$ESTIMATED_LIBRARY_SIZE[1] %>% as.character %>% as.numeric) %>% mutate(UniqueFragNum = MappedFragNum_hg38 * (1-DuplicationRate/100))  %>% rbind(dupResult, .)
}

```
