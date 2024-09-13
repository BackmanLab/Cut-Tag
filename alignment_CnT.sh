#!/bin/bash
##################################################################
# Author: Lucas Carter                                                   
# Email: lucascarter2025@u.northwestern.edu                                                          
# PI: Vadim Backman                                                   
# Description: 
# This script performs QC, trimming, mapping, coverage, and peak 
# calling forCut&Tag generated data. It generates directories, 
# loops through fastqs, and moves results into respective directories.
# This script is written to run on a SLURM HPC. A more thorough 
# description can be found at README.md
################################################################

##--------------------------------------------------------------## Module load

## Customize these for your environment
module load samtools/1.6
module load bowtie2/2.5.4
module load fastqc
module load TrimGalore/0.6.10
module load java/jdk1.8.0_191
module load picard/2.21.4
module load bedtools/2.30.0 
module load deeptools/3.1.1

##--------------------------------------------------------------## Script begin

### Usage function tells users how to run the software
helpFunction()
{
      echo "*********************************** how to use CutTag.sh ***********************************"
      echo -e "Usage: sh $0 -f <Forward Read> -w <Path to Working Directory> -g <Path to Bowtie Indices> -r <Path to Annotation File> -c <chromosome size file name > -t <threads> -m <min MAPQ score> -s <Spike in score (optional)> -p < path to SEACR (optional)> -b < sets MAC2 to --broad (optional)> \n"
      echo -e "Run script in directory where fastqs are located \n"
      echo "This script completes several functions: Read QC, Quality Control, Alignment, Sorting, and read count of features:"
      echo -e "\t -- Quality Control: eliminates the adaptor and low quality reads using trim_galore."
      echo -e "\t -- Mapping: Mapping via Bowtie2/Tophat for all features"
      echo -e "\t -- Sorting: Sort the bam files, mark uniquely mapping reads, generate coverage files "
      echo -e "\t -- Conversion: Convert clean BAMs to BED and BEDGRAPH format "
      echo -e "\t -- Peak Calling: Call peaks using MACS2 or SEACR"
      echo -e "\t -- Read QC: Checks quality of reads using Fastqc \n \n"
      echo -e "\t-h help \n"

      echo "For more detail information, please feel free to contact: lucascarter2025@u.northwestern.edu"
      echo "**************************************************"
   exit 1 # Exit script after printing help
}


##--------------------------------------------------------------## options

while getopts "f:w:g:r:c:t:m:s:p:b" opt
do
   case $opt in
      f) R1=$OPTARG ;; ## forward read | string
      w) wdir=$OPTARG ;; ## working directory | string
      g) gendir=$OPTARG ;; ## path/to/bowtie_indices | string
      r) gtfdir=$OPTARG ;; ## path/to/gtf_file | string
      c) csize=$OPTARG ;; ## name of chrom sizes file | string
      t) threads=$OPTARG ;; ## number of threads to HPC | integer
      m) minScore=$OPTARG ;; ## minimum MAPQ score | integer
      s) spikein=$OPTARG ;; ## spike in value | integer
      p) seacr=$OPTARG ;; ## call peaks with SEACR instead of MACS2 | path to SEACR.sh/.r | string
      b) broad=1 ;; ## broad peaks option for MACS2 | no input
      ?) echo "Unknown argument --${OPTARG}"; helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "${R1}" ] 
then
   echo "*** error: must specify forward read ***";
   helpFunction
fi

if [ -z "${wdir}" ] 
then
   echo "*** error: working directory where FASTQs are located must be provided ***";
   helpFunction
fi

if [ -z "${gendir}" ] 
then
   echo "*** error: path to reference genome indices must be provided ***";
   helpFunction
fi

if [ -z "${gtfdir}" ] 
then
   echo "*** error: path to genome annotation file must be provided ***";
   helpFunction
fi

if [ -z "${csize}" ] 
then
   echo "*** error: must provide chromosome size file ***";
   helpFunction
fi

if [ -z "${threads}" ] 
then
   echo "*** error: number of processors must be provided ***";
   helpFunction
fi

if [ -z "${minScore}" ] 
then
   echo "*** error: a minimum MAPQ score must be provided ***";
   helpFunction
fi

##--------------------------------------------------------------## Paths

# Get absolute path for scripts and results; check if required scripts exist

# Set working directory by moving script to directory where fastqs are located
cd ${wdir}

echo "$(date): Processing alignment of FASTQ sequence files"
echo -e "Generating directories. \n"

# make sure references are prepared correctly
root="$(dirname "$(pwd)")/"

echo "Workign directory is ${root}"
echo "Genome annotation file is located at ${gtfdir}"
echo -e "BOWTIE2 genome assembly is located in ${gendir} \n"

# Make .OUT directory for trimmed reads
outfold="${root}OUT"
[ ! -d $outfold ]&&mkdir $outfold

# make .SAM directory
samfold="${root}SAM"
[ ! -d $samfold ]&&mkdir $samfold

# Make QC directory
qcfold="${root}fastqc"
[ ! -d $qcfold ]&&mkdir $qcfold

# Make coverage directory
covfold="${root}coverage"
[ ! -d $covfold ]&&mkdir $covfold

# Make .BAM directory
bamfold="${root}BAM"
[ ! -d $bamfold ]&&mkdir $bamfold

# Make .BED directory
bedfold="${root}BED"
[ ! -d $bedfold ]&&mkdir $bedfold

# Make .BEDGRAPH directory
bgfold="${root}BEDGRAPH"
[ ! -d $bgfold ]&&mkdir $bgfold

# Make peaks directory
pkfold="${root}peaks"
[ ! -d $pkfold ]&&mkdir $pkfold

##--------------------------------------------------------------## 1. Trim and QC
# begin data processing

cd ${wdir} ## move to fastq directory
R2=$(echo $R1 | sed 's/R1/R2/g')

###
echo "1. Starting the quality control: trimming ${R1} and ${R2}"
echo " **** 1.2 Running Trim Galore on all FASTA files now: "
echo -e " \t trim_galore --paired --retain_unpaired  --dont_gzip -o $outfold ${root}${PWD##*/}/${R1} ${root}${PWD##*/}/${R2}"

trim_galore --paired --retain_unpaired --dont_gzip -o $outfold ${root}${PWD##*/}/${R1} ${root}${PWD##*/}/${R2}

echo -e "${R1} and ${R2} have been processed \n"
###

##--------------------------------------------------------------## 2. Mapping

cd ${outfold} ## move to trimmed fastq folder
T1=$(echo $R1 | sed 's/R1_001.fastq.gz/R1_001_val_1.fq/g')
T2=$(echo $R2 | sed 's/R2_001.fastq.gz/R2_001_val_2.fq/g')
S1=$(echo $T1 | sed 's/R1_001_val_1.fq/R1_001.sam/g')

###
echo "2. Starting the mapping for ${T1} and ${T2}"
echo -e " \t bowtie2 --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p ${threads} -x ${gendir} -1 ${outfold}/${T1} -2 ${outfold}/${T2} -S ${samfold}/${S1} &> ${samfold}/${S1}_bowtie2.txt"
### Mapp reads using Bowtie
bowtie2 --local --very-sensitive --no-mixed --no-discordant --phred33 -I 10 -X 700 -p ${threads} -x ${gendir} -1 ${outfold}/${T1} -2 ${outfold}/${T2} -S ${samfold}/${S1} &> ${samfold}/${S1}_bowtie2.txt

echo -e "Mapping via Bowtie2 complete for ${T1} and ${T2} \n"
###

##--------------------------------------------------------------## 3. Mark duplicates

cd ${samfold} ## move to SAM alignment folder
S2=$(echo $S1 | sed 's/R1_001.sam/R1_001.sort.sam/g')
S3=$(echo $S2 | sed 's/R1_001.sort.sam/R1_001.sort.md.sam/g')

###
echo "3. Sorting ${S1} and marking duplicates in ${S2}"
echo -e " \t picard SortSam I=${samfold}/${S1} O=${samfold}/${S2} SORT_ORDER=coordinate"
echo -e " \t picard MarkDuplicates I=${samfold}/${S2} O=${samfold}/${S3} METRICS_FILE=${samfold}/${R1}_dupMark.txt"

## Sort SAM files
picard SortSam I=${samfold}/${S1} O=${samfold}/${S2} SORT_ORDER=coordinate

## Mark duplicates in sorted SAM files
picard MarkDuplicates I=${samfold}/${S2} O=${samfold}/${S3} METRICS_FILE=${samfold}/${R1}_dupMark.txt

echo -e "Duplicates marked in ${S3} \n"
###

##--------------------------------------------------------------## 4. Sort and statistics

cd ${samfold} ## move to SAM alignment folder
B1=$(echo $S3 | sed 's/_R1_001.sort.md.sam/.bam/g')
B2=$(echo $B1 | sed 's/.bam/.ind.bam/g')

###
echo "4. sorting Bowtie2 mapped reads. Retaining the uniquely mapping reads for ${S3}. Generating aligmnent statistics."

samtools view -bS -F 0x04 -q ${minScore} ${samfold}/${S3} -o ${bamfold}/${B1} ## filter reads by MAPq
samtools sort ${bamfold}/${B1} > ${bamfold}/${B2}
samtools index ${bamfold}/${B2}  ## index reads
samtools flagstat ${bamfold}/${B2} > ${bamfold}/${B2}.FLAGSTAT.txt ## get alignment scores

echo -e "samtools finished processing ${S3} to ${B1} \n"
###

##--------------------------------------------------------------## 5. unscaled BAM coverage

cd ${bamfold}
BW1=$(echo $B2 | sed 's/.ind.bam/.cpm.bigWig/g')
BW2=$(echo $B2 | sed 's/.ind.bam/.rpgc.bigWig/g')

###
echo "4. Generating coverage files for ${B2}"
echo -e "\t bamCoverage --bam ${bamfold}/${B2} --normalizeUsing CPM --outFileName ${covfold}/${BW1} --binSize 1 --numberOfProcessors max"
echo -e "\t bamCoverage --bam ${bamfold}/${B2} --normalizeUsing  RPGC --outFileName ${covfold}/${BW2} --effectiveGenomeSize 2747877702 --binSize 1 --numberOfProcessors max"

# Bam coverage to generate BigWig
bamCoverage --bam ${bamfold}/${B2} --normalizeUsing CPM --outFileName ${covfold}/${BW1} --binSize 1 --numberOfProcessors max
bamCoverage --bam ${bamfold}/${B2} --normalizeUsing  RPGC --outFileName ${covfold}/${BW2} --effectiveGenomeSize 2747877702 --binSize 1 --numberOfProcessors max

echo -e "BAM coverage files generated for ${B2} using CPM and RPGC normalization \n"
### 

##--------------------------------------------------------------## 6. Bedtools conversion

cd ${bamfold} ## move to BAM folder
BD1=$(echo $B2 | sed 's/.ind.bam/.bed/g')
BD2=$(echo $BD1 | sed 's/.bed/.filt.bed/g')
BD3=$(echo $BD2 | sed 's/.filt.bed/.fragments.bed/g')
BD4=$(echo $BD3 | sed 's/.fragments.bed/.bin.bed/g')


###
echo "6. Converting ${B2} to BED files ${BD1}."
echo -e "\t bedtools bamtobed -i ${bamfold}/${B2} -bedpe >${bedfold}/${BD1}"
echo -e "\t awk '$1==$4 && $6-$2 < 1000 {print $0}' ${bedfold}/${BD1} >${bedfold}/${BD2} "
echo -e "\t cut -f 1,2,6 ${bedfold}/${BD2}  | sort -k1,1 -k2,2n -k3,3n  >${bedfold}/${BD3}"
echo -e "\t awk -v w=$binLen '{print $1, int(($2 + $3)/(2*w))*w + w/2}' ${bedfold}/${BD3}  | sort -k1,1V -k2,2n | uniq -c | awk -v OFS="\t" '{print $2, $3, $1}' |  sort -k1,1V -k2,2n  >${bedfold}/${BD4}"

bedtools bamtobed -i ${bamfold}/${B2} -bedpe >${bedfold}/${BD1}

## Keep the read pairs that are on the same chromosome and fragment length less than 1000bp.
awk '$1==$4 && $6-$2 < 1000 {print $0}' ${bedfold}/${BD1} >${bedfold}/${BD2} 

## Only extract the fragment related columns
cut -f 1,2,6 ${bedfold}/${BD2}  | sort -k1,1 -k2,2n -k3,3n  >${bedfold}/${BD3}

## Bin genome for correlation calculation later
binLen=500
awk -v w=$binLen '{print $1, int(($2 + $3)/(2*w))*w + w/2}' ${bedfold}/${BD3}  | sort -k1,1V -k2,2n | uniq -c | awk -v OFS="\t" '{print $2, $3, $1}' |  sort -k1,1V -k2,2n  >${bedfold}/${BD4}

echo -e "Done converting ${B2} to ${BD2} and ${BD3}. Generated binned genome ${BD4} for correlation analysis. \n"
###

##--------------------------------------------------------------## 7. Bedgraph generation

gpath=$(echo $gendir | rev | cut -d'/' -f2- | rev) ## remove last part of directory
chromsize=${gpath}/${csize}
cd ${bedfold}

###
echo "7. Converting ${BD3} to BEDGRAPH for peak calling."
echo -e "\t chromosome size path is ${chromsize}"

if [[ -z "${spikein}" ]] ## If there is no spikein generate bedgraph without
   then

   BG1=$(echo $BD3 | sed 's/.fragments.bed/.bedgraph/g')
	echo -e "\t Generating BEDGRAPH without Spike-In"
   echo -e "\t bedtools genomecov -bg -i ${bedfold}/${BD3}  -g ${chromsize} > ${bgfold}/${BG1}"

   bedtools genomecov -bg -i ${bedfold}/${BD3}  -g ${chromsize} > ${bgfold}/${BG1}
   
   echo -e "Done converting ${BD3} to ${BG1}. \n"
   ###  

   else

   BG1=$(echo $BD3 | sed 's/.fragments.bed/.bedgraph/g')
   BW3=$(echo $BD3 | sed 's/.fragments.bed/.scaled.bigWig/g')

   libSize=`cat ${bedfold}/${BD3}|wc -l`
   scale=`echo "15000000/(${libSize}*${spikein})" | bc -l`
	echo -e "\t Generating BEDGRAPH and BIGWIG with scaled value from Spike-In: ${scale}"
   echo -e "\t bedtools genomecov -bg -scale ${scale} -i ${bedfold}/${BD3}  -g ${chromsize} > ${bgfold}/${BG1}"
   echo -e "\t bamCoverage --bam ${bamfold}/${B2} --scaleFactor ${scale} --normalizeUsing  RPGC --outFileName ${covfold}/${BW3} --effectiveGenomeSize 2747877702 --binSize 1 --numberOfProcessors max"

   bedtools genomecov -bg -scale ${scale} -i ${bedfold}/${BD3}  -g ${chromsize} > ${bgfold}/${BG1}
   bamCoverage --bam ${bamfold}/${B2} --scaleFactor ${scale} --normalizeUsing  RPGC --outFileName ${covfold}/${BW3} --effectiveGenomeSize 2747877702 --binSize 1 --numberOfProcessors max
   
   echo -e "Done converting ${BD3} to ${BG1} and ${BW3}. \n"
   ###   
   
fi

##--------------------------------------------------------------## 8. Peak calling

###

if [ -z "${seacr}" ] ## If there is no path to seacr software, use MACS2
then

cd ${bamfold}
P1=$(echo $BG1 | sed 's/.bedgraph/.macs2/g')
module load MACS2/2.2.9.1

echo "8. Calling MAC2 Peaks on ${B1}"
echo -e "\t macs2 callpeak -t ${bamfold}/${B1} -g hs -f BAMPE -n ${P1} --outdir ${pkfold} -q 0.1 --keep-dup all 2>${pkfold}/${P1}.macs2Peak_summary.txt"

macs2 callpeak -t ${bamfold}/${B1} -g hs -f BAMPE -n ${P1} --outdir ${pkfold} -q 0.1 --keep-dup all 2>${pkfold}/${P1}.macs2Peak_summary.txt

elif [[ $broad -eq 1 ]]
then

cd ${bamfold}
P1=$(echo $BG1 | sed 's/.bedgraph/.macs2/g')
module load MACS2/2.2.9.1

echo "8. Calling MAC2 Peaks on ${B1} with --broad option"
echo -e "\t macs2 callpeak -t ${bamfold}/${B1} -g hs -f BAMPE --broad -n ${P1} --outdir ${pkfold} -q 0.1 --keep-dup all 2>${pkfold}/${P1}.macs2Peak_summary.txt"

macs2 callpeak -t ${bamfold}/${B1} -g hs -f BAMPE --broad -n ${P1} --outdir ${pkfold} -q 0.1 --keep-dup all 2>${pkfold}/${P1}.macs2Peak_summary.txt

else

cd ${bgfold}
P1=$(echo $BG1 | sed 's/.bedgraph/.seacr.peaks/g')

echo "8. Calling SEACR Peaks on ${BG1}"
echo -e "\t bash $seacr ${bgfold}/${BG1} 0.01 non relaxed ${pkfold}/${P1}"

bash $seacr ${bgfold}/${BG1} 0.01 non relaxed ${pkfold}/${P1}

fi

echo -e "Peak calling complete \n"
###

##--------------------------------------------------------------## 9. FASTQC

module unload TrimGalore/0.6.10
module load fastqc
cd ${wdir} 

echo "9. QC on original fastq reads ${R1}"
echo -e "\t fastqc -t ${threads} ${root}${PWD##*/}/${R1} ${root}${PWD##*/}/${R2} -o ${qcfold}"
ls -lh ${root}${PWD##*/}/${R1}
ls -lh ${root}${PWD##*/}/${R2}

fastqc -t ${threads} ${root}${PWD##*/}/${R1} ${root}${PWD##*/}/${R2} -o ${qcfold}
echo -e " ******** All process have been completed ******** \n \n"
###




