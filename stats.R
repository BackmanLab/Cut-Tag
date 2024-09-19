#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#---------------------------------------------------------------#
# Author: Lucas Carter                                                   
# Email: lucascarter2025@u.northwestern.edu                                                          
# PI: Vadim Backman                                                   
# Description: 
# This script generates statistics on alignment, duplication,
# and other facets of Cut&Tag data, post processing. A more 
# thorough description can be found at README.md
#
# Usage: module load R/4.3.0 
# Rscript --vanilla --verbose stats.R <root/file/path> <resultname.txt (optional)>
# Call this script at end of Cut&Tag data processing once all FASTQs are preprocessed
#---------------------------------------------------------------#

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("Path to root analysis directory must be supplied ('path/to/dir').n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "statistics.txt"
}

output.path = file.path(args[1],'results')
cat(paste0("output path: ",output.path, "\n","\n"))

##-------------------------------------------------------## load packages

require(dplyr)
require(tidyr)
require(ggplot2)
require(viridis)
require(ggpubr)
require(corrplot)

##-------------------------------------------------------## Get alignment summary

sam.path <- file.path(args[1],"SAM")
hists <- list.files(sam.path, pattern = ".bowtie2.txt")

alignResult = c()
for(hist in hists){
  alignRes = read.table(file.path(sam.path, hist), header = FALSE, fill = TRUE)
  alignRate = substr(alignRes$V1[6], 1, nchar(as.character(alignRes$V1[6]))-1)
  histInfo = strsplit(hist, "_")[[1]]
  alignResult = data.frame(Target = histInfo[1], Group = histInfo[2], 
                           SequencingDepth = alignRes$V1[1] %>% as.character %>% as.numeric, 
                           MappedFragNum_hg38 = alignRes$V1[4] %>% as.character %>% as.numeric + alignRes$V1[5] %>% as.character %>% as.numeric, 
                           AlignmentRate_hg38 = alignRate %>% as.numeric)  %>% rbind(alignResult, .)
}

##-------------------------------------------------------## Get duplication summary
## Summarize the duplication information from the picard summary outputs.

hists <- list.files(sam.path, pattern = ".dupMark.txt")

dupResult = c()
for(hist in hists){
  dupRes = read.table(file.path(sam.path,hist), header = TRUE, fill = TRUE)
  
  histInfo = strsplit(hist, "_")[[1]]
  dupResult = data.frame(Target = histInfo[1], Group = histInfo[2], MappedFragNum_hg38 = dupRes$READ_PAIRS_EXAMINED[1] %>% as.character %>% as.numeric, DuplicationRate = dupRes$PERCENT_DUPLICATION[1] %>% as.character %>% as.numeric * 100, EstimatedLibrarySize = dupRes$ESTIMATED_LIBRARY_SIZE[1] %>% as.character %>% as.numeric) %>% mutate(UniqueFragNum = MappedFragNum_hg38 * (1-DuplicationRate/100))  %>% rbind(dupResult, .)
}

alignDupSummary = left_join(alignResult, dupResult, by = c("Target", "Group", "MappedFragNum_hg38")) %>% mutate(DuplicationRate = paste0(DuplicationRate, "%"))

write.table(data.frame(alignDupSummary), file=file.path(output.path,args[2]), sep= "\t", row.names=FALSE)

##-------------------------------------------------------## Get fragment size

hists <- list.files(sam.path, pattern = "fragmentLen.txt")

fragLen = c()
for(hist in hists){
  
  histInfo = strsplit(hist, "_")[[1]]
  fragLen = read.table(file.path(sam.path,hist), header = FALSE) %>% mutate(fragLen = V1 %>% as.numeric, fragCount = V2 %>% as.numeric, Weight = as.numeric(V2)/sum(as.numeric(V2)), Target = histInfo[1], Group = histInfo[2], sampleInfo = hist) %>% rbind(fragLen, .) 
}

fragLen$sampleInfo = factor(fragLen$sampleInfo, levels = hists)
fragLen$Target = factor(fragLen$Target, levels = unique(fragLen$Target))

## Generate the fragment size density plot (violin plot)
fig1 = fragLen %>% ggplot(aes(x = sampleInfo, y = fragLen, weight = Weight, fill = Target)) +
  geom_violin(bw = 5) +
  scale_y_continuous(breaks = seq(0, 800, 50)) +
  scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
  theme_bw(base_size = 20) +
  ggpubr::rotate_x_text(angle = 20) +
  ylab("Fragment Length") +
  xlab("")+ theme(legend.position="bottom", plot.title = element_text(size=10), text = element_text(size=10, family="Arial"))

ggsave(file.path(output.path,"fraglen.viol.png"), plot = fig1, width = 20, height = 20, units = 'cm')

fig2 = fragLen %>% ggplot(aes(x = fragLen, y = fragCount, color = Target, group = sampleInfo, linetype = Group)) +
  geom_line(size = 1) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
  theme_bw(base_size = 20) +
  xlab("Fragment Length") +
  ylab("Count") +
  coord_cartesian(xlim = c(0, 500))+  theme(legend.position="bottom", plot.title = element_text(size=10), text = element_text(size=10, family="Arial"))

ggsave(file.path(output.path,"fraglen.lin.png"), plot = fig2, width = 20, height = 20, units = 'cm')

##-------------------------------------------------------## Check correlation between samples

bed.path <- file.path(args[1],"BED")
hists <- list.files(bed.path, pattern = ".bin.bed")

reprod = c()
fragCount = NULL
for(hist in hists){
  
  if(is.null(fragCount)){
    
    fragCount = read.table(file.path(bed.path, hist), header = FALSE) 
    colnames(fragCount) = c("chrom", "bin", hist)
    
  }else{
    
    fragCountTmp = read.table(file.path(bed.path, hist), header = FALSE)
    colnames(fragCountTmp) = c("chrom", "bin", hist)
    fragCount = full_join(fragCount, fragCountTmp, by = c("chrom", "bin"))
    
  }
}

M = cor(fragCount %>% select(-c("chrom", "bin")) %>% log2(), use = "complete.obs") 
write.table(data.frame(M), file=file.path(output.path,"corrmatrix.txt"), sep= "\t", row.names=FALSE)

# Initialize file path
file_path=file.path(output.path,"corrmat.png")
png(height=2000, width=2000, file=file_path, type = "cairo")

corrplot(M, method = "color", outline = T, addgrid.col = "darkgray", order="hclust", addrect = 3, rect.col = "black", rect.lwd = 3,cl.pos = "b", tl.col = "indianred4", tl.cex = 0.5, cl.cex = 0.5, addCoef.col = "black", number.digits = 2, number.cex = 0.5, col = colorRampPalette(c("midnightblue","white","darkred"))(100))

dev.off()

##-------------------------------------------------------## Check peak number

peak.path <- file.path(args[1],"peaks")
hists <- list.files(peak.path, pattern = ".narrowPeak")

peakN = c()
peakWidth = c()
peakType = c("control", "top0.01")
for(hist in hists){
  histInfo = strsplit(hist, "_")[[1]]
  if(histInfo[1] != "IgG"){
    for(type in peakType){
      peakInfo = read.table(file.path(peak.path,hist), header = FALSE, fill = TRUE)  %>% mutate(width = abs(V3-V2))
      peakN = data.frame(peakN = nrow(peakInfo), peakType = type, Histone = histInfo[1], Replicate = histInfo[2]) %>% rbind(peakN, .)
      peakWidth = data.frame(width = peakInfo$width, peakType = type, Histone = histInfo[1], Replicate = histInfo[2])  %>% rbind(peakWidth, .)
    }
  }
}

write.table(data.frame(peakN), file=file.path(output.path,"peakstats.txt"), sep= "\t", row.names=FALSE)
