# 4.1.0.makeQDNAscripts.R
# part of the variant calling best practices
# 
# plot karyotypes from QDNA results for ONE SET
#
###### notes #######
#
# 4.1.0.makeQDNAscripts.R
#         |
#         V
# 4.1.1.plotQDNAresults.R (using segmentations.txt)
#
#
###### begin ########

#get arguments from script
arguments <- commandArgs(trailingOnly = TRUE)

#check number of arguments
if(length(arguments)!=4){
  stop("\n#### arguments > 4.1.0.makeQDNAscripts.R <sample list file> <bamDir> <outDir> <scriptsOut> ####\n")
}

#get sample list infomration
#sampleList <- read.csv(file="~/Projects/glandSeqProject/glandList.LPWGS.csv", header=TRUE, stringsAsFactors=FALSE)
sampleList <- read.csv(file=arguments[1], header=TRUE, stringsAsFactors=FALSE)
sampleList <- sampleList[sampleList[["retain"]]==1, ]

sampleList <- sampleList[sampleList[["sampleInfo"]]=="run11-18", ]


#bamDir <- "/data/BCI-EvoCa-SG/1.3.bamFiles/LPWGS/"
bamDir <- arguments[2]

#outDir <- "/data/BCI-EvoCa-SG/4.1.QDNAseq/LPWGS/"
outDir <- arguments[3]

#scriptsOut <- "~/Projects/glandSeqProject/A.runScripts/"
scriptsOut <- arguments[4]

system(command = paste("mkdir ", scriptsOut, "4.1.QDNAseq/", sep=""))

for(currSet in 1:nrow(sampleList)){
  
  print(paste("####### making QDNA script for sample", sampleList[currSet, "sampleID"], "#########"))
  
  newDir <- paste(outDir, sampleList[currSet, "setID"], "/", sampleList[currSet, "sampleID"], "/", sep="")
  
  RscriptStore <- paste(scriptsOut, "4.1.QDNAseq/runQDNA_", sampleList[currSet, "sampleID"], ".R", sep="")
  outRscript <- paste(outDir, sampleList[currSet, "setID"], "/runQDNA_", sampleList[currSet, "sampleID"], ".R", sep="") 
  
  
  
  ############ make .sh scripts to run Rscripts on the cluster  ############  
  
  outName <- paste(scriptsOut, "4.1.QDNAseq/runQDNA_", sampleList[currSet, "setID"], "_", sampleList[currSet, "sampleID"], ".sh", sep="") 
  totalStrings <- paste("#!/bin/sh
#$ -cwd
#$ -V
#$ -pe smp 1
#$ -l h_rt=100:0:0
#$ -l h_vmem=5G

mkdir ", newDir,"

Rscript ", outRscript,"
", sep="")
      
  #write seq script file
  lapply(totalStrings, write, outName, append=FALSE)

  
  
  ############ make R scripts to run QDNA ############  
  
  outputBam <- paste(bamDir, sampleList[currSet, "setID"], "/", sampleList[currSet, "sampleID"], "/", sampleList[currSet, "sampleID"], ".sorted.bam", sep="")
  graph1 <- paste(newDir, sampleList[currSet, "sampleID"], ".raw_profile.pdf", sep="")
  graph2 <- paste(newDir, sampleList[currSet, "sampleID"], ".isobar.pdf", sep="")
  graph3 <- paste(newDir, sampleList[currSet, "sampleID"], ".noise.pdf", sep="")
  graph4 <- paste(newDir, sampleList[currSet, "sampleID"], ".copy_number_profile.pdf", sep="")
  graph5 <- paste(newDir, sampleList[currSet, "sampleID"], ".segments.pdf", sep="")
  graph6 <- paste(newDir, sampleList[currSet, "sampleID"], ".copy_numbercalls.pdf", sep="")
  
  
  totalStrings <- paste("
# Rscript runs specified bam file through the QDNA pipeline
# QDNA calls CNAs, from low pass whole genome sequencing data

#libraries
library(QDNAseq)

downloadedBins <- getBinAnnotations(binSize=500)

readCounts <- binReadCounts(downloadedBins, bamfiles = \"",outputBam,"\")

#print raw copy profile
pdf(file=\"", graph1,"\" , 5, 5)
  plot(readCounts, logTransform=FALSE, ylim=c(-10, 50), main=\"Raw Profile\")
  highlightFilters(readCounts, logTransform=FALSE, residual=TRUE, blacklist=TRUE)
dev.off()

#print filtering 
readCountsFiltered <- applyFilters(readCounts, residual=TRUE, blacklist=TRUE)
pdf(file=\"", graph2,"\" , 5, 5)
  isobarPlot(readCountsFiltered)
dev.off()
  
#print correction factor
readCountsFiltered <- estimateCorrection(readCountsFiltered)
pdf(file=\"", graph3,"\" , 5, 5)
  noisePlot(readCountsFiltered)
dev.off()
 
#apply the correction for GC content and mappability which we then normalize, smooth outliers, calculate segmentation 
#and plot the copy number profile
copyNumbers <- correctBins(readCountsFiltered)
copyNumbersNormalized <- normalizeBins(copyNumbers)
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized)

pdf(file=\"", graph4,"\" , 5, 5)
  plot(copyNumbersSmooth)
dev.off()

copyNumbersSegmented <- segmentBins(copyNumbersSmooth)
copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)

pdf(file=\"", graph5,"\" , 5, 5)
  plot(copyNumbersSegmented)
dev.off()

copyNumbersCalled <- callBins(copyNumbersSegmented)
pdf(file=\"", graph6,"\" , 5, 5)
  plot(copyNumbersCalled)
dev.off()

", sep="")
  
  #write seq script file
  lapply(totalStrings, write, RscriptStore, append=FALSE)

  
}






