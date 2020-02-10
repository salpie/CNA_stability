# 2.0.1.vepAnnotatePlatypus.R
# part of the variant calling best practices
#
#  1.0.0.makeMutectScripts.R or 1.0.1.makeMutectScripts.byChr.R (calls and merges all vcfs)
#         |
#         V
#  2.0.0.makePlatypusScripts.R -> 2.0.1.vepAnnotatePlatypus.R
#         |
#         V
#  2.1.0.processVCFplatypus.R (local script: convert merged GATK vcf file to usable .txt table)
#         |
#         V
#  2.2.0.scoreVariantsVEP-VEST.R (local script: annotate variant using annovar and add VEP information)
#         |
#         V
#  2.3.0.filterVariants.R
#         |
#         V
#  2.4.0.getSummaryAndPlot.R
#         |
#         V
#  2.5.0.compileDrivers.R
#         |
#         V
#        ...
#
################### notes ###################
#
#
#
################### libraries ###################

###################   begin   ###################

#arguments <- commandArgs(trailingOnly = TRUE)

#check number of arguments
#if(length(arguments)!=4){
#  stop("\n#### arguments > 2.0.0.makePlatypusScripts.R <sample list file> <platDir> <bcfDir> <scriptsOut> ####\n")
#}

#get sample list infomration
sampleList <- read.csv(file="~/mount3/glandList.genomes.csv", header=TRUE, stringsAsFactors=FALSE)
#sampleList <- read.csv(file="~/projects/gastricProject/sampleList_GC-CK-7901.csv", header=TRUE, stringsAsFactors=FALSE)

#sampleList <- read.csv(file=arguments[1], header=TRUE, stringsAsFactors=FALSE)
#sampleList <- sampleList[sampleList[["retain"]]==2, ]
setNames <- unique(sampleList[["setID"]])

vcfDir <- "/data/BCI-EvoCa-SG/3.0.platypusCalls/genomes/"
#vcfDir <- "/data/BCI-EvoCa2/marnix/data/3.0.platypusCalls/exomes/"
#vcfDir <- arguments[2]

mutectDir <-  "/data/BCI-EvoCa-SG/3.0.mutectCalls/exomes/"
#mutectDir <- "/data/BCI-EvoCa2/marnix/data/3.0.mutectCalls/exomes/"
#mutectDir <- arguments[3]

scriptsOut <- "~/projects/glandSeqProject/A.runScripts/"
#scriptsOut <- "~/projects/gastricProject/A.runScripts/"
#scriptsOut <- arguments[4]

system(command = paste("mkdir ", scriptsOut, "3.0.runPlatypus/", sep=""))

for(currSet in 1:length(setNames)){
  
  subSample <- sampleList[sampleList[["setID"]]==setNames[currSet], ]
  
  system(command = paste0("mkdir ", scriptsOut, "3.0.runPlatypus/", setNames[currSet]))
  
  #output bam file list
  bamFileVect <- paste(subSample[1, "FileHolder2"], "1.3.bamFiles/exomes/", setNames[currSet], "/", subSample[["sampleID"]], "/", subSample[["sampleID"]], ".clipped.bam", sep="")
  bamFileVect <- unique(bamFileVect)
  bamListFileName <- paste(scriptsOut, "3.0.runPlatypus/", setNames[currSet], ".bamList.txt", sep="")
  write.table(matrix(bamFileVect, nrow = length(bamFileVect), ncol=1), file = bamListFileName, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  vcfFiles <- c()
  vcfCounter <- 1
  
  #for each chromosome make a platypus calling script
  for(currChr in c(1:22, "X", "Y")){
    platLogFile <- paste(vcfDir, setNames[currSet], "/", setNames[currSet], "_platypus.chr", currChr,".log", sep="")
    platVCFFile <- paste(vcfDir, setNames[currSet], "/", setNames[currSet], "_platypus.chr", currChr,".vcf", sep="")
    vcfFiles[vcfCounter] <- platVCFFile
    vcfCounter <- vcfCounter + 1
    mutectFile <- paste(mutectDir, setNames[currSet], "/", setNames[currSet], ".mutect.vcf.gz", sep="")
    
    if(subSample[1, "regions"] == "WGS"){
      regionsString <- paste(" \\
--regions=chr", currChr, sep="")
    }else{
      regionsString <- paste(" \\
--regions=", subSample[1, "regions"], ".chr", currChr,".txt", sep="")
    }

    outName <- paste(scriptsOut, "3.0.runPlatypus/", setNames[currSet], "/run_Platypus_", setNames[currSet], "_chr", currChr, ".sh", sep="") 
    totalStrings <- paste("#!/bin/sh
#!/bin/sh
#$ -cwd
#$ -V
#$ -pe smp 1
#$ -l h_rt=100:0:0
#$ -l h_vmem=4G

module load anaconda2
source activate platypus

#--minPosterior=0 \\

mkdir ", vcfDir, setNames[currSet],"/

platypus callVariants \\
--bamFiles=", vcfDir, setNames[currSet],".bamList.txt \\
--output=", platVCFFile," \\
--refFile=", subSample[1, "reference"]," \\
--nCPU=1 \\
--minReads=3 \\
--maxVariants=100 \\
--minMapQual=20 \\
--mergeClusteredVariants=1 \\
--genIndels=1 \\
--genSNPs=1 \\
--filterReadsWithUnmappedMates=1 \\
--trimOverlapping=0 \\
--logFileName=", platLogFile," \\
--source=", mutectFile, regionsString,"

", sep="")
    lapply(totalStrings, write, outName, append=FALSE)
  }
    
    
  #make vcf merging script
  mergedVCF <- paste(vcfDir, setNames[currSet], "/", setNames[currSet], ".merged.vcf", sep="")
  annoVCF <- paste(vcfDir, setNames[currSet], "/", setNames[currSet], ".merged.vepAnno.vcf", sep="")
  
  
  outNameMerge <- paste(scriptsOut, "3.0.runPlatypus/merge_finalVCF_", setNames[currSet], ".sh", sep="") 
  totalStringsMerge <- paste("#!/bin/sh
#$ -cwd
#$ -V
#$ -pe smp 1
#$ -l h_rt=10:0:0      # Request 100 hour runtime
#$ -l h_vmem=2G

module load vcftools

vcf-concat ", paste(vcfFiles, collapse = " "), " > ", mergedVCF,"
  
module load use.dev ensembl-vep/92
  
vep -i ", mergedVCF," \\
-o ", annoVCF," \\
--force_overwrite -format vcf -offline -cache -cache_version 93 dir_cache ~/.vep -assembly GRCh37 -stats_text -stats_text -variant_class -sift b -polyphen b -tab -nearest gene -symbol

", sep="")   
  
  #write seq script file
  lapply(totalStringsMerge, write, outNameMerge, append=FALSE)
    
}

