# 6.0.0.makeSequenzaScripts.R
# part of the variant calling best practices
#
#
#
################### notes ###################
#
# run sequenza locally
#
# 6.0.0.makeSequenzaScripts.R
#           |
#           V
# 6.0.1.analyseSequenza.R
#           |
#           V
# 6.0.2.plotKaryotypes.R
#           |
#           V
# 6.1.0.runSequenzaJointly.R
#
#
################# libraries #################


################# subroutines #################


############### main program ################


#get arguments from script
arguments <- commandArgs(trailingOnly = TRUE)

#check number of arguments
if(length(arguments)!=4){
  stop("\n#### arguments > 6.0.0.makeSequenzaScripts.R <sample list file> <bamFiles> <seqFiles> <scriptOut> ####\n")
}

#get sample list infomration
sampleList <- read.csv(file="~/mount3/glandList.genomes.csv", header=TRUE, stringsAsFactors=FALSE)
sampleList <- read.csv(file=arguments[1], header=TRUE, stringsAsFactors=FALSE)
sampleList <- sampleList[sampleList[["retain"]]!=3 & sampleList[["retain"]]!=0,]
setNames <- unique(sampleList[[1]])

bamFiles <- "/data/BCI-EvoCa-SG/1.3.bamFiles/genomes/"
#bamFiles <- "/data/BCI-EvoCa2/marnix/data/1.3.bamFiles/"
bamFiles <- arguments[2]

seqLocal <- "~/projects/glandSeqProject/6.0.0.sequenza/"

seqFiles <- "/data/BCI-EvoCa-SG/6.0.0.sequenza/"
seqFiles <- arguments[3]

scriptDir <- "~/projects/glandSeqProject/A.runScripts/6.0.0.sequenza/"
scriptDir <- arguments[4]

GCcontentFile <- "/data/BCI-EvoCa2/wchc/referenceFiles/hg19.gc50Base.txt.gz"
refFile <- "/data/BCI-EvoCa2/wchc/referenceFiles/ucsc.hg19.fasta"

system(paste("mkdir ", scriptDir, sep=""))

#make sequenza prep scripts (for each sample)
for(currSam in 1:nrow(sampleList)){
  print(paste("#### making sequenza script for sample ", sampleList[currSam, "sampleID"], " ####",sep=""))
  
  sampleID <- sampleList[currSam, "sampleID"]
  
  #current sample
  setID <- sampleList[currSam, "setID"]
  
  #current normal
  subSample <- sampleList[sampleList[["setID"]]==setID, ]
  
  #normal name
  normalName <- subSample[1,"normalID"]
  
  #bam file names
  tumourBam <- paste(bamFiles, setID, "/", sampleID, "/", sampleID, ".mkdub.bam", sep="")
  normBam <- paste(bamFiles, setID, "/", normalName, "/", normalName, ".mkdub.bam", sep="")
  
  #sequenza files
  seqFile <- paste(seqFiles, setID, "/", sampleID, ".seqz.gz", sep="")
  binnedSeqFile <- paste(seqFiles, setID, "/", sampleID, ".seqz.binned.gz", sep="")
  
  localSeqFile <- paste(seqLocal, setID, "/", sampleID, sep="")
  
  if(sampleID == normalName){
    next()
  }else if(dir.exists(localSeqFile)){
    print(paste("#### seqs file already exists for sample ", sampleList[currSam, "sampleID"], " ####",sep=""))
    #next()
  }
  
  #prep .sh file name
  SHfile <- paste(scriptDir, "sequenzaPrep.", setID, ".", sampleID, ".sh", sep="")
  
  prepSHstrings <- paste("#!/bin/sh
#$ -cwd
#$ -V
#$ -pe smp 1        # Request 1 CPU core1
#$ -l h_vmem=6G    # Request 20GB RAM / core, i.e. 24GB total
#$ -l h_rt=120:0:0   # Request 48 hour runtime

#. /data/home/mpx155/sequenza/bin/activate

#make analysis dir
mkdir ", seqFiles, setID,"

#chmod 777 ", seqFiles, setID,"

#get .seqz file
./bin/sequenza-utils.py bam2seqz -gc ", GCcontentFile," --fasta ", refFile," -n ", normBam," -t ", tumourBam," | gzip > ", seqFile,"

#bin .seqz file to shorten analysis time
./bin/sequenza-utils.py seqz-binning -w 250 -s ", seqFile," | gzip > ", binnedSeqFile,"

#chmod 777 ", seqFile,"
#chmod 777 ", binnedSeqFile,"


#run pre-processed files using sequenza in R locally
", sep="")

  #write prep .sh file
  lapply(prepSHstrings, write, SHfile, append=FALSE)
     
}

