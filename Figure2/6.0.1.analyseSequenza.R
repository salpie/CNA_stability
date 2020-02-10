# 6.0.1.analyseSequenza.R
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
library("sequenza")

################# subroutines #################


############### main program ################


#get arguments from script
#arguments <- commandArgs(trailingOnly = TRUE)

#check number of arguments
#if(length(arguments)!=2){
#  stop("\n#### arguments > 6.0.1.analyseSequenza.R <sample list file> <seqFiles> ####\n")
#}

#get sample list infomration
sampleList <- read.csv(file="~/mount3/glandList.exomes.csv", header=TRUE, stringsAsFactors=FALSE)
#sampleList <- read.csv(file=arguments[1], header=TRUE, stringsAsFactors=FALSE)
sampleList <- sampleList[sampleList[["retain"]]!=3 & sampleList[["retain"]]!=0,]
#sampleList <- sampleList[sampleList[[1]]=="C319_R", ]
setNames <- unique(sampleList[[1]])

seqFiles <- "~/projects/glandSeqProject/6.0.0.sequenza/"
#seqFiles <- arguments[2]

#firstRunSwitch <- 1
firstRunSwitch <- 0

#sampleList <- sampleList[sampleList[[1]]=="C267", ]
#sampleList <- sampleList[sampleList[[2]] %in% c("CE2_20", "CE2_18", "CE2_P1_GEC16" , "CE2_P1_GEC19", "CE1_P4_21", "CE1_4_19", "CE1_P3_2", "CE1_P3_1"), ]

for(currSam in 1:nrow(sampleList)){
  
  if(sampleList[currSam, "retain"]==0){
    next()
  }
  
  #current set
  currSet <- sampleList[currSam, "setID"]
  currID <- sampleList[currSam, "sampleID"]
  
  subSample <- sampleList[sampleList[[1]]==currSet, ]
  normalName <- subSample[1, "normalID"]
  
  if(currID == normalName){
    next()
  }
  
  #read data from tumour seqz file name
  seqDirName <- paste(seqFiles, currSet, "/", currID, sep="")
  # 
  # if(firstRunSwitch == 0){
  #   if(currID == normalName | dir.exists(seqDirName)){
  #     next()
  #   }
  #   seqFileName <- paste(seqFiles, currSet, "/", currID, ".seqz.binned.gz", sep="")
  #   if(!file.exists(seqFileName)){
  #     print(paste("#### seqz file missing for ", sampleList[currSam, "sampleID"], " ####",sep=""))
  #     next()
  #   }else{
  #     print(paste("#### performing sequenza analysis for ", sampleList[currSam, "sampleID"], " ####",sep=""))
  #     #stop()
  #   }
  #   newDir <- paste(seqFiles, currSet, "/", currID, "/", sep="")
  #   system(paste("mkdir", newDir))
  #   system(paste("mv", seqFileName, newDir))
  # 
  # }else{
  #   if(currID == normalName){
  #     next()
  #   }
  #   seqFileName <- paste(seqFiles, currSet, "/", currID, "/", currID, ".seqz.binned.gz", sep="")
  #   if(!file.exists(seqFileName)){
  #     print(paste("#### seqz file missing for ", sampleList[currSam, "sampleID"], " ####",sep=""))
  #     next()
  #   }else{
  #     print(paste("#### performing sequenza analysis for ", sampleList[currSam, "sampleID"], " ####",sep=""))
  #     #stop()
  #   }
  # }
  # 
  # now set parameters based on sample (bulk gland)
  seqFileName <- paste(seqFiles, currSet, "/", currID, "/", currID, ".seqz.binned.gz", sep="")
  
  if(sampleList[currSam, "retain"] == "B"){
    gammaParam <- 400
    cellParam <- seq(0.05, 1, 0.01)
  }else if(currSet == "C267" | currSet == "C267_AD" | currSet == "C319" | currSet == "C278"){
    gammaParam <- 80
    cellParam <- seq(0.5, 1, 0.01)
  }else{
    gammaParam <- 35
    cellParam <- seq(0.75, 1, 0.01)
  }

  #analyze data using sequenza
  chromosomes <- paste("chr", c(1:22), sep="")
  seqzExt <- sequenza.extract(file=seqFileName, chromosome.list = chromosomes, gamma = gammaParam, kmin = 25, min.reads.baf = 20)
  
  #infer cellularity and ploidy
  paraSpace <- sequenza.fit(seqzExt, cellularity = cellParam, ploidy = seq(0.8,5,0.1))
  
  #sequenza analysis
  sequenza.results(sequenza.extract = seqzExt, cp.table = paraSpace, sample.id = sampleList[currSam, "sampleID"], out.dir=paste(seqFiles, currSet, "/", currID, sep="") )
  
}




