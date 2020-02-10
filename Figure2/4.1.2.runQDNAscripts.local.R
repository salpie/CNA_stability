# 4.1.2.runQDNAscripts.local.R
# part of the variant calling best practices
# makes QDNA analysis scripts
#
###### notes #######
#
# example: 
# Rscript ~/Code/bestPractices/4.1.2.runQDNAscripts.local.R ~/Projects/glandSeqProject/sampleList.csv /data/BCI-EvoCa2/wchc/glandSeqProject/2.processedBams/ ~/Projects/glandSeqProject/runScripts/
#
# note requires bam files to be organised into a cancer/sample/file dir structre
#
#
#
###### begin ########

library(QDNAseq)
#library(robustbase)


#get arguments from script
#arguments <- commandArgs(trailingOnly = TRUE)

#check number of arguments
#if(length(arguments)!=3){
#  stop("\n#### arguments > 4.1.2.runQDNAscripts.local.R <sample list file> <bamDir> <outDir> ####\n")
#}

#get sample list infomration
#sampleList <- read.csv(file="~/projects/glandSeqProject/glandList.LPWGS.csv", header=TRUE, stringsAsFactors=FALSE)
sampleList <- read.csv(file="~/mount3/glandList.LPWGS.csv", header=TRUE, stringsAsFactors=FALSE)
sampleList <- sampleList[sampleList[[1]]=="C278", ]

#sampleList <- read.csv(file=arguments[1], header=TRUE, stringsAsFactors=FALSE)
#sampleList <- sampleList[sampleList[["retain"]]==1, ]
#sampleList <- sampleList[sampleList[["setID"]]=="C277", ]

runIDs <- unique(sampleList[["sampleInfo"]])
#runIDs <- "run11-18"

bamDir <- "~/mount3/1.3.bamFiles/LPWGS/"
#bamDir <- "/Volumes/EvoCaBackup/dataStore/LPWGS/"
#bamDir <- "projects/singleCellLPWGS/dataDump/"
#bamDir <- arguments[2]

outDir <- "~/projects/glandSeqProject/4.1.QDNAseq/LPWGS/"
#outDir <- "projects/singleCellLPWGS/4.1.QDNAseq/LPWGS/"
#outDir <- arguments[3]

#set binsize and download
#binSize <- 5000
#binSize <- 50
#downloadedBins <- getBinAnnotations(binSize = 50, genome = "hg19")
downloadedBins <- readRDS(file="~/projects/referenceFiles/refQDNAseq")

for(currRun in 1:length(runIDs)){
  
  runID <- runIDs[currRun]
  sampleListSub <- sampleList[sampleList[["sampleInfo"]]==runID, ]
  
  setNames <- unique(sampleListSub[["setID"]])
  
  for(currSet in 1:length(setNames)){
    print(paste("####### performing QDNA analysis for sample", setNames[currSet], " batch ", runID, " #########"))
    
    #begin processing
    system(paste("mkdir ", outDir, setNames[currSet] , sep=""))
    
    subSample <- sampleListSub[sampleListSub[["setID"]]==setNames[currSet], ]
    
    progVect <- seq(1, nrow(subSample), length.out = 9)
    progVect <- round(progVect, digits = 0)
    
    print(paste("####### getting bams #########"))
    bamDataFile <- paste(outDir, setNames[currSet], "/", setNames[currSet], ".", runID, ".bams.Rdata", sep="")
    if(!file.exists(bamDataFile)){
      binCountsList <- as.list(NA)
      for(currSam in 1:nrow(subSample)){
        if(currSam %in% progVect){
          print(paste("#### ", which(currSam == progVect), "0 % complete ####", sep=""))
        }
        
        #save read counts to table
        inputBam <- paste(bamDir, setNames[currSet], "/", subSample[currSam, "sampleID"], "/", subSample[currSam, "sampleID"], ".bam", sep="")
        #sampleList[sampleList[["sampleID"]]==subSample[currSam, "sampleID"], "reads"] <- system(command=paste("~/Software/samtools-0.1.19/samtools view -c -q 20 ", outputBam, sep=""))
        #inputBam <- paste(bamDir, subSample[currSam, "sampleID"], "_dedup.bam", sep="")
        
        binCountsList[[currSam]] <- binReadCounts(downloadedBins, bamfiles=inputBam, isPaired = TRUE)
      }
      
      #save bin counts in Robject
      save(binCountsList, file = paste(outDir, setNames[currSet], "/", setNames[currSet], ".", runID, ".bams.Rdata", sep=""))

    }else{
      load(file = bamDataFile)
    }
    
    
    #print filtered reads as isobar
    readCountsFiltered <- as.list(NA)
    graph2 <- paste(outDir, setNames[currSet], "/", setNames[currSet], ".", runID, ".isobar.pdf", sep="")
    pdf(file=graph2, width=5, height=5)
    for(currSam in 1:length(binCountsList)){
      readCountsFiltered[[currSam]] <- applyFilters(binCountsList[[currSam]], residual=TRUE, blacklist=TRUE, mappability = 60)
      isobarPlot(readCountsFiltered[[currSam]])
    }
    dev.off()
    
    
    #plot raw filtered read counts
    graph3 <- paste(outDir, setNames[currSet], "/", setNames[currSet], ".", runID, ".raw.bins.pdf", sep="")
    pdf(file=graph3, width=15, height=5)
    par(ps=5, cex=2, mar=c(3,3,3,3))
    for(currSam in 1:length(binCountsList)){
      plot(readCountsFiltered[[currSam]], log=FALSE)
    }
    dev.off()
    
    
    #apply the correction for GC content and mappability, then normalize, smooth outliers, calculate segmentation 
    #and plot the copy number profile
    graph4 <- paste(outDir, setNames[currSet], "/", setNames[currSet], ".", runID, ".normalised.bins.pdf", sep="")
    copyNumbers <- as.list(NA)
    copyNumbersNormalized <- as.list(NA)
    copyNumbersSmooth <- as.list(NA)
    pdf(file=graph4, width=15, height=5)
    par(ps=5, cex=2, mar=c(3,3,3,3))
    for(currSam in 1:length(binCountsList)){
      tempData <- estimateCorrection(readCountsFiltered[[currSam]], span=0.65,  family="symmetric", adjustIncompletes=TRUE)
      copyNumbers[[currSam]] <- correctBins(tempData)
      copyNumbersNormalized[[currSam]] <- normalizeBins(copyNumbers[[currSam]])
      copyNumbersSmooth[[currSam]] <- smoothOutlierBins(copyNumbersNormalized[[currSam]], logTransform = FALSE)
      plot(copyNumbersSmooth[[currSam]], log=TRUE)
    }
    dev.off()
    
    
    
    #plot final segmented graph
    graph5 <- paste(outDir, setNames[currSet], "/", setNames[currSet], ".", runID, ".segments.pdf", sep="")
    copyNumbersSegmented <- as.list(NA)
    copyNumbersSeg <- as.list(NA)
    pdf(file=graph5, width=15, height=6)
    for(currSam in 1:length(binCountsList)){
      par(ps=4, cex=2.2, mar=c(3,3,3,3))
      copyNumbersSegmented[[currSam]] <- segmentBins(copyNumbersSmooth[[currSam]], alpha=0.01, smoothBy=FALSE, transformFun="log2")
      #copyNumbersSeg[[currSam]] <- callBins(copyNumbersSegmented[[currSam]], organism="human", method = "cutoff")
      
      #copyNumbersSeg[[currSam]] <- callBins(copyNumbersSegmented[[currSam]], organism="human")
      plot(copyNumbersSegmented[[currSam]], main=subSample[currSam, "sampleID"], ylim=c(0,4), logTransform = FALSE)
    }
    dev.off()
    
    
    #export segmentation table
    system(paste("mkdir ", outDir, setNames[currSet], "/segTables/", sep=""))
    for(currSam in 1:length(copyNumbersSegmented)){
      fileName <- paste(outDir, setNames[currSet], "/segTables/", subSample[currSam, "sampleID"], ".segmentations.txt", sep="")
      exportBins(copyNumbersSegmented[[currSam]], file=fileName, format = "igv", type="segments", filter=TRUE, logTransform = FALSE)
    }
    
    #merge rows in segmentation tables to get true segmentation
    for(currSam in 1:length(copyNumbersSegmented)){
      fileName <- paste(outDir, setNames[currSet], "/segTables/", subSample[currSam, "sampleID"], ".segmentations.txt", sep="")
      tempTab <- read.table(file = fileName, header = TRUE, stringsAsFactors = FALSE, sep = "\t")
      names(tempTab)[5] <- "ratio"
      
      rowCounter <- 1
      remRow <- c()
      for(i in 1:nrow(tempTab)){
        if(i == nrow(tempTab)){
          break()
        }else{
          if(tempTab[i, "ratio"] == tempTab[i+1, "ratio"] & tempTab[i, "chromosome"] == tempTab[i+1, "chromosome"]){
            remRow[rowCounter] <- i+1
            rowCounter <- rowCounter + 1
          }
        }
      }
      #get annotation rows and populate new table
      newTab <- data.frame(matrix(NA, nrow = 0, ncol = 5))
      names(newTab) <- names(tempTab)
      
      annoRows <- c(1:nrow(tempTab))
      annoRows <- annoRows[!(annoRows %in% remRow)]
      for(add in 1:length(annoRows)){
        newTab[add, "chromosome"] <- tempTab[annoRows[add], "chromosome"]
        newTab[add, "start"] <- tempTab[annoRows[add], "start"]
        newTab[add, "ratio"] <- tempTab[annoRows[add], "ratio"]
        if(add == length(annoRows)){
          newTab[add, "end"] <- tempTab[nrow(tempTab), "end"]
        }else{
          newTab[add, "end"] <- tempTab[(annoRows[add+1]-1), "end"]
        }
        newTab[add, "feature"] <- paste(newTab[add, "chromosome"], ":", newTab[add, "start"], "-", newTab[add, "end"], sep="")
      }
      fileName <- paste(outDir, setNames[currSet], "/segTables/", subSample[currSam, "sampleID"], ".merged.segs.txt", sep="")
      write.table(newTab, file = fileName, sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
      
    }
    
  }
}



#write.table(sampleList, file="~/Projects/glandSeqProject/LPWGS.summary.txt", sep="\t", header=TRUE, quote = FALSE)


