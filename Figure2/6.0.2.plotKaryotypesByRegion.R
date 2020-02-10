# 6.0.2.plotKaryotypes.R
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
#library("sequenza")

################# subroutines #################

getCol <- function(subVect){
  if(subVect[1, "A"]==1 & subVect[1, "B"]==1){
    colSub <- "grey95"
  }else if(subVect[1, "A"]==2 & subVect[1, "B"]==0){
    colSub <- "steelblue"
  }else if(subVect[1, "A"]==1 & subVect[1, "B"]==0){
    colSub <- "lightblue"
  }else if(subVect[1, "A"]==2 & subVect[1, "B"]==1){
    colSub <- "red"
  }else if(subVect[1, "A"]==2 & subVect[1, "B"]==2){
    colSub <- "darkred"
  }else{
    colSub <- "brown"
  }
  return(colSub)
}

#tempTabSub <- tempTab; chromPosSub <- chromPos
getGenomeCoord <- function(tempTabSub, chromPosSub){
  chrVect <- tempTabSub[["chr"]]
  chrRows <- chromPosSub[match(chrVect, chromPosSub[["chromosome"]]), "start"]
  tempTabSub["start"] <- tempTabSub[["start"]] + chrRows
  tempTabSub["end"] <- tempTabSub[["end"]] + chrRows
  return(tempTabSub)
}



############### main program ################


#get sample list infomration
exomeList <- read.csv(file="~/mount3/7.0.finalFigures/tables/glandList.exomes.csv", header=TRUE, stringsAsFactors=FALSE)
exomeList <- exomeList[exomeList[["retain"]]!=3 & exomeList[["retain"]]!=0, ]
setNames <- unique(exomeList[["setID"]])

lpwgsList <- read.csv(file="~/mount3/glandList.LPWGS.csv", header=TRUE, stringsAsFactors=FALSE)
lpwgsList <- lpwgsList[lpwgsList[["retain"]]==1, ]

vcfDir <- "~/projects/glandSeqProject/3.0.platypusCalls/exomes/"

#cnaExomesDir <- "~/projects/glandSeqProject/6.0.0.sequenza/"
cnaExomesDir <- "~/projects/glandSeqProject/6.1.0.sequenzaJointCalls/"

rerunFlag <- ""
#rerunFlag <- ".rerun"

lpwgsDir <- "~/projects/glandSeqProject/4.1.QDNAseq/LPWGS/"

#get driver information
driverInfo <- read.csv(file="~/projects/glandSeqProject/COSMIC_Census_04_2018_allEpithelial_Tier1.csv", stringsAsFactors = FALSE, header = TRUE)

workingDir <- "~/projects/glandSeqProject/6.0.2.plotKaryotypes/"

#chromosome positions
chromPos <- read.table(file="~/projects/referenceFiles/chromosomeStartEnd.csv", stringsAsFactors = FALSE, header=TRUE, sep=",")
chromPos["chromosome"] <- paste("chr", chromPos[["chromosome"]], sep="")
chromPos[23, "chromosome"] <- "chrX"
chromPos[24, "chromosome"] <- "chrY"

centroPos <- read.table(file="~/projects/referenceFiles/centromerePositions.csv", stringsAsFactors = FALSE, header=TRUE, sep=",")
chromPos <- cbind(chromPos, centroPos)
chromPos["CENTROSTART"] <- chromPos[["CENTROSTART"]] + chromPos[["start"]]
chromPos["CENTROEND"] <- chromPos[["CENTROEND"]] + chromPos[["start"]]

for(currSet in 1:length(setNames)){
  
  subLpwgs <- lpwgsList[lpwgsList[["setID"]]==setNames[currSet],]
  
  if(nrow(subLpwgs)==0 | currSet == 5){
    next()
  }
  
  #subsample list
  subExomes <- exomeList[exomeList[["setID"]]==setNames[currSet], ] 
  subExomes <- subExomes[subExomes[["retain"]]!="B", ]
  bioListEx <- subExomes[["sampleID"]]
  bioListEx <- bioListEx[bioListEx!=subExomes[1, "normalID"]]
  
  bioListLp <- subLpwgs[["sampleID"]]
  bioListLp <- bioListLp[bioListLp!=subLpwgs[1, "normalID"]]
  
  bioListTotal <- c(bioListEx, bioListLp)
  
  #read in plot table
  plotList <- read.csv(file=paste(workingDir, setNames[currSet], "/", setNames[currSet], "_glandlist.csv", sep=""), stringsAsFactors = FALSE)
  plotList <- plotList[plotList[["ID"]] %in% bioListTotal, ]
  
  
  stateList <- read.csv(file=paste(workingDir, setNames[currSet], "/", setNames[currSet], "_copyStateRef.csv", sep=""), stringsAsFactors = FALSE, header = TRUE)
  #stateList <- stateList[stateList[["type"]]=="monosomy" | stateList[["type"]]=="trisomy",]
  
  #popluate segmentation list
  segmentationList <- as.list(NA)
  segCounter <- 1
  
  gradFit <- as.list(NA)
  dataStore <- as.list(NA)
  gradCounter <- 1
  
  for(currBio in 1:nrow(plotList)){
    currID <- plotList[currBio, "ID"]
    
    if(plotList[currBio, "type"] == "bulk"){
      #read in sequenza table
      cnaIn <- read.table(file=paste(cnaExomesDir, setNames[currSet], "/", currID, "/", currID, "_segments.txt", sep=""), header = TRUE, sep = "\t") 
      tempTab <- cnaIn[c("chromosome", "start.pos", "end.pos", "A", "B")]
      names(tempTab) <- c("chr", "start", "end", "A", "B")
      
      #convert coordinates to genome
      tempTabGen <- getGenomeCoord(tempTab, chromPos)
      
      segmentationList[[segCounter]] <- tempTabGen
      names(segmentationList)[[segCounter]] <- currID
      segCounter <- segCounter + 1
      
    }else if(plotList[currBio, "type"] == "exome"){
      #read in sequenza table
      cnaIn <- read.table(file=paste(cnaExomesDir, setNames[currSet], "/", currID, "/", currID, "_segments.txt", sep=""), header = TRUE, sep = "\t") 
      tempTab <- cnaIn[c("chromosome", "start.pos", "end.pos", "A", "B")]
      names(tempTab) <- c("chr", "start", "end", "A", "B")
      
      #convert coordinates to genome
      tempTabGen <- getGenomeCoord(tempTab, chromPos)
      
      segmentationList[[segCounter]] <- tempTabGen
      names(segmentationList)[[segCounter]] <- currID
      segCounter <- segCounter + 1
      
    }else if(plotList[currBio, "type"] == "lpwgs"){
      
      cnaInFile <- paste(lpwgsDir, setNames[currSet], "/segTables/", currID, ".merged.segs.txt", sep="")
      if(!file.exists(cnaInFile)){
        next()
      }
      
      #read in QDNAseq table
      cnaIn <- read.table(file=cnaInFile, header = TRUE, sep = "\t") 
      cnaIn["chromosome"] <- paste("chr", cnaIn[["chromosome"]], sep="")
      cnaIn["cn"] <- NA
      
      #get gain ratio
      cnaRatio <- data.frame(matrix(NA, nrow = nrow(stateList), ncol = 2))
      names(cnaRatio) <- c("ratio", "state")
      for(i in 1:nrow(cnaRatio)){
        #if several segmentations are available remove those that are < 1Mb and calculate mean
        tempValues <- cnaIn[cnaIn[["chromosome"]] == stateList[i, "chromosome"], ]
        cnaLength <- tempValues[["end"]] - tempValues[["start"]]
        cnaRatio[i, "ratio"] <- mean(tempValues[cnaLength>5000000, "ratio"])
        cnaRatio[i, "state"] <- stateList[i, "state"]
      }
      
      #get gradient of fit
      gradFit[[gradCounter]] <- lm(data = cnaRatio, formula = state ~ ratio)
      
      for(curr in 1:nrow(cnaIn)){
        cnaIn[curr, "cn"] <- round(cnaIn[curr, "ratio"] * gradFit[[gradCounter]]$coefficients[2] + gradFit[[gradCounter]]$coefficients[1], digits = 0)
      }
      if(sum(cnaIn[["ratio"]]) == 0){
        next()
      }
      cnaIn <- cnaIn[cnaIn[["cn"]] > 0 & cnaIn[["cn"]] < 8, ]
      
      #exclude samples where the ratio range is too low (i.e cellularity is low)
      #if(mean(gainRatio) < 0.05){
      #  next()
      #}
      
      #save table
      newOut <- paste(lpwgsDir, setNames[currSet], "/segTables/", currID, ".merged.segs.cn.txt", sep="")
      write.table(x = cnaIn, file = newOut, sep="\t", quote = FALSE, row.names = FALSE)
      
      dataStore[[gradCounter]] <- cnaIn
      names(dataStore)[[gradCounter]] <- currID
      gradCounter <- gradCounter + 1
      
      cnaIn["A"] <- NA
      cnaIn["B"] <- NA
      for(i in 1:nrow(cnaIn)){
        if(cnaIn[i, "cn"] == 2){
          cnaIn[i, "A"] <- 1
          cnaIn[i, "B"] <- 1 
        }else if(cnaIn[i, "cn"] == 1){
          cnaIn[i, "A"] <- 1
          cnaIn[i, "B"] <- 0
        }else if(cnaIn[i, "cn"] == 3){
          cnaIn[i, "A"] <- 2
          cnaIn[i, "B"] <- 1
        }else if(cnaIn[i, "cn"] == 4){
          cnaIn[i, "A"] <- 2
          cnaIn[i, "B"] <- 2
        }else{
          cnaIn[i, "A"] <- 2
          cnaIn[i, "B"] <- 3
        }
      }
      tempTab <- cnaIn[c("chromosome", "start", "end", "A", "B")]
      names(tempTab) <- c("chr", "start", "end", "A", "B")
      tempTabGen <- getGenomeCoord(tempTab, chromPos)
      
      segmentationList[[segCounter]] <- tempTabGen
      names(segmentationList)[[segCounter]] <- currID
      segCounter <- segCounter + 1
    }
  }  
  
  ########## plot linear model fit and store gradients ######
  resTable <- data.frame(matrix(NA, ncol = 3, nrow = length(gradFit)))
  names(resTable) <- c("sampleID", "p.value", "gradient")
  resTable[1] <- names(dataStore)

  # plot models of depth ratio 
  pdf(file=paste(workingDir, setNames[currSet], "/", setNames[currSet], rerunFlag,".lmDepthRatio.LPGWS.pdf", sep=""), width = 5, height = 5)
  par(xpd=FALSE)
    for(i in 1:length(gradFit)){
      plot(x = dataStore[[i]][["ratio"]], y = dataStore[[i]][["cn"]], xlab="depth ratio sample vs normal", ylab="chomosome copy state", ylim=c(0,4), xlim=c(0,2), col="grey85", pch=20, main=paste("copy state vs depth ratio ", names(dataStore)[[i]], sep=""))
      abline(gradFit[[i]], col="red")
      mtext(paste("p =", round(summary(gradFit[[i]])[[4]][2,4], digits = 5)))
      mtext(paste("gradient =", round(summary(gradFit[[i]])[[4]][2,1], digits = 2)), line = -1)
      
      resTable[i, "p.value"] <- round(summary(gradFit[[i]])[[4]][2,4], digits = 5)
      resTable[i, "gradient"] <- round(summary(gradFit[[i]])[[4]][2,1], digits = 1)
    }
  dev.off() 
  
  
  # filter results by gradient (must be > 10)
  filterIDList <- names(segmentationList)[which((names(segmentationList) %in% resTable[resTable[["gradient"]] > 0.5 & resTable[["gradient"]] <= 20, "sampleID"] ))]
  #filterIDList <- names(segmentationList)
  filterIDList <- c(filterIDList, bioListEx)
  
  plotList <- plotList[plotList[["ID"]] %in% filterIDList, ]
  colVect <- plotList[["type"]]
  colVect[colVect == "exome"]  <- "black"
  colVect[colVect == "lpwgs"]  <- "grey65"
  plotIndexList <- plotList[["ID"]]
  
  plotListSub <- plotList[plotList[["type"]]!="bulk", ]
  write.csv(x = plotListSub, file = paste0(workingDir, setNames[currSet], "/", setNames[currSet], rerunFlag, "_glandlist.plot.csv"), quote = FALSE)
  
  
  ###########  3) plot CNAs  ##########
  pdf(file=paste(workingDir, setNames[currSet], "/", setNames[currSet], rerunFlag, ".segSummary.pdf", sep=""), height = 0.2*length(plotIndexList), width = 15)
    par(xpd=TRUE, mar=c(4,8,2,2))
    plot(1, 1, col="white", axes=F, xlim=c(0, chromPos[22,"end"]), ylim=c(0,length(plotIndexList)), xlab="", ylab="", main="")
    
    plotCounter <- 1
    for(currSam in 1:length(plotIndexList)){
      currindex <- plotIndexList[currSam]
      corrSegs <- segmentationList[[currindex]]
      
      if(is.null(corrSegs)){
        next
      }
      
      #walk through genome and plot 
      for(currSeg in 1:nrow(corrSegs)){
        if(is.na(corrSegs[currSeg, "A"])){
          next
        }
        currCol <- getCol(corrSegs[currSeg, ])
        rect(xleft = corrSegs[currSeg, "start"], xright = corrSegs[currSeg, "end"], 
             ybottom = plotCounter-0.45, ytop = plotCounter+0.45, col = currCol, border = NA)
      }
      text(x = -10, y = plotCounter, labels = currindex, pos = 2, col = colVect[currSam])
      plotCounter <- plotCounter + 1
    }
    axis(side = 1, at = c(chromPos[1:22,"CENTROSTART"]), labels = chromPos[1:22, "chromosome"], las=2, lwd = 0)
    for(i in 1:nrow(chromPos[1:22,])){
      #chromosome boundries
      lines(x = c(chromPos[i,"end"], chromPos[i,"end"]), y = c(-1,(length(plotIndexList)+1)), lty=1)
      
      #centromeres
      lines(x = c(chromPos[i,"CENTROSTART"], chromPos[i,"CENTROSTART"]), y = c(0,(length(plotIndexList))+1), lty=2, col="grey55")
      lines(x = c(chromPos[i,"CENTROEND"], chromPos[i,"CENTROEND"]), y = c(0,(length(plotIndexList))+1), lty=2, col="grey55")
      
    }
  dev.off()
  
  
}




