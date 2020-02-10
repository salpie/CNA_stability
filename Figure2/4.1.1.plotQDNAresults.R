# 4.1.1.plotQDNAresults.R
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

#sample lists and details
subSample <- read.csv(file="~/Projects/glandSeqProject/4.1.QDNAseq/glandList.LPWGS.C274.csv", header=TRUE, stringsAsFactors=FALSE)
setID <- subSample[1,"setID"]

#reference material
wholeGenData <- read.table(file = "~/Projects/ReferenceGenome/chromosomeStartEnd.csv", 
                           header = TRUE, sep=",", stringsAsFactors = FALSE)
wholeGenData["gPos"] <- 0
wholeGenData[1, "gPos"] <- wholeGenData[1, "size"]
for(reg in 2:nrow(wholeGenData)){
  wholeGenData[reg, "gPos"] <- wholeGenData[(reg-1), "gPos"] + wholeGenData[reg, "size"]
}
wholeGenData[23, "chromosome"] <- "X"
wholeGenData[24, "chromosome"] <- "Y"


#centromere data
centromereData <- read.table(file = "~/Projects/ReferenceGenome/centromerePositions.csv", 
                           header = TRUE, sep=",", stringsAsFactors = FALSE)
centromereData["CENTROSTART"] <- centromereData[["CENTROSTART"]] + wholeGenData[["start"]]
centromereData["CENTROEND"] <- centromereData[["CENTROEND"]] + wholeGenData[["start"]]

#data dirs
segDir <- "/Users/cross01/projects/glandSeqProject/4.1.QDNAseq/LPWGS/"
segName <- ".segmentations.txt"

sequenzaDir <- "/Users/cross01/projects/glandSeqProject/6.sequenzaCalls/2.mergedSegs/"
seqName <- "_segments.txt"

#QDNAseq reference positions
dataRef <- read.table(file = paste(segDir, setID, "/segTables/", subSample[56, "sampleID"], segName, sep=""), header = TRUE, sep="\t", stringsAsFactors = FALSE)
dataRef <- dataRef[1:4]





########## make main plotting table to populate ##############


#get cumulative positions across the genome


plotTab <- data.frame(matrix(NA, nrow = nrow(subSample), ncol = nrow(dataRef) ))
row.names(plotTab) <- subSample[["sampleID"]]
names(plotTab) <- dataRef[["feature"]]

#get number of tumour regions to plot

sidePieces1 <- subSample[subSample[["sideNo"]]==1, ]
piecesSide1 <- unique(sidePieces1[c(3:4,6)])

sidePieces2 <- subSample[subSample[["sideNo"]]==2, ]
piecesSide2 <- unique(sidePieces2[c(3:4,6)])

totalpieces <- rbind(piecesSide1, piecesSide2)

#populate main plotting table
for(currRow in 1:nrow(plotTab)){
  
  #read in data
  dataTemp <- read.table(file = paste(segDir, setID, "/segTables/", subSample[currRow, "sampleID"], segName, sep=""), 
                         header = TRUE, sep="\t", stringsAsFactors = FALSE)
  plotTab[currRow, ] <- dataTemp[[5]]
 
}

plotTab[plotTab < -1] <- -1
plotTab[plotTab > 1] <- 1

#get sequenza plots (region archetypes)
seqList <- as.list(NA)
for(tumPiece in 1:nrow(totalpieces)){
  if(is.na(totalpieces[tumPiece, "archetype"])){
    seqList[[tumPiece]] <- NA
  }else{
    seqList[[tumPiece]] <- read.table(file=paste(sequenzaDir, setID, "/", totalpieces[tumPiece, "archetype"], "/", totalpieces[tumPiece, "archetype"], seqName, sep=""), header = TRUE, sep="\t", stringsAsFactors = FALSE)
  }
}


############### plot data ##################

#colour table
colourTab <- colorRampPalette(c("blue", "white", "red"), alpha = TRUE)(21)
names(colourTab) <- seq(-1, 1, 0.1)
  
#record piece label position
totalpieces["coord"] <- NA


graphName <- paste("~/Projects/glandSeqProject/4.1.QDNAseq/LPWGS/", subSample[1, "setID"], "/", subSample[1, "setID"], ".CNAsummary.pdf", sep="")
pdf(file=graphName, width = 10, height = 4)
  par(mfrow=c(4,4,6,6))
  plot(1, 1, col="white", axes=F, xlim=c(0, wholeGenData[(nrow(wholeGenData)-2), "gPos"]), ylim=c(0,(nrow(plotTab) + nrow(totalpieces))), xlab="", ylab="", main="")
  
  archCounter <- 1
  plotRowCounter <- 1
  for(currPlot in 1:nrow(plotTab)){
    
    #get current piece
    currArch <- strsplit(rownames(plotTab)[currPlot], "_")[[1]][c(2:3)]
    currArch[1] <- strsplit(currArch[1], "E")[[1]][2]
    currArch[2] <- strsplit(currArch[2], "P")[[1]][2]
    assArch <- which(totalpieces[["sideNo"]]==currArch[1] & totalpieces[["pieceNo"]]==currArch[2] )
    
    ###### plot archetype karyotype for this region (if needed) #####
    if(currPlot == 1 | archCounter != assArch){
      tempArch <- seqList[[archCounter]]
      if(!is.data.frame(tempArch)){
        #no archetype available (add blank line)
        rect(xleft = 0, xright = wholeGenData[(nrow(wholeGenData)-2), "gPos"], ybottom = (plotRowCounter-0.2), ytop = (plotRowCounter+0.2), col = "white", border = FALSE)
      }else{
        #plot architype
        archData <- seqList[[archCounter]]
        
        #plot each segmentation
        for(seg in 1:nrow(archData)){
          #get seg colour
          if(archData[seg, "CNt"] > 2){
            segCol <- "red"
          }else if(archData[seg, "CNt"] < 2){
            segCol <- "blue"
          }else{
            segCol <- "white"
          }
          
          #plot cumulative (absolute genome) position
          currChrom <- strsplit(archData[seg, "chromosome"], split = "chr")[[1]][2]
          cumPos <- wholeGenData[wholeGenData[["chromosome"]]==currChrom, "gPos"] - wholeGenData[wholeGenData[["chromosome"]]==currChrom, "size"]
          startPos <- archData[seg, "start.pos"] + cumPos
          endPos <- archData[seg, "end.pos"] + cumPos
          rect(xleft = startPos, xright = endPos, ybottom = (plotRowCounter-1), ytop = (plotRowCounter+1), col = segCol, border = FALSE)
        }
        # increment plot row counter
        points(x = 0, y=plotRowCounter, pch = 8, cex=0.25)
      }
      #lines(x = c(0, wholeGenData[(nrow(wholeGenData)-2), "gPos"]), y = c((currPlot-1),(currPlot-1)), pty=0.2)
      abline(h=(plotRowCounter-1), lwd=0.3, lty=1)
      plotRowCounter <- plotRowCounter + 2
      totalpieces[assArch, "coord"] <- plotRowCounter
      archCounter <- assArch
    }
    
    ##### plot current gland segmentations ####
    for(seg in 1:ncol(plotTab)){
      #get seg colour
      plotRatio <- plotTab[currPlot, seg]
      ratioDist <- abs(as.numeric(names(colourTab)) - plotRatio)
      segCol <- colourTab[which(min(ratioDist) == ratioDist)]
      
      #plot cumulative (absolute genome) position
      currChrom <- dataRef[seg, "chromosome"]
      cumPos <- wholeGenData[wholeGenData[["chromosome"]]==currChrom, "gPos"] - wholeGenData[wholeGenData[["chromosome"]]==currChrom, "size"]
      startPos <- dataRef[seg, "start"] + cumPos
      endPos <- dataRef[seg, "end"] + cumPos
      rect(xleft = startPos, xright = endPos, ybottom = (plotRowCounter-0.2), ytop = (plotRowCounter+0.2), col = segCol, border = FALSE)
    }
    # increment plot row counter
    plotRowCounter <- plotRowCounter + 1
  }
  #plot axis and labels
  
  #add chromosome labels
  for(chr in 1:(nrow(wholeGenData)-2)){
    if(chr > 18){
      sizeTemp <- 5
    }else{
      sizeTemp <- 10
    }
    
    labPos <- (wholeGenData[chr, "gPos"] - (wholeGenData[chr, "size"] / 2))
    par(ps=sizeTemp)
    axis(side=1, at = labPos, labels = wholeGenData[chr, "chromosome"], tick = FALSE)
  }
  
  # add chromosome boundries
  abline(v=wholeGenData[1:22,"start"], lwd=0.1, lty=1, col="grey60")
  abline(v=centromereData[1:22,"CENTROSTART"], lwd=0.2, lty=1, col="grey60")
  abline(v=centromereData[1:22,"CENTROEND"], lwd=0.2, lty=1, col="grey60")
  
  
  #add piece information
  par(ps=6)
  #totalpieces[8, "coord"] <- 100
  for(pieceInfo in 1:nrow(totalpieces)){
    tempInfo <- totalpieces[pieceInfo, "coord"]
    axis(side=4, at = tempInfo, line = -1.8, labels = paste("side:", totalpieces[pieceInfo, "sideNo"], " piece:", totalpieces[pieceInfo, "pieceNo"], sep=""), tick = FALSE, las=2)
  }
  
dev.off()  

