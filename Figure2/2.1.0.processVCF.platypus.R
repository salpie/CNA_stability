# 2.1.0.processVCFplatypus.R
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
#  2.2.1.scoreVariantsVEP-VEST.R  <- 2.2.0.getVestScores.R
#         |
#         V
#  2.3.0.getSummaryAndPlot.R
#         |
#         V
#  2.4.0.filterVariants.R
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

#samples <- subSample; prepend <- ".somatic.tsv"; output <- ".somatic"; dirSub <- platDir
makeAnnovar <- function(samples, prepend, output, dirSub){
  sampleNames <- unique(samples[[1]])
  namePrepended <- prepend
  
  #now process .vcf file and make annoVar input 
  for(j in 1:length(sampleNames)){
    print(paste("#### making table for sample ", sampleNames[j], " ####",sep=""))
    
    #subset main list
    subSample <- subset(samples, samples[1]==sampleNames[j])
    #print(subSample)
    
    setName <- unique(subSample[[1]])
    
    #setup input/output names
    dataIn <- read.table(file=paste(dirSub, setName,"/", setName, namePrepended, sep=""), sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
    
    dataIn <- as.data.frame(append(dataIn, list(to = dataIn[[2]]), after = 2), stringsAsFactors = FALSE)
    names(dataIn)[1:6] <- c("chr", "pos", "end", "ref", "alt", "filter")
    
    # correct for indel format 
    indelVect <- c()
    counter <- 0
    unknownVect <- c()
    unknownCounter <- 0
    normalBases <- c("A", "C", "G", "T")
    for(row in 1:nrow(dataIn)){
      if(dataIn[row, 4] %in% normalBases & dataIn[row, 5] %in% normalBases){
        next()
      }else{
        # this is an indel, reformat
        tempPos <- dataIn[row, "pos"] 
        tempRef <- as.character(dataIn[row, "ref"])
        tempAlt <- as.character(dataIn[row, "alt"])
        if(nchar(tempRef)!=1 & nchar(tempAlt)==1){
          mutsize <- (nchar(tempRef))-1
          
          #variant is a deletion
          dataIn[row, "ref"] <- substr(tempRef, start = 2, stop = nchar(tempRef))
          dataIn[row, "alt"] <- "-"
          dataIn[row, "pos"] <- tempPos + 1
          dataIn[row, "end"] <- tempPos + mutsize
          
        }else if(nchar(tempRef)==1 & nchar(tempAlt)!=1){
          #variant is a insertion
          mutsize <- nchar(tempAlt)
          
          #variant is a deletion
          dataIn[row, "ref"] <- "-"
          dataIn[row, "alt"] <- substr(tempAlt, start = 2, stop = mutsize)
          dataIn[row, "end"] <- tempPos
        }else{
          unknownCounter <- unknownCounter + 1
          unknownVect[unknownCounter] <- row
        }
        counter <- counter + 1
        indelVect[counter] <- row
      }
    }
    
    unknownIndels <- dataIn[unknownVect, ]
    unknownFile <- paste(dirSub, setName,"/", setName, output, ".unknownIndelTypes.tsv", sep="")
    write.table(unknownIndels, file= unknownFile, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
    
    dataIn <- dataIn[-unknownVect, ]
    annoInput <- paste(dirSub, setName,"/", setName, output, ".annoVarInput.tsv", sep="")
    write.table(dataIn, file= annoInput, sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)
    
    annoOut <- paste(dirSub, setName,"/", setName, output,".annoVar", sep="") 
    
    #run annovar program on new file
    system(command=paste("~/bin/annovar/annotate_variation.pl -out " ,annoOut ," -build hg19 ", annoInput, " ~/bin/annovar/humandb/",sep=""))
    
    #annovar output files
    annoExoFunct <- paste(dirSub, setName,"/", setName, output,".annoVar.exonic_variant_function", sep="")
    annoExoFunctNew <- paste(dirSub, setName,"/", setName, output,".annoVar.exonic_variant_function.tsv", sep="")
    
    annoVarFunct <- paste(dirSub, setName,"/", setName, output,".annoVar.variant_function", sep="")
    annoVarFunctNew <- paste(dirSub, setName,"/", setName, output,".annoVar.variant_function.tsv", sep="")
    
    
    #rename files with .tsv prepend
    system(command=paste("mv ", annoExoFunct, " ", annoExoFunctNew, sep=""))
    system(command=paste("mv ", annoVarFunct, " ", annoVarFunctNew, sep=""))
    
    system(command=paste("rm ", annoInput, sep=""))
    
    # read in and name headers
    bioNamesSub <- read.table(file=paste(dirSub, setName, "/", setName, ".bioIDorder.tsv", sep=""), sep = "\t", header = TRUE)
    
    dataIn <- read.table(file = annoExoFunctNew, sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
    names(dataIn) <- c("line", "mutType", "anno", "chr", "pos", "end", "ref", "alt", "filter", paste0(bioNamesSub[[1]], ".DP"), paste0(bioNamesSub[[1]], ".NV"))
      
    write.table(dataIn, file= annoExoFunctNew, sep="\t", quote=FALSE, row.names=FALSE)
    
    
    dataIn <- read.table(file = annoVarFunctNew, sep="\t", header=FALSE, fill=TRUE, stringsAsFactors = FALSE)
    names(dataIn) <- c("region", "gene", "chr", "pos", "end", "ref", "alt", "filter", paste0(bioNamesSub[[1]], ".DP"), paste0(bioNamesSub[[1]], ".NV"))
    
    write.table(dataIn, file= annoVarFunctNew, sep="\t", quote=FALSE, row.names=FALSE)
    
  }
}



###################   begin   ###################

#get arguments from script
arguments <- commandArgs(trailingOnly = TRUE)

#check number of arguments
if(length(arguments)!=3){
  stop("\n#### arguments > 2.1.0.processVCFplatypus.R <sample list file> <platDir> <outAnno> ####\n")
}

#get sample list infomration
# sampleList <- read.csv(file="~/mount3/glandList.genomes.csv", header=TRUE, stringsAsFactors=FALSE)
sampleList <- read.csv(file=arguments[1], header=TRUE, stringsAsFactors=FALSE)
sampleList <- sampleList[sampleList[["retain"]]!=3,]
setNames <- unique(sampleList[[1]])

# select final sets
#setNames <- setNames[c(1:9)]

platDir <- "~/projects/glandSeqProject/3.0.platypusCalls/genomeFinal/"
platDir <- arguments[2]

vcfName <- ".merged.vcf"

#outAnno <- "~/Projects/glandSeqProject/4.1.0.annoPlatypusCalls/exomes/"
outAnno <- "~/projects/glandSeqProject/3.0.platypusCalls/genomeFinal/"
outAnno <- arguments[3]

#concatenate names in table and delete unwanted columns
sampleList["sampleInfo"] <- paste(platDir, sampleList[[1]], "/", sampleList[[1]], vcfName, sep="")


#perform pipeline for each sample
for(j in 1:length(setNames)){
  
  currSetId <- setNames[j]
	print(paste("#### filtering sample ", currSetId, " ####",sep=""))
	
	#setup input/output names
  subSample <- sampleList[sampleList[[1]]==currSetId, ]
  
	#### these file are temporary and later deleted ####
	FILTName <- paste(platDir, currSetId,"/", currSetId, ".FILT.vcf", sep="")
	VARName <- paste(platDir, currSetId,"/", currSetId, ".VAR.vcf", sep="")
	SOMAName <- paste(platDir, currSetId,"/", currSetId, ".SOMA.vcf", sep="")
	
	#somatic files
	confTotalName <- paste(platDir, currSetId,"/", currSetId, ".somatic.vcf", sep="")
	confTotalOutput <- paste(platDir, currSetId,"/", currSetId, ".somatic.tsv", sep="")
	
	#SNP (germline variant) files
	germTotalName <- paste(platDir, currSetId,"/", currSetId, ".germline.vcf", sep="")
	germTotalOutput <- paste(platDir, currSetId,"/", currSetId, ".germline.tsv", sep="")
	
	#biopsy order lists and therefore indexes to remove
	biopList <- system(command = paste("grep \"#CHROM\" ", subSample[1,"sampleInfo"], sep=""), intern = TRUE, wait = TRUE)
	biopList <- strsplit(biopList, split = "\t")
	biopList <- biopList[[1]]
	biopList <- biopList[10:length(biopList)]
	write.table(biopList, file=paste(platDir, currSetId, "/", currSetId, ".bioIDorder.tsv", sep=""), sep = "\t", row.names = FALSE, quote = FALSE)
	normalIndex <- (which(biopList==subSample[1, "normalID"]))-1
	
	#subset exclusion list (set to never remove sample)
	colIndexes <- c(0:(length(biopList)-1))
	remList <- NULL
	if(length(remList)==0){
	  keepList <- colIndexes
	}else{
	  keepList <- colIndexes[-(which(biopList %in% remList))]
	  biopListTemp <- biopList[-(which(biopList %in% remList))]
	}
	
	#remove normal column from keep list
	keepListNoNorm <- keepList[-which(keepList == normalIndex)]

	#prepare indexes for total variant output
	counterTemp <- 1
	totalIndexStrings <- as.list(NA)
	for(k in keepListNoNorm){
		totalIndexStrings[[counterTemp]] <- paste(" ((GEN[", k,"].NR > 9) & (GEN[", k,"].NV > 0)) | ", sep="")
		counterTemp <- counterTemp +1
	}
	totalIndexStrings[[length(totalIndexStrings)]] <- substr(totalIndexStrings[[length(totalIndexStrings)]], 1,  (nchar(totalIndexStrings[[length(totalIndexStrings)]]) - 2 ))
	
	
	#prepare indexes for extractFields command
	counterTemp <- 1
	extractStringsNR <- as.list(NA)
	for(k in keepList){
		extractStringsNR[[counterTemp]] <- paste(" \"GEN[", k,"].NR\" ", sep="")
		counterTemp <- counterTemp +1

	}
	counterTemp <- 1
	extractStringsNV <- as.list(NA)
	for(k in keepList){
		extractStringsNV[[counterTemp]] <- paste(" \"GEN[", k,"].NV\" ", sep="")
		counterTemp <- counterTemp +1
	}
	

	# 1 .filter by FILTER field
	filtVarCommand <- paste("cat ", subSample[1, "sampleInfo"], " | java -jar ~/bin/SnpSift.jar filter \"( ( (FILTER ='PASS') | (FILTER ='alleleBias') ) )\" > ", FILTName, sep="")
	system(command=filtVarCommand)
	
	# 2. annotate variant types
	annoCommand <- paste("java -jar ~/bin/SnpSift.jar varType ", FILTName ," > ", VARName, sep="")
	system(command=annoCommand)
  
  
	#### somatic files ####
	
	# 3. filter for somatic mutations (not in normal)
  somaticCommand <- paste("cat ", VARName, " | java -jar ~/bin/SnpSift.jar filter \"( ( (GEN[", normalIndex,"].NV < 2) & (GEN[", normalIndex,"].NR < 100) & (GEN[", normalIndex,"].NR > 9) ) | ( (GEN[", normalIndex,"].NV < 3) & (GEN[", normalIndex,"].NR > 99) ) )\" > ", SOMAName, sep="")
	system(command=somaticCommand)
	
	
	# 4. filter by read depth (>9X in ANY samples), this is the somatic variants vcf
	depthCommand5 <- paste("cat ", SOMAName, " | java -jar ~/bin/SnpSift.jar filter \"( ( ", paste(totalIndexStrings, collapse=" ") ,"))\" > ", confTotalName, sep="")
	system(command=depthCommand5)
	

	#### germline files ####

	# 4b. produce non-fitered germline file, this is the germline vcf
	depthCommand3 <- paste("cat ", VARName, " | java -jar ~/bin/SnpSift.jar filter \"( (GEN[", normalIndex,"].NV > 0) & (GEN[", normalIndex,"].NR > 9) )\" > ", germTotalName, sep="")
	system(command=depthCommand3)
	
	
	
	#### vcf to text file conversion ####
	
	extractCommand2 <- paste("java -jar ~/bin/SnpSift.jar extractFields ", confTotalName, " \"CHROM\" \"POS\" \"REF\" \"ALT\" \"FILTER\" ", paste(extractStringsNR, collapse=" "), " ", paste(extractStringsNV, collapse=" "), " > ", confTotalOutput, sep="")
	system(command=extractCommand2)
		
	extractCommand5 <- paste("java -jar ~/bin/SnpSift.jar extractFields ", germTotalName, " \"CHROM\" \"POS\" \"REF\" \"ALT\" \"FILTER\" ", paste(extractStringsNR, collapse=" "), " ", paste(extractStringsNV, collapse=" "), " > ", germTotalOutput, sep="")
	system(command=extractCommand5)

	#tidy up
	system(command=paste("rm ", FILTName, sep=""))
	system(command=paste("rm ", VARName, sep=""))
	system(command=paste("rm ", SOMAName, sep=""))

	
  #annotate somatic file with annoVar
  makeAnnovar(subSample, ".somatic.tsv", ".somatic", platDir)
  makeAnnovar(subSample, ".germline.tsv", ".germline", platDir)
	
  system(command=paste("rm ", confTotalOutput, sep=""))
  system(command=paste("rm ", germTotalOutput, sep=""))
  
}



