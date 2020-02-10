# 4.4.phylogeneticsPrep.R
# part of the variant calling best practices
# produces a phylogenetic tree construction files (for both SNVs and indels)
#
###### notes #######
#
# example: 
# Rscript ~/Code/bestPractices/4.4.phylogeneticsPrep.R ~/Projects/glandSeqProject/sampleList.csv ~/Projects/glandSeqProject/1.platypusCalls/ ~/Code/bestPractices/4.3.driverMutations.table.csv
#
# note requires fastq files to be organised into read1 and read2 dirs
#
#
#
###### libraries #####


###### begin ########
#get arguments from script
arguments <- commandArgs(trailingOnly = TRUE)

#check number of arguments
if(length(arguments)!=5){
  stop("\n#### arguments > 4.4.phylogeneticsPrep.R <sample list file> <platDir> <phyloDir> <driverFile> <namePrep> ####\n")
}

#get sample list infomration
sampleList <- read.csv(file="~/mount3/7.0.finalFigures/tables/glandList.exomes.csv", header=TRUE, stringsAsFactors=FALSE)
#sampleList <- read.csv(file=arguments[1], header=TRUE, stringsAsFactors=FALSE)
sampleList <- sampleList[sampleList[["retain"]]!=3, ]
sampleList <- sampleList[sampleList[["retain"]]!=0, ]

setNames <- unique(sampleList[["setID"]])
#setNames <- setNames[c(1,2,3,4,5,6,7,8)]


platDir <- "~/projects/glandSeqProject/3.0.platypusCalls/exomesFinal/"
platDir <- arguments[2]

phyloDir <- "~/projects/glandSeqProject/5.0.phylogenetics/exomeFinal/"
phyloDir <- arguments[2]

apocritaDir <- "/data/BCI-EvoCa-SG/5.0.phylogenetics/exomeFinal/"

namePrep <- ".variant_function.VEPannotated.final.tsv"
#namePrep <- ".indel.annoVar.txt"
namePrep <- arguments[5]


#perform pipeline for each sample
for(j in 1:length(setNames)){
  currSetId <- setNames[j]
	print(paste("#### analysing sample ", currSetId, " ####",sep=""))
	
	system(command = paste("mkdir ", phyloDir, currSetId, "/", sep=""))
	
	#setup input/output names
  subSample <- sampleList[sampleList[[1]]==currSetId, ]
  normID <- subSample[1, "normalID"]
  
  biopList <- subSample[subSample[["retain"]]==1, "sampleID"]
  normName <- subSample[1, "normalID"]
  bioTotal <- biopList[-(which(biopList == normName))]
  noSamples <- length(bioTotal)
  
  annoTable <- read.table(file=paste(platDir, currSetId, "/", currSetId, namePrep, sep=""), sep="\t", header = TRUE, stringsAsFactors = FALSE)
  
  # remove indels
  annoTable <- annoTable[annoTable[["ref"]] %in% c("A", "C", "G", "T") & annoTable[["alt"]] %in% c("A", "C", "G", "T"), ]
  
  
  #mark each variant by phylogenetic locations and VAF
  annoTable[ncol(annoTable)+1] <- 0
  names(annoTable)[ncol(annoTable)] <- "phyloLoc"
  
  annoTable[ncol(annoTable)+1] <- 0
  names(annoTable)[ncol(annoTable)] <- "VAForder"
  
  for(currRow in 1:nrow(annoTable)){
    #get phylogenetic location
    assessRow <- annoTable[currRow, paste(bioTotal, ".VAF", sep="")] > 0
    assessTable <- table(assessRow)
    
    if(!(TRUE %in% names(assessTable))){
      annoTable[currRow, "phyloLoc"] <- 0
      next()
    }
    
    #add VAF information (sum)
    tempVal <- annoTable[currRow, paste(bioTotal, ".VAF", sep="")]
    tempVal <- tempVal[!is.na(tempVal)]
    tempVal <- tempVal[tempVal!=0]
    annoTable[currRow, "VAForder"] <- sum(tempVal)
  
    #if trunkal mark as 1 etc
    if(as.integer(assessTable["TRUE"])==noSamples){
      annoTable[currRow, "phyloLoc"] <- 1
    }
    if(as.integer(assessTable["TRUE"])==1){
      annoTable[currRow, "phyloLoc"] <- 3
    }
    if(as.integer(assessTable["TRUE"])!=1 & as.integer(assessTable["TRUE"])!=noSamples){
      annoTable[currRow, "phyloLoc"] <- 2
    }
  }
  
  #annoTable <- annoTable[annoTable[["phyloLoc"]]!=0, ]
  sortTable <- annoTable[order( annoTable[["phyloLoc"]], -annoTable[["VAForder"]] ), ]
  
  #output sorted table
  write.table(sortTable, file=paste(phyloDir, currSetId, "/", currSetId, ".phyloLoc.sorted.tsv", sep=""), quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  
  
  ########## prepare data for phylogenetic tree construction
  
  #convert VAFs to binary table
  confData <- sortTable[, paste(biopList, ".VAF", sep="")]
  confData[confData > 0] <- 1

  # set normal genotypes to all zero (removing very small number of non-ref reads)
  confData[paste0(normID, ".VAF")] <- 0
  
  #store strings with names in list
  sampleData <- as.list(0)
  for(z in 1:ncol(confData)){	
    sampleData[[z]] <- confData[[z]]
    names(sampleData)[[z]] <- strsplit(names(confData)[z], split = "[.]")[[1]][1]
  }
  noSamples <- noSamples + 1
  counter<-1
  asStrings <- as.list(0)
  for(k in seq(1,(2*(noSamples)),2)){
    asStrings[[k]] <- names(sampleData[counter])
    asStrings[[k+1]] <- paste(sampleData[[counter]], collapse='')
    counter <- counter+1
  }
  
  #get assign all file locations and parsomony parameters
  
  noTaxa <- noSamples
  taxaNames <- paste0(unlist(strsplit(names(confData), split = "[.]"))[seq(1, (2*ncol(confData)), 2)], collapse = " ")
  noCharacters <- nchar(asStrings[[2]][[1]])
  matrixStrings <- paste(asStrings, collapse="\n")
  ascestState <- 0
  for(w in 1:(noCharacters-1)){ascestState <- paste(ascestState,0,sep="")}
  
  upperHIvalue <- round(noCharacters + ((noCharacters / 100) * 55), 0)
  
  #dir names for nexux control file
  outputLoc <- paste(apocritaDir, currSetId, "/", currSetId, sep="")
  
  treeFile <- paste(outputLoc, ".tre", sep="")
  treeFileAll <- paste(outputLoc, ".allTreeSearch.tre", sep="")
  logFile <- paste(outputLoc, ".log", sep="")
  logFileAll <- paste(outputLoc, ".allTrees.log", sep="")
  logFileBoot <- paste(outputLoc, ".boot.log", sep="")
  histogramFile <- paste(outputLoc, ".hist.txt", sep="")
  
  nexusTreeFileLoc <- paste(outputLoc, ".nex", sep="")
  nexusBootFileLoc <- paste(outputLoc, ".boot.nex", sep="")
  nexusAllTreesFile <- paste(outputLoc, ".allTrees.nex", sep="")
  
  outputLocNex <- paste(phyloDir, currSetId, "/", currSetId, sep="")
  
  #put all strings into nexus format for tree production
  totalStringsTree <- paste("#NEXUS
begin paup;
set criterion=parsimony autoclose=yes warntree=no warnreset=no monitor=no warntsave=no maxtrees=1000 increase=no;
log start file= ", logFile,"  replace;

begin taxa;
dimensions ntax=", noTaxa,";
taxlabels
", taxaNames,";
end;

begin characters;
dimensions nchar= ", noCharacters," ;
format symbols = \"01\";
matrix\n", matrixStrings,";

end;

begin assumptions;
options deftype = unord;

usertype statetransitionsmatrix (stepmatrix)= 2
0 1 
. i 
1 . 
;

end;

hsearch addseq=simple nreps=20 hold=1 rearrlimit=10000000 nbest=1000;

outgroup ", normName," /only;

roottrees outroot=monophyl rootmethod=outgroup;

describetrees /plot=both root=outgroup brlens=yes chglist=yes apolist=yes diag=yes homoplasy=yes xout=both;

savetrees from=1 to=1000 file= ", treeFile,"  savebootp=Both format=altnex brlens=yes root=yes append=yes;

log stop;

end;

quit;", sep="")


  #write nexus file
  lapply(totalStringsTree, write, paste(outputLocNex,".nex",sep=""), append=FALSE)
  
  
  #put all strings into nexus format for all tree production
  totalStringsTreeAll <- paste("#NEXUS
begin paup;
set criterion=parsimony autoclose=yes warntree=no warnreset=no monitor=no warntsave=no maxtrees=1000 increase=no;
log start file= ", logFileAll,"  replace;

begin taxa;
dimensions ntax=", noTaxa,";
taxlabels
", taxaNames,";
end;

begin characters;
dimensions nchar= ", noCharacters," ;
format symbols = \"01\";
matrix\n", matrixStrings,";

end;

begin assumptions;
options deftype = unord;

usertype statetransitionsmatrix (stepmatrix)= 2
0 1 
. i 
1 . 
;

end;

alltrees fdfile=", histogramFile," keep=", upperHIvalue,";

outgroup ", normName," /only;

roottrees outroot=monophyl rootmethod=outgroup;

describetrees /plot=both root=outgroup brlens=yes chglist=yes apolist=yes diag=yes homoplasy=yes xout=both;

savetrees from=1 to=1000 file= ", treeFileAll,"  savebootp=Both format=altnex brlens=yes root=yes append=yes;

log stop;

end;

quit;", sep="")


  #write nexus file
  lapply(totalStringsTreeAll, write, paste(outputLocNex,".allTrees.nex",sep=""), append=FALSE)
  
  
  totalStringsBoot <- paste("#NEXUS
begin paup;
set criterion=parsimony autoclose=yes warntree=no warnreset=no monitor=no warntsave=no maxtrees=1000 increase=no;
log start file= ", logFileBoot,"  replace;

begin taxa;
dimensions ntax=", noTaxa,";
taxlabels
", taxaNames,";
end;

begin characters;
dimensions nchar= ", noCharacters," ;
format symbols = \"01\";
matrix\n",matrixStrings,";

end;

begin assumptions;
options deftype = unord;

usertype statetransitionsmatrix (stepmatrix)= 2
0 1 
. i 
1 . 
;

end;

hsearch addseq=simple nreps=20 hold=1 rearrlimit=10000000;

outgroup ", normName," /only;

roottrees outroot=monophyl rootmethod=outgroup;

bootstrap nreps=10000;

log stop;

end;

quit;", sep="")

	
	#write nexus file
	lapply(totalStringsBoot, write, paste(outputLocNex,".boot.nex",sep=""), append=FALSE)

	#create shell scripts to run nexus files
	shellStrings <- paste("#!/bin/sh
#$ -cwd
#$ -V
#$ -pe smp 1            # Request 2 CPU cores
#$ -l h_rt=48:0:0      # Request 48 hour runtime
#$ -l h_vmem=2G    # Request 12GB RAM / core, i.e. 24GB total

./bin/paup4a166_centos64 -n ", nexusTreeFileLoc, "
./bin/paup4a166_centos64 -n ",nexusBootFileLoc, "
./bin/paup4a166_centos64 -n ", nexusAllTreesFile
, sep="")

	#write nexus file
	lapply(shellStrings, write, paste(phyloDir, currSetId,".runPAUP.sh",sep=""), append=FALSE)
	
	
}

