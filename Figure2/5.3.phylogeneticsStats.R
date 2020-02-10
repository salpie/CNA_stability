# script reads in list of tree files and outputs table of stats

################  notes  ##################
#
# get most parsimonious tree from given sets 
#
##############   funtions   ################

############## main program ################

library(apTreeshape)
library(ape)
#library(phangorn)
#library(irr)

#input sampleList from commandline arguments
arguments <- commandArgs(trailingOnly = TRUE)
if(length(arguments)!=3){
	stop("\n#### please use syntax > Rscript 4.6.phylogeneticStats.R < sample list file > < holding directory > < prepended.name > ####\n")
}
sampleList <- read.csv(file=arguments[1], header=FALSE, stringsAsFactors=FALSE)
holdingDir <- arguments[2]
namePrepended <- arguments[3]

sampleList <- read.csv(file="~/mount3/glandList.exomes.csv", header=TRUE, stringsAsFactors=FALSE)
sampleList <- sampleList[sampleList[["retain"]]==1, ]


holdingDir <- "~/projects/glandSeqProject/5.0.phylogenetics/exomeFinal/"
namePrepended <- ".tre"
namePrepended2 <- ".allTreeSearch.tre"

#get sample names
sampleNames <- unique(sampleList[[1]])
#sampleNames <- sampleNames[c(2,3,5)]
noSets <- length(sampleNames)


######## get the most parsimonious and shape stats ########
for(i in 1:noSets){
  subList <- subset(sampleList, sampleList[1]== sampleNames[i])
  
  treFile <- paste(holdingDir, subList[1,1],"/", subList[1,1], namePrepended, sep="")
  treFileTotal <- paste(holdingDir, subList[1,1],"/", subList[1,1], namePrepended2, sep="")
  
  normalName <- subList[1, "normalID"]
  histFile <- paste(holdingDir, subList[1,1],"/", subList[1,1], ".hist.txt", sep="")
  
  treeList <- as.list(0)
  
  #get specific tree set from file
  if(file.exists(treFileTotal)){
    treeList <- read.nexus(file=treFileTotal)
    #histTab <- read.table(file=histFile, header=FALSE, sep="\t")
  }else if(file.exists(treFile)){
    treeList <- read.nexus(file=treFile)
    histTab <- NA
  }else{
    next
  }
  
  #no of trees
  if(length(is.binary.tree(treeList)) > 1){
    notrees <- length(treeList)
  }else{
    print(paste("#### only one tree for set", sampleNames[i],"####"))
    next()
  }
  
  #setup results table
  resultsTable <- as.data.frame(matrix(NA, nrow=notrees, ncol=3))
  names(resultsTable) <- c("tree length", "colless test (yule)", "colless test (PDA)")
  row.names(resultsTable) <- paste("tree", c(1:notrees))
  
  #for each tree get stats
  for(j in 1:notrees){
    print(paste("calculating stats for tree", j, "of set", subList[1,1])) 
    
    #convert to single tree object
    phyloTree <- treeList[[j]]
    
    #resolve polytomies
    phyloTree <- multi2di(phyloTree, random = TRUE)
    
    #get tree length
    resultsTable[j,1] <- sum(phyloTree$edge.length)
    
    #drop normal sample
    treeTemp <- drop.tip(phyloTree, normalName)
    
    #convert to treeshape object
    phyloShape <- as.treeshape(treeTemp)
    
    #perform stats (silently)
    dummyFile <- file()
    sink(file=dummyFile)
    if(length(phyloShape$names) < 5){
       yuleTemp <- colless.test(phyloShape, model="yule", n.mc=1000, alternative="greater")
    }else{
       yuleTemp <- likelihood.test(phyloShape, model="yule", alternative = "greater")
    }
    pdaTemp <- colless.test(phyloShape, model="pda", n.mc=1000, alternative="greater")
    sink()
    close(dummyFile)
    
    #populate table
    resultsTable[j,2] <- yuleTemp$p.value
    resultsTable[j,3] <- pdaTemp$p.value
  }
  
  #reorder table
  resultsTable <- resultsTable[order(resultsTable[1]), ]
  #parseTreeList[i, treeTypes[eventType]] <- row.names(resultsTable)[1]
  
  write.table(resultsTable, file=paste(holdingDir, subList[1,1], "/", subList[1,1], ".treeStats.txt", sep=""), sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
  
  #output most parsimonous tree
  parseFile <- paste(holdingDir, subList[1,1],"/", subList[1,1], ".parseTree.tre", sep="")
  mostParseTree <- as.numeric(strsplit(rownames(resultsTable)[1], " ")[[1]][2])
  parseTree <- treeList[[mostParseTree]]
  write.tree(parseTree, file=parseFile)
  
  #plot tree distribution file
  #if(length(histTab)>1){
  #  pdf(file=paste(holdingDir, subList[1,1],"/", subList[1,1], ".tree_dist.pdf", sep=""), onefile=TRUE, width=5, height=5)
  #    barplot(histTab[[2]], main="tree length distribution", xlab="tree length", ylab="frequency", xaxt="n", space=0)
  #    lengthInt <- histTab[[1]]
  #    axis(side=1, at=c(0.5:(length(histTab[[1]])-0.5)), labels=lengthInt, las=2, tick=FALSE)
  #  dev.off()
  #}
  
  #plot graphs of distribution
  #pdf(file=paste(holdingDir, subList[1,1],"/", subList[1,1], ".tree_p.value.pdf", sep=""), onefile=TRUE, width=5, height=5)
  #  plot(resultsTable[[2]], resultsTable[[1]], main="tree length vs p-value from yule", xlab="p-value from yule", ylab="tree length")
  #dev.off()
}



