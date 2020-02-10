# 2.2.0.getVestScores.R
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
# this script makes a single "VEST formatted" table
#
#
################### libraries ###################

###### begin ########
#get arguments from script
arguments <- commandArgs(trailingOnly = TRUE)

#check number of arguments
if(length(arguments)!=4){
  stop("\n#### arguments > 4.4.compileDrivers.R <sample list file> <SNVDir> <indelDir> <outDir> ####\n")
}

#get sample list infomration
sampleList <- read.csv(file="~/mount3/glandList.exomes.csv", header=TRUE, stringsAsFactors=FALSE)
sampleList <- read.csv(file=arguments[1], header=TRUE, stringsAsFactors=FALSE)
sampleList <- sampleList[sampleList[["retain"]]!=3,]
setNames <- unique(sampleList[[1]])
#setNames <- setNames[c(1,2,3,4,9,10,11)]

SNVDir <- "~/projects/glandSeqProject/3.0.platypusCalls/exomes/"
SNVDir <- arguments[2]

outDir <- "~/projects/glandSeqProject/3.1.pathogenicMutations/"
outDir <- arguments[3]

driverInfo <- read.csv(file="~/projects/glandSeqProject/COSMIC_Census_04_2018_allEpithelial_Tier1.csv", stringsAsFactors = FALSE, header = TRUE)

mainindelTable <- data.frame(matrix(NA, ncol = 9, nrow = 0))
names(mainindelTable) <- c("line","mutType","anno","chr","pos","end","ref","alt", "setID")
  
#perform pipeline for each sample
for(j in 1:length(setNames)){
  
  ######### 1. subset indels from main mutation table ###########
  currSetId <- setNames[j]
  
  #### SNV total file ####
  fileName <- paste(SNVDir, currSetId, "/", currSetId, ".somatic.annoVar.exonic_variant_function.tsv", sep="")
  annoData <- read.table(file = fileName, sep = "\t", header = TRUE, stringsAsFactors = FALSE)
  
  annoData <- annoData[1:9]
  annoData[9] <- currSetId
  names(annoData)[9] <- "setID" 
  
  annoData <- annoData[annoData[["mutType"]]=="frameshift deletion" | 
                         annoData[["mutType"]]=="frameshift insertion" |
                         annoData[["mutType"]]=="nonframeshift deletion" | 
                         annoData[["mutType"]]=="nonframeshift insertion", ]
  
  mainindelTable <- rbind(mainindelTable, annoData)
  
}  

# reformat and write out table
mainindelTable <- mainindelTable[c("line", "chr", "pos", "end", "ref", "alt", "setID")]
mainindelTable[1] <- paste0("ID_", c(1:nrow(mainindelTable)))
mainindelTable[4] <- "+"

#write driver table
write.table(mainindelTable, file="~/projects/glandSeqProject/3.1.pathogenicMutations/VESTinput.indel.tsv", sep="\t", quote = FALSE, row.names = FALSE, col.names = FALSE)



