# 2.2.1.scoreVariantsVEP-VEST.R
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
# 2.2.0.getVestScores has to be run on all indels before annotations can be included
#
#
################### libraries ###################


###################   begin   ###################

arguments <- commandArgs(trailingOnly = TRUE)

#check number of arguments
#if(length(arguments)!=3){
#  stop("\n#### arguments > 3.4.0.getVariantMetrics.R <sample list file> <platDir> <scriptsOut> ####\n")
#}

#get sample list infomration
sampleList <- read.csv(file="~/mount3/glandList.genomes.csv", header=TRUE, stringsAsFactors=FALSE)
#sampleList <- read.csv(file=arguments[1], header=TRUE, stringsAsFactors=FALSE)
sampleList <- sampleList[sampleList[["retain"]]!=3, ]
setNames <- unique(sampleList[["setID"]])
#setNames <- setNames[c(1,2,3,5,9,10,11)]


vcfDir <- "~/projects/glandSeqProject/3.0.platypusCalls/genomeFinal/"
#vcfDir <- arguments[2]

#vestScores <- read.table(file="~/projects/glandSeqProject/3.1.pathogenicMutations/Variant.Result.tsv", header = TRUE, stringsAsFactors = FALSE, sep="\t", fill = TRUE)
vestScores <- read.csv(file="~/projects/glandSeqProject/3.1.pathogenicMutations/vestAnalysis/Variant.Result.tsv", header = TRUE, stringsAsFactors = FALSE, fill = TRUE, sep="\t")

# make log table to check file contents
logTable <- data.frame(matrix(NA, nrow = length(setNames), ncol = 7))
names(logTable) <- c("setName", "PASSflags", "alleleBiasFlags", "noVarsVCF", "annoVar", "exomeVars", "notes")
logTable["setName"] <- setNames 


for(currSet in 1:length(setNames)){
  
  print(paste0("######### scoring variants for set: ", setNames[currSet], " #########"))
  
  vcfFileName <- paste(vcfDir, setNames[currSet], "/", setNames[currSet], ".somatic.vcf", sep="")
  mainDat <- read.table(file=vcfFileName, header = FALSE)
  
  logTable[currSet, "PASSflags"] <- nrow(mainDat[mainDat[[7]]=="PASS", ])
  logTable[currSet, "alleleBiasFlags"] <- nrow(mainDat[mainDat[[7]]=="alleleBias", ])
  logTable[currSet, "noVarsVCF"] <- nrow(mainDat)
  
  exomeFileName <- paste(vcfDir, setNames[currSet], "/", setNames[currSet], ".somatic.annoVar.exonic_variant_function.tsv", sep="")
  exomeIn <- read.table(file=exomeFileName, sep="\t", stringsAsFactors = FALSE, header = TRUE)
  
  vcfFileName <- paste(vcfDir, setNames[currSet], "/", setNames[currSet], ".somatic.annoVar.variant_function.tsv", sep="")
  mutIn <- read.table(file=vcfFileName, sep="\t", stringsAsFactors = FALSE, header = TRUE)
  
  logTable[currSet, "annoVar"] <- nrow(mutIn)
  
  mutIn <- as.data.frame(append(mutIn, values = NA, after = 2), stringsAsFactors = FALSE)
  names(mutIn)[3] <- "type"
  
  mutIn <- as.data.frame(append(mutIn, values = NA, after = 3), stringsAsFactors = FALSE)
  names(mutIn)[4] <- "AAsub"
  
  mutIn <- as.data.frame(append(mutIn, values = NA, after = 4), stringsAsFactors = FALSE)
  names(mutIn)[5] <- "SIFTsc"
  
  mutIn <- as.data.frame(append(mutIn, values = NA, after = 5), stringsAsFactors = FALSE)
  names(mutIn)[6] <- "SIFTanno"
  
  mutIn <- as.data.frame(append(mutIn, values = NA, after = 6), stringsAsFactors = FALSE)
  names(mutIn)[7] <- "POLYPHENsc"
  
  mutIn <- as.data.frame(append(mutIn, values = NA, after = 7), stringsAsFactors = FALSE)
  names(mutIn)[8] <- "POLYPHENanno"
  
  mutIn <- as.data.frame(append(mutIn, values = NA, after = 8), stringsAsFactors = FALSE)
  names(mutIn)[9] <- "vestFDR"
  
 
  bioList <- read.table(file=paste(vcfDir, setNames[currSet], "/", setNames[currSet], ".bioIDorder.tsv", sep=""), sep = "\t", stringsAsFactors = FALSE, header = TRUE)

  # add VAF columns
  for(i in 1:nrow(bioList)){
    mutIn[paste0(bioList[i,1], ".VAF")] <- mutIn[[paste0(bioList[i,1], ".NV")]] / mutIn[[paste0(bioList[i,1], ".DP")]]
  }
  
  # get VEP scores and format to get pvalues 
  vepFile <- paste(vcfDir, setNames[currSet], "/", setNames[currSet], ".merged.vepAnno.vcf", sep="")
  vepIn <- read.table(file=vepFile, sep="\t", stringsAsFactors = FALSE, header = FALSE)
  names(vepIn) <- c("Uploaded_variation",	"Location",	"Allele",	"Gene",	"Feature",	"Feature_type",	"Consequence",	"cDNA_position",	"CDS_position",	"Protein_position",	"Amino_acids",	"Codons",	"Existing_variation",	"IMPACT",	"DISTANCE",	"STRAND",	"FLAGS",	"VARIANT_CLASS",	"SYMBOL",	"SYMBOL_SOURCE",	"HGNC_ID",	"NEAREST",	"SIFT",	"PolyPhen")
  vepIn["chrom"] <- unlist(strsplit(vepIn[["Location"]], split = ":"))[seq(1,(nrow(vepIn)*2),2)]
  vepIn["pos"] <- unlist(strsplit(vepIn[["Location"]], split = ":"))[seq(2,(nrow(vepIn)*2),2)]
  
  impactOrder <- c("HIGH", "MODERATE", "MODIFIER", "LOW")
  impactOrderSift <- c("deleterious", "deleterious_low_confidence", "tolerated_low_confidence", "tolerated")
  impactOrderPoly <- c("probably_damaging", "possibly_damaging", "benign", "unknown")
  
  vepIn[vepIn == "-"] <- "()"
  vepIn[1] <- unlist(strsplit(vepIn[["SIFT"]], split = "\\("))[seq(1,(nrow(vepIn)*2),2)]
  vepIn[2] <- unlist(strsplit(vepIn[["SIFT"]], split = "\\("))[seq(2,(nrow(vepIn)*2),2)]
  vepIn[2] <- substr(vepIn[[2]], start = 1, stop = nchar(vepIn[[2]])-1)
  vepIn[3] <- unlist(strsplit(vepIn[["PolyPhen"]], split = "\\("))[seq(1,(nrow(vepIn)*2),2)]
  vepIn[4] <- unlist(strsplit(vepIn[["PolyPhen"]], split = "\\("))[seq(2,(nrow(vepIn)*2),2)]
  vepIn[4] <- substr(vepIn[[4]], start = 1, stop = nchar(vepIn[[4]])-1)
  
  #order effects by severity
  vepIn <- vepIn[order(vepIn[["chrom"]], vepIn[["pos"]], match(vepIn[["IMPACT"]], impactOrder), match(vepIn[[1]], impactOrderSift), match(vepIn[[3]], impactOrderPoly)), ]
  
  
  # get exome variants
  #vcfFileName <- paste(vcfDir, setNames[currSet], "/", setNames[currSet], ".somatic.annoVar.exonic_variant_function.tsv", sep="")
  #exomeIn <- read.table(file=vcfFileName, sep="\t", stringsAsFactors = FALSE, header = TRUE)
  
  newTab <- mutIn
  
  maxIter <- nrow(newTab)
  reportIter <- seq(0, maxIter, length.out = 11)
  reportIter <- round(reportIter, digits = 0)
  names(reportIter) <- seq(0,100,10)
  
  for(currMut in 1:nrow(newTab)){
    if(currMut %in% reportIter){
      percCom <- names(reportIter)[reportIter==currMut]
      print(paste("parsing", percCom, "% complete"))
    }
    
    if(newTab[currMut, "region"] %in% c("exonic", "exonic;splicing")){
      currPos <- newTab[currMut, c("chr", "pos")]
      
      tempAnno <- vepIn[vepIn[["chrom"]]==as.character(currPos[1,1]) & vepIn[["pos"]]==as.numeric(currPos[1,2]), ]
      if(nrow(tempAnno)==0){
        #stop()
        tempAnno <- exomeIn[exomeIn[["chr"]]==as.character(currPos[1,1]) & exomeIn[["pos"]]==as.numeric(currPos[1,2]), ]
        if(nrow(tempAnno)!=0){
          newTab[currMut, "type"] <- tempAnno[1, "mutType"]
          animoStr <- tempAnno[1, "anno"]
          animoStr <- strsplit(animoStr, split = ":")[[1]][5]
          if(length(grep(pattern = ",", x = animoStr))!=0){
            animoStr <- strsplit(animoStr, split = ",")[[1]][1]
          }
          newTab[currMut, "AAsub"] <- animoStr
        }else{
          next()
        }
      }else{
        #add information to row
        newTab[currMut, "type"] <- tempAnno[1, "Consequence"]
        
        if(tempAnno[1, "Amino_acids"] != "()"){
          animoStr <- tempAnno[1, "Amino_acids"]
          animoStr <- strsplit(animoStr, split = "\\/")[[1]]
          newTab[currMut, "AAsub"] <- paste(animoStr[1], tempAnno[1, "Protein_position"], animoStr[2], sep="")
          
          newTab[currMut, "SIFTsc"] <- tempAnno[1, 1]
          newTab[currMut, "SIFTanno"] <- tempAnno[1, 2]
          newTab[currMut, "POLYPHENsc"] <- tempAnno[1, 3]
          newTab[currMut, "POLYPHENanno"] <- tempAnno[1, 4]
        }
        
      }
    }else{
      next()
    }
    
  }
  
  
  # copy over VEST p-values
  vestScoresTemp <- vestScores[vestScores[["Sample.ID"]]==setNames[currSet], ]
  for(i in 1:nrow(newTab)){
    if((newTab[i, "ref"] == "-" | newTab[i, "alt"] == "-")){
      rowIndex <- which(vestScoresTemp[["Chromosome"]] == newTab[i, "chr"] & as.numeric(vestScoresTemp[["Position"]]) == as.numeric(newTab[i, "pos"]))
      if(length(rowIndex) == 0){
        next()
      }else{
        tempTab <- vestScoresTemp[rowIndex, ]
        newTab[i, "vestFDR"] <- tempTab[1,  "VEST.FDR"]
      }
    }
  }
  
  
  logTable[currSet, "exomeVars"] <- nrow(newTab[newTab[["region"]]=="exonic", ])
  
  #write annotated sequencing file
  write.table(newTab, file=paste(vcfDir, setNames[currSet], "/", setNames[currSet], ".variant_function.VEPannotated.tsv", sep=""), sep = "\t", quote = FALSE, row.names = FALSE)
  
}


write.csv(logTable, file="~/projects/glandSeqProject/3.0.platypusCalls/genome/2.2.1.logTab.csv", quote = FALSE, row.names = FALSE)



