# 1.0.0.makeMutectScripts.R
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
#  2.2.0.scoreVariantsVEP-VEST.R (local script: annotate variant using annovar and add VEP information)
#         |
#         V
#  2.3.0.filterVariants.R
#         |
#         V
#  2.4.0.getSummaryAndPlot.R
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

###################   begin   ###################

arguments <- commandArgs(trailingOnly = TRUE)

#check number of arguments
if(length(arguments)!=3){
  stop("\n#### arguments > 1.0.0.makeMutectScripts.R <sample list file> <platDir> <scriptsOut> ####\n")
}

#get sample list infomration
sampleList <- read.csv(file="~/projects/glandSeqProject/glandList.exomes.csv", header=TRUE, stringsAsFactors=FALSE)
sampleList <- read.csv(file=arguments[1], header=TRUE, stringsAsFactors=FALSE)

#retain bulks
sampleList <- sampleList[sampleList[["retain"]]==2, ]
setNames <- unique(sampleList[["setID"]])

vcfDir <- "/data/BCI-EvoCa-SG/3.0.mutectCalls/exomes/"
vcfDir <- arguments[2]

bamDir <- "/data/BCI-EvoCa-SG/1.3.bamFiles/exomes/"
bamDir <-  arguments[3]

scriptsOut <- "~/projects/glandSeqProject/A.runScripts/3.0.runMutect/"
scriptsOut <- arguments[4]

mutectJar <- "/data/home/mpx155/bin/gatk-4.0.1.0/gatk-package-4.0.1.0-local.jar"

system(command = paste("mkdir ", scriptsOut, sep=""))

vcfNormals <- c()
for(currSet in 1:length(setNames)){
  
  subSample <- sampleList[sampleList[["setID"]]==setNames[currSet], ]
  
  #output bam file list
  biopsyNames <- subSample[["sampleID"]]
  normalID <- subSample[1,"normalID"]
  biopsyNames <- biopsyNames[biopsyNames!=normalID]
  
  normalBam <- paste(bamDir, setNames[currSet], "/", normalID, "/", normalID, ".clipped.bam", sep="")
  vcfOutPON <- paste(vcfDir, setNames[currSet], "/", normalID, ".normalArtifacts.vcf", sep="")
  vcfNormals[currSet] <- vcfOutPON
  
  vcfFiles <- c()
  vcfCounter <- 1
  
  system(paste("mkdir ", scriptsOut, setNames[currSet], sep=""))
  for(currBio in 1:length(biopsyNames)){
    currID <- biopsyNames[currBio]
    bamName <- paste(bamDir, setNames[currSet], "/", biopsyNames[currBio], "/", biopsyNames[currBio], ".clipped.bam", sep="")
    vcfOut <- paste(vcfDir, setNames[currSet], "/", biopsyNames[currBio], ".mutect.vcf", sep="")
    
    vcfOutMark <- paste(vcfDir, setNames[currSet], "/", biopsyNames[currBio], ".mutect.mark.vcf", sep="")
    vcfOutFilt <- paste(vcfDir, setNames[currSet], "/", biopsyNames[currBio], ".mutect.filt.vcf", sep="")
    
    
    vcfFiles[vcfCounter] <- vcfOutFilt
    vcfCounter <- vcfCounter + 1
    
    if(subSample[1, "regions"] == "WGS"){
      regionsString <- ""
    }else{
      regionsString <- paste("-L ", subSample[1, "regions"], ".bed", sep="")
    }
    
    outName <- paste(scriptsOut, setNames[currSet], "/runMutect.", setNames[currSet], ".", biopsyNames[currBio], ".sh", sep="") 
    totalStrings <- paste("#!/bin/sh
#!/bin/sh
#$ -cwd
#$ -pe smp 10
#$ -l h_rt=240:0:0
#$ -l h_vmem=6G

mkdir ", vcfDir, setNames[currSet],"/

#index bam file
samtools index ", bamName,"


#run mutect 2 (with PON)
java -Xmx4g -jar ", mutectJar," Mutect2 \\
-R ", subSample[currBio, "reference"], " \\
-I ", bamName," \\
-I ", normalBam," \\
-tumor ", biopsyNames[currBio]," \\
-normal ", normalID," \\
-pon /data/BCI-EvoCa-SG/3.0.mutectCalls/exomes/PON.vcf.gz \\
--dont-use-soft-clipped-bases TRUE \\
-O ", vcfOut," \\
", regionsString,"

#use GATK to filter variants
java -Xmx4g -jar ", mutectJar," FilterMutectCalls \\
-V ", vcfOut," \\
-O ", vcfOutMark," \\
--normal-artifact-lod 0.3

#use GATK to filter vcf file by \"PASS\" flag
java -Xmx4g -jar /data/home/mpx155/bin/gatk-4.0.1.0/gatk-package-4.0.1.0-local.jar SelectVariants \\
-V ", vcfOutMark," \\
-O ", vcfOutFilt," \\
-select 'vc.isNotFiltered()'

", sep="")
    lapply(totalStrings, write, outName, append=FALSE)
  }
 
  
  
  #make vcf merging script
  mergedVCF <- paste(vcfDir, setNames[currSet], "/", setNames[currSet], ".mutect.vcf", sep="")
  annoVCF <- paste(vcfDir, setNames[currSet], "/", setNames[currSet], ".mutect.vepAnno.txt", sep="")
  
  outNameMerge <- paste(scriptsOut, "mergeAndAnnoMutectVCFs.", setNames[currSet], ".sh", sep="") 
  totalStringsMerge <- paste("#!/bin/sh
#$ -cwd
#$ -V
#$ -pe smp 1
#$ -l h_rt=10:0:0      # Request 100 hour runtime
#$ -l h_vmem=10G

# combine all varaints (merging same ID genotype fields)
java -Xmx8g -jar ~/Software/GenomeAnalysisTK.jar -T CombineVariants \\
-R /data/BCI-EvoCa2/wchc/referenceFiles/ucsc.hg19.fasta \\
-genotypeMergeOptions UNSORTED \\
", paste("--variant", vcfFiles, collapse = " "), " \\
-o ", mergedVCF, "


#annotate merged vcf with VEP (getting sift and polyphen scores) 
module load ensembl-vep

vep -i ", mergedVCF," \\
-o ", annoVCF," \\
-format vcf -offline -cache -cache_version 93 dir_cache ~/.vep \\
-assembly GRCh37 -stats_text -stats_text -variant_class -sift b -polyphen b -tab -nearest gene -symbol

", sep="")
  
  
  #write seq script file
  lapply(totalStringsMerge, write, outNameMerge, append=FALSE)
  
  
  
  #make vcf pannel of normals (PON) script

  outNamePON <- paste(scriptsOut, "runMutect.makePON.", setNames[currSet], ".sh", sep="") 
  totalStringsMerge <- paste("#!/bin/sh
#$ -cwd
#$ -V
#$ -l h_rt=120:0:0
#$ -pe smp 1
#$ -l h_vmem=6G

######### run mutect2 on normal ##########
java -Xmx4g -jar /data/home/mpx155/bin/gatk-4.0.1.0/gatk-package-4.0.1.0-local.jar Mutect2 \\
-R ", subSample[currBio, "reference"], " \\
-I ", normalBam," \\
-tumor ", normalID," \\
--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \\
-O ", vcfOutPON," \\
", regionsString, sep="")

  
  #write seq script file
  lapply(totalStringsMerge, write, outNamePON, append=FALSE)
}







############  make mutect PON merge script #############


vcfNormals <- paste("-vcfs", vcfNormals)

outNamePONmerge <- paste(scriptsOut, "runMutect.makePON.createPannel.sh", sep="") 
totalStringsMerge <- paste("#!/bin/sh
#$ -cwd
#$ -V
#$ -l h_rt=120:0:0
#$ -pe smp 1
#$ -l h_vmem=6G

######### run mutect2 on normal ##########
java -Xmx4g -jar /data/home/mpx155/bin/gatk-4.0.1.0/gatk-package-4.0.1.0-local.jar CreateSomaticPanelOfNormals \\
", paste(vcfNormals, collapse=" "),"
-O /data/BCI-EvoCa-SG/3.0.mutectCalls/exomes/PON.vcf.gz                           

", sep="")

#write script
lapply(totalStringsMerge, write, outNamePONmerge, append=FALSE)

