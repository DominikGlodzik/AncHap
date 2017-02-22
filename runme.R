library(anchap)

 # inputs
genotypesFile <- "data/Orca/orca749chr2/orca749chr2"; # path to bin, bed, fam files - binary PLINK format

# HapMap genetic map available to install from
# http://hapmap.ncbi.nlm.nih.gov/downloads/recombination/2011-01_phaseII_B37/genetic_map_HapMapII_GRCh37.tar.gz
genMap <- list('data/genetic_map_HapMapII_GRCh37/genetic_map_GRCh37_chr', '.txt'); # for build 37 only

# parameters
IBDTHRESH <- 5 # Stage I threshold, in centiMorgans, used 3 cM in Orcades study
TRIMMARGIN <- 100; #  how many markers trimmed from IBD segments on each side, default 100
OVERLAPT <- 10 # stitiching parameter, default 10 [SNPs]
MATCHT <- 0.02 # stitching parameter, default 0.02
RESULT.BASE.FOLDER<- 'results/runme/' 
# Stage III parameters
IBDTHRESHSCAN2 <- 2 # in centiMorgans - set to NULL if round 1 only
LENGTHTHRESH <- 200 # minimum number of phased SNPs

######################################################################################################
chr2Data <- loadData2(genotypesFile = genotypesFile, geneticSnpMapFile = genMap, build= 37) # load data

# prepare output paths
outputPaths <- list()                               
outputPaths$surrmatFile <- paste(RESULT.BASE.FOLDER, "Surrs-IBDTHRESH-",IBDTHRESH,".Rdata", sep=""); # path, where the surrogacy matrix from Stage I is stored                           
outputPaths$plotFolder <-  paste(RESULT.BASE.FOLDER,  sep="");  # path, to where the plots illustrating stitches and inconsistences go                                       
outputPaths$phasingOutputFolder <- paste(RESULT.BASE.FOLDER, sep=""); # path where the phased haplotypes will be saved                                       
outputPaths$trioResultPath <- RESULT.BASE.FOLDER # path where the results of the trio analysis will go
outputPaths$resultPath <- paste(RESULT.BASE.FOLDER, "result-IBDTHRESH-",IBDTHRESH,"-trim-",TRIMMARGIN,"-overlapT-", OVERLAPT, "-matchT-",MATCHT,"-lengthThreshold-",LENGTHTHRESH,".RData", sep=""); # where all results combined are saved
outputPaths$summaryPath <- NULL # file path for summary of results
outputPaths$paramEvalPath <- NULL # file path for another summary of results

# STAGES I,II of ANCHAP
result <- phaseGenDist(
                data = chr2Data, # data in the ANCHAP format, from the loadData2 function
                outputPaths = outputPaths,
                individualsPhased = 1:5, # for testing, only phase first 5 individuals, otherwise remove from parameter list (comment out this line)
                genThresh = IBDTHRESH, # Stage I threshold, in centiMorgans
                trimMargin=TRIMMARGIN, # comment out for default
                overlapT = OVERLAPT,  # comment out for default
                matchT=MATCHT, # comment out for default 
                lengthThresh=LENGTHTHRESH, # comment out for default 
               )

# selection of samples for re-sequencing
reseq.selection <- pickReseqIndsFast(result$surrMatDf, 2, plotPath='results/reseq.pdf', chr2Data, appFactor=10)
   



# Stage III of ANCHAP
if (is.null(IBDTHRESHSCAN2)) {
    save(result, 'results/runme/result.RData')  # save results of Stages I, II only  
} else {
    surrMatHapDf <- scanSurrsHapDf(result=result, lengthThresh=LENGTHTHRESH, chromEnds=chr2Data$chromEnds, geneticMap=chr2Data$geneticSnpMap, gDistT=IBDTHRESHSCAN2)
  surrHapDensity <- NA;
  if ((dim(surrMatHapDf)[1])>0) { # if found any sharing in stage III
    surrHapDensity <- computeSurrDensity(surrMatHapDf)}
    settingsString <- paste('scan2thresh', IBDTHRESHSCAN2, sep='-')
    fileName <- paste(RESULT.BASE.FOLDER,"surrMatHapDf", '-', settingsString,".RData",sep="");
    save(result, surrMatHapDf, surrHapDensity, file= fileName); # save results of Stages I, II, III 
}  


