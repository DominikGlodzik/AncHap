# loads
# map - map of the SNPs' locations on the genome
# chromEnds - the index of the last SNP on each chromosome, indexed by chromosomes
# cumEnds - the positions of chromosomes on the genome, expressed in bp, used for plotting

# genotypeMatrix - a matrix of alleles, indexed y individuals, and SNPs
# iGene - genotypes of individuals stored as integers, for easier passing to C; a liste indexed by individuals

# probands - indices of individuals in the cohort
# orcaIds, idNames - IDs of the individuals in the cohort
# pedigree - pedigree for the cohort - rows of individuals whose parents are specified
# can be null
# chromList - a 1-D representation of genotypes as integers; a list indexed by chromosomes
# noChroms=22,
# noInds=NULL

loadDataGenABEL <- function(genotypeFileS,
                            phenotypesFile=NULL,
                     geneticSnpMapFile = NULL,
                     noChroms=22,
                     noInds=NULL,
                     ids = NULL,
                     iGeneCreate=TRUE,
                     iGeneFile=NULL,
                     geneticMap=NULL,
                     snpMapBiomart=NULL
                     ) {


  readGAFile <- function(genotypesFile) {
    if ((substr(genotypesFile, nchar(genotypesFile)-3, nchar(genotypesFile))==".raw") && !(is.null(phenotypesFile))) {
                                        # read in the data
      alldata <- load.gwaa.data(phe= phenotypesFile, gen=genotypesFile, force=T, makemap=F, sort=F)
    } else if (substr(genotypesFile, nchar(genotypesFile)-5, nchar(genotypesFile))==".RData") {
                                        # VDR data
      vNames <- load(genotypesFile)
      eval(parse(text = paste('alldata <-', vNames[1]))) 
    }
    alldata
  }

  mergeGwaa <- function(x,y) {
    snpdata <- merge(x@gtdata,y@gtdata)$data
    #phdata <- merge(x@phdata,y@phdata,all=T)
    #rownames(phdata) <- phdata$id
    #phdata <- phdata[snpdata@idnames,]
    out <- new("gwaa.data",phdata=data.frame(),gtdata=snpdata)
    out
  }

  dataList <- list()    
  noFiles <- length(genotypeFileS)   
  for (f in 1:noFiles) {
    dataList[[f]] <- readGAFile(genotypeFileS[f])
  }

  df <- dataList[[1]]
  if (noFiles>1) {
    for (f in 2:noFiles) {
      # do the merging
      df <- mergeGwaa(df, dataList[[f]])
    }
  }

  # update SNP positions to ones from an external file
  if (!is.null(snpMapBiomart)) {
    mapTable <- read.csv(snpMapBiomart, header=TRUE)
    genabelSnpNames <- colnames(df@gtdata)
    genabelNewMap<- merge(data.frame(Name=genabelSnpNames),mapTable, all.x=TRUE, sort=FALSE, by.x='Name', by.y='Name' )

    df@gtdata@map<-genabelNewMap$Position
  }
  
  # reorder the SNPs according to their map position 
  snpMapOrder <- order(df@gtdata@map)
  orderedData <- df[,snpMapOrder]
  df@gtdata@map <- df@gtdata@map[snpMapOrder]
  alldata <- orderedData
  
  # here would stop and save the data
  # save.gwaa.data(alldata, genofile = "../Data/data/VDR/VDR_region_JOINT.raw", human = FALSE)
  # export them to PLINK
  # export.plink(alldata, filebasename = "../Data/data/VDR/VDR_region_JOINT", phenotypes=NULL)
  
    #newOrder <- order(dfVDR@phdata$Disease)  
    #orderedData <- dfVDR[newOrder,]
  
    
  

  
  # converting the genetic map
  map <- data.frame(position=as.numeric(alldata@gtdata@map), chromosome=as.integer(as.character(alldata@gtdata@chromosome)))
  # converting the genotype matrix
  genotypeMatrix <- as.numeric(alldata@gtdata) + 1;
  genotypeMatrix[is.na(genotypeMatrix)] <- 0

  # does the order of SNPs in the genotype matrix and 
  matrixSnps <- colnames(genotypeMatrix)
  mapSnps <- colnames(alldata@gtdata@map)
  mismatches <- sum(matrixSnps!=mapSnps)
  
  
  # sort the genetic map
  newOrder <- order(map$chromosome, map$position)
  # was the genetic map ordered in the first place?
  wasOrdered <- sum(newOrder!=(1:length(newOrder)))==0
  map <- map[newOrder,]
  genotypeMatrix <- genotypeMatrix[,newOrder]

  chromEnds <- unlist(extractEndpointsSnp(map));
  cumEnds <- cumsum(as.numeric(extractEndpoints(map)))

  
  #####
  pedigree <- NULL;
   

  # subset the genotype matrix
  if (!is.null(ids)) {
    indsIn <- which(rownames(genotypeMatrix) %in% ids)
    genotypeMatrix <- genotypeMatrix[indsIn,]
    # print(dim(genotypeMatrix))
  }
  # count the number of individuals
 
 
  
  if (is.null(noInds)) {noInds <- dim(genotypeMatrix)[1]};
  probands <- 1:noInds

  orcaIds <- rownames(genotypeMatrix)
  idNames <- rownames(genotypeMatrix)

  cat('Converting data ... \n')
  if (iGeneCreate) {
    iGene <- list()
    for (i in probands) { if ((i %% 10)==0) {cat(paste('Individual',i, 'out of', length(probands)  ,'\n'))}
                          iGene[[i]] <- as.integer(genotypeMatrix[i, ]) }
    if (!is.null(iGeneFile)) {save(iGene, file=iGeneFile )} }
  else {load(iGeneFile)}


  gm<-NULL
  geneticSnpMap <- NULL
  if (!is.null(geneticSnpMapFile)) {
    gm <- readGenMap(filePath=geneticSnpMapFile);
  }
  if (!is.null(geneticMap)) {
    gm <- geneticMap
  }

  if(!is.null(gm)) {
  
    
    # a vector for each SNP
    geneticSnpMap <- vector()
    # SNP index
    gI <- 1;

    # loop over the chromosomes
    cat("Localising the SNPs on the genetic map \n")
    for (chrom in 1:22) {

          print(paste("chromosome", chrom))
      
     chromMap <- gm$mapCm[[chrom]]
     mapPosBp <- gm$mapPosBpRaw[[chrom]]
     # SNP map on the chromosome
     snpMap <- map[map$chromosome==chrom,]

     # loop over the SNPs on the chromosome
     
     if (dim(snpMap)[1]>0) {
       for (s in 1:((dim(snpMap)[1]))) {
                                        # which entry on the genetic map is closest to the SNP
         mapEntry1 <- which.min(abs(mapPosBp - snpMap$position[s]))
         
         geneticSnpMap[gI] <-chromMap[mapEntry1];
         gI <- gI + 1;
       }
     }
     # end loop over SNPs
   }
   # end loop over chromosomes
  }
  # end if genetic map file specified

  
  print("building the chromosome data structure ... ")
  chromList <<- vector("list",noChroms)
  for (chromosomeNo in 1:noChroms) {
    print(paste("chromosome", chromosomeNo))

    if (chromosomeNo>1) {firstS <- chromEnds[chromosomeNo-1]+1} else {firstS <- 1;}
    lastS <- chromEnds[chromosomeNo]

    if (lastS>firstS) {
    snpsAnalysed <- firstS:lastS;
    noSnpsChrom <- as.integer(length(snpsAnalysed))
    noGenes <- noInds * noSnpsChrom

    genesVec <<- integer(noGenes)
    for (i in 1:noInds) {
      startIndx <- (i-1)*noSnpsChrom + 1;
      endIndx <- i*noSnpsChrom;
      genesVec[startIndx:endIndx] <- iGene[[i]][snpsAnalysed]   
    }
    chromList[[chromosomeNo]] <- genesVec;
    rm(genesVec)
    }
  }

  data <- list();
  data$iGene <- iGene;
  data$map <- map;
  data$chromEnds <- chromEnds
  data$cumEnds <-cumEnds;
  data$genotypeMatrix <- genotypeMatrix
  data$probands <- probands
  data$idNames <- idNames;
  data$pedigree <- pedigree;
  data$chromList <- chromList;
  data$orcaIds <- orcaIds
  data$geneticSnpMap <-geneticSnpMap
  data
}
