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

loadData <- function(mapFile,
                     genotypesFile,
                     pedigreeFilepath,
                     geneticSnpMapFile = NULL,
                     iGeneFile = NULL,
                     iGeneCreate=TRUE,
                     noChroms=22,
                     noInds=NULL,
                     ids = NULL
                     ) {

  read.pedfile.map <- function (file) {
    result <- read.table(file, col.names = c("chromosome", "snp.names", 
                                 "genetic.distance", "position"), as.is = TRUE, fill = TRUE, 
                         flush = TRUE)
    rownames(result) <- result$snp.names
    result[, c("snp.names", "position", "chromosome")]
  }
  
  map <- read.pedfile.map(mapFile);
 
  chromEnds <- unlist(extractEndpointsSnp(map));

  cumEnds <- cumsum(as.numeric(extractEndpoints(map)))
  genotypeMatrix <- read.plink(genotypesFile);
  
  g1 <- genotypeMatrix[246,]
   g2 <- genotypeMatrix[192,]
  # subset the genotype matrix
  if (!is.null(ids)) {
    indsIn <- which(rownames(genotypeMatrix) %in% ids)
    genotypeMatrix <- genotypeMatrix[indsIn,]
    print(dim(genotypeMatrix))
  }
  # count the number of individuals
 
 
  
  if (is.null(noInds)) {noInds <- dim(genotypeMatrix)[1]};
  probands <- 1:noInds

  orcaIds <- rownames(genotypeMatrix)
  idNames <- rownames(genotypeMatrix)

  pedigree <- NULL;
  if (!is.null(pedigreeFilepath))
    pedigree <- read.xls(pedigreeFilepath)
  end;

    
  
  if (iGeneCreate) {
    iGene <- list()
    for (i in probands) { print(i)
                          iGene[[i]] <- as.integer(genotypeMatrix[i, ]) }
    if (!is.null(iGeneFile)) {save(iGene, file=iGeneFile )} }
  else {load(iGeneFile)}


  geneticSnpMap <- NULL
  if (!is.null(geneticSnpMapFile)) {

    print("Localising the SNPs on the genetic map")
    gm <- readGenMap(filePath=geneticSnpMapFile);

    # a vector for each SNP
    geneticSnpMap <- vector()
    # SNP index
    gI <- 1;

    # loop over the chromosomes
    for (chrom in 1:22) {

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
         if ((gI %% 1000)==0) {print(gI)}
       }
     }
     # end loop over SNPs
   }
   # end loop over chromosomes
  }
  # end if genetic map file specified

  
  print("building the genotype data structure ... ")
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

  
