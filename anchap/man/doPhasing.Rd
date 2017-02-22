\name{doPhasing}
\alias{doPhasing}
\title{
Phases a genotype of an individual.
}
\description{
Phases a genotype of an individual, to give two haplotypes, using the IBD sharing between the individual, and others. 
}
\usage{
doPhasing(probandId, iGene, chromList, selectedSequencesH1, selectedSequencesH2, alsoPlot, plotFolder, outputFolder, cumEnds, map, chromEnds)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probandId}{
The index of the individual to be phased.
}
  \item{iGene}{
  Genotypes of individuals in the cohort. A list indexed by individuals, whose entry contains an individual's genotype, represented by an array of integers.
}
  \item{chromList}{
   A 1-D representation of genotypes as integers; a list indexed by chromosomes.
}
  \item{selectedSequencesH1}{
The trimmed IBD sequences of the individual classified to group 1.
}
  \item{selectedSequencesH2}{
The trimmed IBD sequences of the individual classified to group 2.
}
  \item{alsoPlot}{
TRUE/FALSE - whether plots illustrating stitching and phasing should be generated and saved to disk.
}
  \item{plotFolder}{
The desination folder for the plots.
}
  \item{outputFolder}{
Folder where the recovered haplotypes will be saved as RData. 
}
  \item{cumEnds}{
Positions of the last SNP of each chromosome on the genome, needed for plotting.
}
  \item{map}{
  Map of the SNPs' locations on the genome.
}
  \item{chromEnds}{
  The index of the last SNP on each chromosome, indexed by chromosomes.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 \item{h1}{Haplotype 1 of the phased individual - an array of the same length as the genotype. }
 \item{h2}{Haplotype 2 of the phased individual - an array of the same length as the genotype.}
 \item{m1}{Most common allele at homozygous loci among the surrogate parents in group 1 - an array of the same length as the genotype. Used for detecting between haplotype inconsistencies.}
 \item{m2}{Most common allele at homozygous loci among the surrogate parents in group 2.}
 \item{noSurrs1}{Number of surrogates in group 1 at a locus.}
 \item{noSurrs2}{Number of surrogates in group 2 at a locus.}
 \item{i1}{Number of within-haplotype inconsistencies for group 1.}
 \item{i2}{Number of within-haplotype inconsistencies for group 2.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (probandId, iGene, chromList, selectedSequencesH1, selectedSequencesH2, 
    alsoPlot, plotFolder, outputFolder, cumEnds, map, chromEnds) 
{
    noSnps <- length(iGene[[1]])
    shiftedCumEnds <- c(0, cumEnds)
    iMap <- map[, "position"] + shiftedCumEnds[map[, "chromosome"]]
    doChromPhasing <- function() {
        r <- .C("doFastPhasing", noSnpsChrom, firstS, lastS, 
            firstS, lastS, probandId, noInds, genesVec, noH1Seqs, 
            h1Starts, h1Ends, h1Ids, noH2Seqs, h2Starts, h2Ends, 
            h2Ids, noH1atLocus, noH2atLocus, haplotype1, haplotype2, 
            majority1, majority2, genotype, incons1, incons2)
        noH1atLocus <<- r[[17]]
        noH2atLocus <<- r[[18]]
        majority1 <<- r[[21]]
        majority2 <<- r[[22]]
        haplotype1 <<- r[[19]]
        haplotype2 <<- r[[20]]
        genotype <<- r[[23]]
        incons1 <<- r[[24]]
        incons2 <<- r[[25]]
        rm(r)
        gc()
    }
    getChromStart <- function(chromosomeNo) {
        if (chromosomeNo == 1) {
            value <- as.integer(1)
        }
        else {
            value <- as.integer(chromEnds[chromosomeNo - 1] + 
                1)
        }
        value
    }
    getChromEnd <- function(chromosomeNo) {
        as.integer(chromEnds[chromosomeNo])
    }
    getChromNo <- function(firstS, lastS) {
        chromMap <- map$chromosome[firstS:lastS]
        chromNo <- max(chromMap)
        if (chromNo != min(chromMap)) 
            print("Surrogate seq on different chromosomes!")
        chromNo
    }
    probandId <<- as.integer(probandId)
    noInds <<- length(iGene)
    noH1Seqs <<- dim(selectedSequencesH1)[1]
    h1Starts <<- as.numeric(selectedSequencesH1$start)
    h1Ends <<- as.numeric(selectedSequencesH1$end)
    h1Ids <<- selectedSequencesH1$id
    noH2Seqs <<- dim(selectedSequencesH2)[1]
    h2Starts <<- as.numeric(selectedSequencesH2$start)
    h2Ends <<- as.numeric(selectedSequencesH2$end)
    h2Ids <<- selectedSequencesH2$id
    h1 <<- integer(length(iGene[[1]]))
    h2 <<- integer(length(iGene[[1]]))
    m1 <<- integer(length(iGene[[1]]))
    m2 <<- integer(length(iGene[[1]]))
    noSurrs1 <<- integer(length(iGene[[1]]))
    noSurrs2 <<- integer(length(iGene[[1]]))
    i1 <<- integer(length(iGene[[1]]))
    i2 <<- integer(length(iGene[[1]]))
    t1 <- proc.time()
    genesVec <<- vector()
    for (chromosomeNo in 1:22) {
        genesVec <- chromList[[chromosomeNo]]
        firstS <<- getChromStart(chromosomeNo)
        lastS <<- getChromEnd(chromosomeNo)
        if (lastS > firstS) {
            snpsAnalysed <<- firstS:lastS
            noSnpsChrom <<- as.integer(length(snpsAnalysed))
            noGenes <<- noInds * noSnps
            noH1atLocus <<- integer(lastS - firstS + 1)
            noH2atLocus <<- integer(lastS - firstS + 1)
            haplotype1 <<- integer(lastS - firstS + 1)
            haplotype2 <<- integer(lastS - firstS + 1)
            majority1 <<- integer(lastS - firstS + 1)
            majority2 <<- integer(lastS - firstS + 1)
            genotype <<- integer(lastS - firstS + 1)
            incons1 <<- integer(lastS - firstS + 1)
            incons2 <<- integer(lastS - firstS + 1)
            n1 <<- integer(lastS - firstS + 1)
            n2 <<- integer(lastS - firstS + 1)
            doChromPhasing()
            genesVec = NULL
            h1[firstS:lastS] <<- haplotype1
            h2[firstS:lastS] <<- haplotype2
            m1[firstS:lastS] <<- majority1
            m2[firstS:lastS] <<- majority2
            i1[firstS:lastS] <<- incons1
            i2[firstS:lastS] <<- incons2
            noSurrs1[firstS:lastS] <<- noH1atLocus
            noSurrs2[firstS:lastS] <<- noH2atLocus
            noH1atLocus <<- NULL
            noH2atLocus <<- NULL
            haplotype1 <<- NULL
            haplotype2 <<- NULL
            majority1 <<- NULL
            majority2 <<- NULL
            genotype <<- NULL
            incons1 <<- NULL
            incons2 <<- NULL
            n1 <<- NULL
            n2 <<- NULL
            gc()
            mismatches <- sum(genotype[1:length(genotype)] != 
                iGene[[probandId]][firstS:lastS])
            if (mismatches > 0) {
                print("indexing error")
            }
        }
    }
    print(proc.time() - t1)
    hapMismatches1 <- sum((h1[1:length(h1)] != 3) & (iGene[[probandId]][1:length(h1)] == 
        3))
    hapMismatches2 <- sum((h2[1:length(h1)] != 1) & (iGene[[probandId]][1:length(h1)] == 
        1))
    if ((hapMismatches1 + hapMismatches1) > 0) {
        print("haplotyping error error")
    }
    betHapMismatches <<- (iGene[[probandId]][1:length(h1)] == 
        2) & (m1[1:length(h1)] != 0) & (m1[1:length(h1)] == m2[1:length(h1)])
    selectedSequences <- rbind(selectedSequencesH1, selectedSequencesH2)
    coverageAll <- calcCoverage(allSurrogates, noSnps)
    coverageSel <- calcCoverage(selectedSequences, noSnps)
    noSurrAtLoc <- mean(noSurrs1 + noSurrs2)
    allSeqsCoverage <- sum(allSurrogates$end - allSurrogates$start)/noSnps
    selectedH1Coverage <- sum(selectedSequencesH1$end - selectedSequencesH1$start)/noSnps
    selectedH2Coverage <- sum(selectedSequencesH2$end - selectedSequencesH2$start)/noSnps
    phasingPercentage <- sum(h1 != 0)/length(h1)
    noIncons1 <- sum(i1 > 0)
    noIncons2 <- sum(i2 > 0)
    betIncons <- sum(betHapMismatches > 0)
    aResult <- data.frame(id = probandId, phased = phasingPercentage, 
        noSurrsAtLocus = noSurrAtLoc, noStitches = as.numeric(coverageSel[2]), 
        allSeqsCoveragePerc = as.numeric(coverageAll[1]), selectedSeqsCoveragePerc = as.numeric(coverageSel$coverage), 
        allSeqsCoverage = allSeqsCoverage, selectedH1Coverage = selectedH1Coverage, 
        selectedH2Coverage = selectedH2Coverage, inconsH1 = noIncons1, 
        inconsH2 = noIncons2, betweenHincons = betIncons)
    save(h1, h2, m1, m2, i1, i2, noSurrs1, noSurrs2, aResult, 
        selectedSequencesH1, selectedSequencesH2, rejectedSequences, 
        file = paste(outputFolder, "results", probandId, ".Rdata", 
            sep = ""))
    if (alsoPlot) {
        plotPhasing(probandId, plotFolder, map, cumEnds, iMap, 
            noInds)
    }
    result <- list()
    result[["h1"]] <- h1
    result[["h2"]] <- h2
    result[["m1"]] <- m1
    result[["m2"]] <- m2
    result[["noSurrs1"]] <- noSurrs1
    result[["noSurrs2"]] <- noSurrs2
    result[["i1"]] <- i1
    result[["i2"]] <- i2
    result
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
