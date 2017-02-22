\name{doStitching}
\alias{doStitching}
\title{Clustering all the IBD sequences of an individual into two groups.}
\description{
For and individual and its IBD sequences, clusters the sequences into two groups. 
The sequences in the first group share 1 with haplotype 1 of the proband,
while the sequences from the second group share with haplotype 2. 
}
\usage{
doStitching(individualPhased, surrMat, iGene, chromList, surrogateMargin, excludedParents, overlapT, matchT, map, chromEnds)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{individualPhased}{
The individuals whose IBD sequences are to be clustered.
}
  \item{surrMat}{
The surrogacy matrix showing all pair-wise sharing in a cohort. An i,j entry gives sharing between individuals i and j.
}
  \item{iGene}{
  Genotypes of individuals in the cohort. A list indexed by individuals, whose entry contains an individual's genotype, represented by an array of integers.
}
  \item{chromList}{
   A 1-D representation of genotypes as integers; a list indexed by chromosomes.
}
  \item{surrogateMargin}{
  Number of SNPs trimmed from the IBD sequences at their extremes.
}
  \item{excludedParents}{
  Identifiers of individuals excluded from the anlysis for an individual - typically indicies of their choldren.
}
  \item{overlapT}{
  Number of consectutive SNPs of a smallest overlap, when two sequences can be stitched.
}
  \item{matchT}{
  Maximum allowed ratio of opposing homozygotes between and IBD sequence, and a haplotype, over the number of heterozygous loci in the overlap region.
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
\item{summary}{Statistics of how the IBD sequences matched each other.}
\item{selectedSequencesH1}{The trimmed sequences classified to group 1.}
\item{selectedSequencesH2}{The trimmed sequences classified to group 2.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Dominik Glodzik
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
function (individualPhased, surrMat, iGene, chromList, surrogateMargin, 
    excludedParents, overlapT, matchT, map, chromEnds) 
{
    getSeqsWithin <- function(seq, selectedSequences) {
        sStart <- seq$start
        sEnd <- seq$end
        chosenSeqs <- selectedSequences$start <= sEnd & selectedSequences$end >= 
            sStart
        chosenSeqs <- selectedSequences[chosenSeqs, ]
    }
    countOppHomoR <- function(s1, s2) {
        if (length(s1) != length(s2)) {
            stop()
        }
        zeros1 <- (s1 == 0)
        zeros2 <- (s2 == 0)
        oppHomos <- (abs(s1 - s2) == 2) & (!zeros1) & (!zeros2)
        sum(oppHomos)
    }
    getOverlapLength <- function(individualPhased, firstH, secondH) {
        overlapStart <- max(firstH$start, secondH$start)
        overlapEnd <- min(firstH$end, secondH$end) - 1
        noHeteroInOverlap <- sum(iGene[[individualPhased]][overlapStart:overlapEnd] == 
            2)
        noHeteroInOverlap
    }
    assignToHaplotypes <- function(allSurrogate) {
        countEmpty <- 0
        countMatchBoth <- 0
        countMatchH1 <- 0
        countMatchH2 <- 0
        countMatchNone <- 0
        countNoOverlap <- 0
        hapMatching <- data.frame()
        allSurrCount <- dim(allSurrogates)[1]
        for (s in 1:allSurrCount) {
            if ((s\%\%50) == 0) {
                print(paste("surrogate sequence ", s, " of ", 
                  allSurrCount))
            }
            cSurr <- allSurrogates[s, ]
            seqsH1 <- getSeqsWithin(cSurr, selectedSequencesH1)
            seqsH2 <- getSeqsWithin(cSurr, selectedSequencesH2)
            allSeqs <- rbind(seqsH1, seqsH2)
            if (dim(allSeqs)[1] == 0) {
                selectedSequencesH1 <<- rbind(selectedSequencesH1, 
                  cSurr)
                countEmpty <- countEmpty + 1
            }
            else {
                overlapLength <- -1
                overlapSufficient <- FALSE
                for (q in 1:dim(allSeqs)[1]) {
                  aSurr <- allSeqs[q, ]
                  r <- getOverlapLength(individualPhased, aSurr, 
                    cSurr)
                  if (r >= overlapT) {
                    overlapSufficient <- TRUE
                    overlapLength <- r
                  }
                }
                if (overlapSufficient) {
                  haplotypes <- doPartialPhasing(individualPhased, 
                    cSurr$start, cSurr$end, selectedSequencesH1, 
                    selectedSequencesH2, chromEnds, map, iGene)
                  genotype <- iGene[[individualPhased]][cSurr$start:(cSurr$end - 
                    1)]
                  keyLoci <- sum((genotype == 2) & (haplotypes[[1]] != 
                    0))
                  oppHomo1 <- countOppHomoR(iGene[[cSurr$id]][cSurr$start:(cSurr$end - 
                    1)], haplotypes[[1]])
                  matchH1 <- ((oppHomo1/keyLoci) < matchT)
                  oppHomo2 <- countOppHomoR(iGene[[cSurr$id]][cSurr$start:(cSurr$end - 
                    1)], haplotypes[[2]])
                  matchH2 <- ((oppHomo2/keyLoci) < matchT)
                  if (keyLoci > overlapT) {
                    aMatch <- data.frame(matchH1Rate = (oppHomo1/keyLoci), 
                      matchH2Rate = (oppHomo2/keyLoci), matchH1 = matchH1, 
                      matchH2 = matchH2, keyLoci = keyLoci, overlapLength = overlapLength, 
                      start = cSurr$start, end = cSurr$end)
                    hapMatching <- rbind(hapMatching, aMatch)
                    if (matchH1 & matchH2) {
                      countMatchBoth <- countMatchBoth + 1
                      matchingBothSequences <<- rbind(matchingBothSequences, 
                        cSurr)
                    }
                    else {
                      if (matchH1) {
                        selectedSequencesH1 <<- rbind(selectedSequencesH1, 
                          cSurr)
                        countMatchH1 <- countMatchH1 + 1
                      }
                      else if (matchH2 & (keyLoci > overlapT)) {
                        selectedSequencesH2 <<- rbind(selectedSequencesH2, 
                          cSurr)
                        countMatchH2 <- countMatchH2 + 1
                      }
                      else {
                        rejectedSequences <<- rbind(rejectedSequences, 
                          cSurr)
                        countMatchNone <- countMatchNone + 1
                        cSeqH1 <- extractIBDnoChrom(iGene[[cSurr$id]][cSurr$start:(cSurr$end - 
                          1)], haplotypes[[1]], 100)
                        if (dim(cSeqH1)[1] > 0) {
                          cSeqH1$start <- cSeqH1$start + cSurr$start - 
                            1
                          cSeqH1$end <- cSeqH1$end + cSurr$start - 
                            1
                          cSeqH1 <- cbind(cSeqH1, id = cSurr$id)
                          partialSequencesH1 <<- rbind(partialSequencesH1, 
                            cSeqH1)
                        }
                        cSeqH2 <- extractIBDnoChrom(iGene[[cSurr$id]][cSurr$start:(cSurr$end - 
                          1)], haplotypes[[2]], 100)
                        if (dim(cSeqH2)[1] > 0) {
                          cSeqH2$start <- cSeqH2$start + cSurr$start - 
                            1
                          cSeqH2$end <- cSeqH2$end + cSurr$start - 
                            1
                          cSeqH2 <- cbind(cSeqH2, id = cSurr$id)
                          partialSequencesH2 <<- rbind(partialSequencesH2, 
                            cSeqH2)
                        }
                      }
                    }
                  }
                  else {
                    countNoOverlap <- countNoOverlap + 1
                    noOverlapSequences <<- rbind(noOverlapSequences, 
                      cSurr)
                  }
                }
                else {
                  countNoOverlap <- countNoOverlap + 1
                  noOverlapSequences <<- rbind(noOverlapSequences, 
                    cSurr)
                }
            }
        }
        result <- list()
        result$matchingSummary <- data.frame(allSurrCount = allSurrCount, 
            countEmpty = (countEmpty/allSurrCount), countMatchBoth = (countMatchBoth/allSurrCount), 
            countMatchH1 = (countMatchH1/allSurrCount), countMatchH2 = (countMatchH2/allSurrCount), 
            countMatchNone = (countMatchNone/allSurrCount), countNoOverlap = (countNoOverlap/allSurrCount))
        result$sequenceMatching <- hapMatching
        if ((dim(hapMatching)[1]) > 0) {
            rejected <- (!hapMatching$matchH1) & (!hapMatching$matchH2)
            result$rejectedSeqs <- hapMatching[rejected, ]
        }
        result
    }
    selectedSequencesH1 <<- data.frame()
    selectedSequencesH2 <<- data.frame()
    rejectedSequences <<- data.frame()
    matchingBothSequences <<- data.frame()
    noOverlapSequences <<- data.frame()
    partialSequencesH1 <<- data.frame()
    partialSequencesH2 <<- data.frame()
    source("utils/getAllSurrogatesOf.R")
    allSurrogates <<- getAllSurrogatesOf(surrMat, individualPhased, 
        surrogateMargin)
    notExcludedParents <- !(allSurrogates$id \%in\% excludedParents)
    allSurrogates <<- allSurrogates[notExcludedParents, ]
    sortedSurrogates <- data.frame()
    if ((dim(allSurrogates)[1]) > 0) {
        sortedSurrogates <- allSurrogates[order(-allSurrogates$end + 
            allSurrogates$start), ]
    }
    result <- list()
    result[["summary"]] <- assignToHaplotypes(sortedSurrogates)
    result[["selectedSequencesH1"]] <- selectedSequencesH1
    result[["selectedSequencesH2"]] <- selectedSequencesH2
    result
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
