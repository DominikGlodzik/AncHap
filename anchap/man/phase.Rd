\name{phase}
\alias{phase}
\title{
Phasing function, with an IBD threshold specified as a count of SNPs without opposing homozygotes
}
\description{
First establishes pair-wise IBD sharing between individuals in the cohort. 
Then for each individual performs stitching and phasing. 
Stitching is alignment of the shared sequences for an individual into two categories, 
so that in each there are sequences sharing with one haplotype. 
Phasing uses the shared sequences, and their genotype, to recover haplotypes.
Finally, the procedure is evaluated using the parent-offspring trios present in the cohort.
}
\usage{
phase(data, outputPaths=NULL, individualsPhased = ( 1:(dim(data$genotypeMatrix)[1])), ibdThresh=1000, onlySurrs=FALSE, excludeChildren = TRUE, excludeParents = FALSE,  trimMargin=20, overlapT=50)
}

\arguments{
  \item{data}{Data for the cohort, in standard Anchap data format.}
  \item{outputPaths}{File paths for saving results and plots. See details.}
  \item{individualsPhased}{Indices of the individual to be phased.}
  \item{ibdThresh}{IBD threshold, expressed as a count of SNPs.}
  \item{onlySurrs}{Whether the method should find surrogate parents only, without phasing.}
  \item{excludeChildren}{Whether children should be excluded as surrogate parents. Avoids spurious detections.}
  \item{excludeParents}{Whether real parents should be excluded. Useful for trio testing.}
  \item{trimMargin}{Count of markers at the borders of shared regions that are trimmed.}
  \item{overlapT}{Overlap threshold between two shared sequences, expressed as a count of consecutive markers. Used for alignment of shared sequences.}
}

\details{
    The argument outputPaths is a list, and should have the following elements: surrmatFile, plotFolder, phasingOutputFolder, trioResultPath. If outputPaths is NULL, no results will be saved.}


\value{
    A list with the following fields:
  \item{surrMat}{The surrogacy matrix showing all pair-wise sharing in a cohort. An i,j entry gives sharing between individuals i and j.}
  \item{stitchingResults }{The results summarsing the stitching procedure, for all individuals.}
  \item{phasingResults}{The results summarsing phasing of all individuals. Contrains the actual haplotypes, as well as amounts of inconsistencies.}
  \item{trioResults}{Summary of the evaluation if phasing using parent-offspring trios.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
    Dominik Glodzik
}
\note{
}
\seealso{
}
\examples{
}
\keyword{ phasing }