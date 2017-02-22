\name{phaseGenDist}
\alias{phaseGenDist}
\title{
Stages I and II of ANCHAP.
}
\description{
First establishes pair-wise IBD sharing between individuals in the cohort (Stage I).
Then for each individual performs stitching and phasing (Stage II). 
Stitching is alignment of the shared sequences for an individual into two categories, 
so that in each there are sequences sharing with one haplotype. 
Phasing uses the shared sequences, and their genotype, to recover haplotypes.
Finally, the procedure is evaluated using the parent-offspring trios present in the cohort, where the pedigree data is available.
}
\usage{
phaseGenDist(data, outputPaths, individualsPhased, genThresh=3, onlySurrs=FALSE,
                  excludeChildren = TRUE, excludeParents = FALSE,  trimMargin=100, overlapT=10, matchT=0.02,  lengthThresh= 200, trioAnalysis=FALSE, surrMatFile=NULL)
}
\arguments{
  \item{data}{Data for the cohort, in standard ANCHAP data format. As prepared by the loadData2 function.}
  \item{outputPaths}{File paths for saving results and plots. A list of character vectors with following fields: surrmatFile, plotFolder, phasingOutputFolder, trioResultPath, resultPath, summaryPath, paramEvalPath. Can be set to NULL if a particular type of result is not meant to be saved.}
  \item{individualsPhased}{Indices of the individual to be phased. For example, for testing, one may use only first 20 individuals, ie 1:20. If not specified, all samples would be used.}
  \item{genThresh}{IBD threshold, expressed in cM. A numeric value. Whether the method should find surrogate parents only, without phasing. TRUE or FALSE.}
  \item{onlySurrs}{Whether the method should find surrogate parents only, without phasing. Only works when pedigree is available. TRUE/FALSE.}
  \item{excludeChildren}{Whether children should be excluded as surrogate parents. Avoids spurious detections.}
  \item{excludeParents}{Whether real parents should be excluded. Useful for trio testing. Only works when pedigree is available. TRUE/FALSE.}
  \item{trimMargin}{Count of markers at the borders of shared regions that are trimmed. Number of markers trimmed at each border of each IBD segment. Default 100 SNPs.}
  \item{overlapT}{Overlap threshold between two shared sequences, expressed as a count of consecutive markers. Used for alignment of shared sequences. Default 10 SNPs.}
  \item{matchT}{Minimum compliance between stitched IBD segments. Minimum complience between stitched IBD segments. Default 0.02.}
  \item{trioAnalysis}{Whether to evaluate the resulting haplotypes by trio analysis. TRUE/FALSE}
  \item{surrMatFile}{Location of the file for surrogacy matrix, if this is meant to be saved separately. Character array.}
}
\details{
}
\value{
  \item{result}{A list with all results of ANCHAP.}
}
  
\author{
	Dominik Glodzik
}
\seealso{
}
\examples{
}
\keyword{ phasing }
\keyword{ genetic }
\keyword{ distance }