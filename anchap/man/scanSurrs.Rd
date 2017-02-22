\name{scanSurrs}
\alias{scanSurrs}
\title{
Scanning for surrogate parents
}
\description{
Establishes pair-wise IBD sharing between individuals in the cohort. 

}
\usage{
scanSurrs(iGene,chromEnds, destinationFile=NULL, lengthThresh, probands, trimMargin, geneticMap=NULL, gDistT=NULL)
}
\arguments{
  \item{iGene}{A list of genotypes (one for each individual), casted to integers.}
  \item{chromEnds}{An array of integers, where ith entry gives indices of the last SNP on each chromosome.}
  \item{destinationFile}{ Path to file where the surrogacy matrix will be saved.}
  \item{lengthThresh}{IBD threshold, expressed as a number of consecutive SNPs without opposing homozygotes.}
  \item{probands}{An array with indices of individuals whose genotypes are to compared to others.}
  \item{trimMargin}{Number of SNPs cut at the borders.}
  \item{geneticMap}{The genetic map, as in Anchap's standard data structure.}
  \item{gDistT}{ IBD threshold expressed in cM. Will be applied jointly with lengthThresh.}
}
\details{
}
\value{
  \item{surrMat}{The surrogacy matrix showing all pair-wise sharing in a cohort. An i,j entry gives sharing between individuals i and j.}
}
\references{
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ surrogate}
\keyword{ parents }
\keyword{ surrogacy }
\keyword{ matrix }