\name{scanIBD2}
\alias{scanIBD2}
\title{
Scanning for IBD 2 sharing
}
\description{
Establishes pair-wise IBD 2 sharing between individuals in the cohort. 
}
\usage{
scanIBD2(iGene,chromEnds, destinationFile, lengthThresh, probands, trimMargin)
}
\arguments{
  \item{iGene}{A list of genotypes (one for each individual), casted to integers.}
  \item{chromEnds}{An array of integers, where ith entry gives indices of the last SNP on each chromosome.}
  \item{destinationFile}{ Path to file where the surrogacy matrix will be saved.}
  \item{lengthThresh}{IBD threshold, expressed as a number of consecutive SNPs without opposing homozygotes.}
  \item{probands}{An array with indices of individuals whose genotypes are to compared to others.}
  \item{trimMargin}{Number of SNPs cut at the borders.}
}
\details{
}
\value{
  \item{surrMat}{The surrogacy matrix showing all pair-wise sharing in a cohort. An i,j entry gives IBD 2 sharing between individuals i and j.}
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
\keyword{ IBD2 }
\keyword{ sharing }
