\name{loadData2}
\alias{loadData2}
\title{
Loading data required by ANCHAP.
}
\description{
Loads all the data for anchap: genetic data in plink format, optionally with genetic map.
}
\usage{
loadData2(genotypesFile=pathToPlinkFiles, geneticSnpMapFile=pathToGeneticMap, build=37)
}
\arguments{
  \item{genotypesFile}{Path to Plink binary files (with extentions bim, bed ans fam) with the genotypes for the cohort.}
  \item{geneticSnpMapFile}{Path to a HapMap genetic map. Download it, unpack and point to the folder.}
}
\details{
}
\value{
A list with the following fields.
item{iGene}{Genotypes of individuals in the cohort. A list indexed by individuals, whose entry contains an individual's genotype, represented by an array of integers.}
\item{map}{  Map of the SNPs' locations on the genome.}
\item{chromEnds}{The index of the last SNP on each chromosome, indexed by chromosomes.}
\item{cumEnds}{Positions of the last SNP of each chromosome on the genome, needed for plotting.}
\item{genotypeMatrix}{A matrix storing all the genotypes for the cohort. A matrix indexed by individuals and loci.}
\item{probands}{Indices of the individuals present in the cohort, in the same order as genotypeMatrix.}
\item{idNames}{Array with identifiers of samples, in the same order as genotypeMatrix.}
\item{pedigree}{Pedigree for the cohort - rows of individuals whose parents are specified.}
\item{chromList}{A 1-D representation of genotypes as integers; a list indexed by chromosomes.}
\item{geneticSnpMap}{Position of each SNP in centiMorgans, starting from 0 at each chromosome.}
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
\keyword{data}
\keyword{loading}
