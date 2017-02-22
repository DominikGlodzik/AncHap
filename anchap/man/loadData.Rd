\name{loadData}
\alias{loadData}
\title{
Loading data required by anchap.
}
\description{
Loads all the data for anchap.
}
\usage{
loadData(mapFile, genotypesFile, pedigreeFilepath, geneticSnpMapFile = NULL)
}
\arguments{
  \item{mapFile}{
  Path for the file storing the mao of all the SNPs, in the PLINK .map format.
}
  \item{genotypesFile}{
  Paths to Plink binary files (bim, bed, fam) with the genotypes for the cohort.
}
  \item{pedigreeFilepath}{
  Path to a pedigree file for the cohort, given as an xsl spread-sheet.
}
  \item{geneticSnpMapFile}{
  Path to a genetic map. A map for each chromosome is typically in separate files. This argument takes a form of a 2-element list, with first element of the list being the prefix of the files, and the second element is to postfix. Between these two parts, the files should have integer 1-22.
}
}
\details{
}
\value{
All the data for a cohort required by anchap.
item{iGene}{Genotypes of individuals in the cohort. A list indexed by individuals, whose entry contains an individual's genotype, represented by an array of integers.}
\item{map}{  Map of the SNPs' locations on the genome.}
\item{chromEnds}{The index of the last SNP on each chromosome, indexed by chromosomes.}
\item{cumEnds}{Positions of the last SNP of each chromosome on the genome, needed for plotting.}
\item{genotypeMatrix}{A matrix storing all the genotypes for the cohort. A matrix indexed by individuals and loci.}
\item{probands}{Indices of the individuals present in the cohort.}
\item{pedigree}{Pedigree for the cohort - rows of individuals whose parents are specified.}
\item{chromList}{A 1-D representation of genotypes as integers; a list indexed by chromosomes.}
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
\keyword{ data}
\keyword{loading}
