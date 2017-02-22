\name{doHapAssocTest}
\alias{doHapAssocTest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Statistical test of haplotypes at a locus and a phenotype
}
\description{
Given shared haplotypes at a locus, and continous value phenotypes, 
checks if haplotypes are associated with the outcome.
}
\usage{
doHapAssocTest(clusterAssignmentH1, clusterAssignmentH2, clusterSizes, phenotypes, MIN.CLUSTER.SIZE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{clusterAssignmentH1}{
  Haplotype id for haplotype 1 of all individuals. An array indexed by individuals.
}
  \item{clusterAssignmentH2}{
  Haplotype id for haplotype 2 of all individuals. An array indexed by individuals.
}
  \item{clusterSizes}{
  Frequency of each haplotype at a locus. An array indexed by haplotypes.
}
  \item{phenotypes}{
  Phenotypes for all individuals - a continous trait. An array indexed by individuals.
}
  \item{MIN.CLUSTER.SIZE}{	
  The minimum size of a cluster to be considered in the test - an integer.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
log of the P value
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
function(
                           clusterAssignmentH1,
                           clusterAssignmentH2,
                           clusterSizes,
                           phenotypes,
                           MIN.CLUSTER.SIZE) {

  clusterSizes10 <- clusterSizes[clusterSizes>=MIN.CLUSTER.SIZE]

  if (length(clusterSizes10)>0) {
 
  topClusters <- as.integer(unlist(dimnames(clusterSizes10)));

 haplotypeId1ofN <- matrix(0,length(clusterAssignmentH1), length(topClusters));
  colnames(haplotypeId1ofN ) <-  c(topClusters) 
  for (indiv in 1:(length(clusterAssignmentH1))) {
    
    hapId1 <- clusterAssignmentH1[indiv];    
    if (hapId1 \%in\% topClusters) {    haplotypeId1ofN[indiv,as.character(hapId1)] <- haplotypeId1ofN[indiv,as.character(hapId1)] + 1; }
    hapId2 <- clusterAssignmentH2[indiv];    
    if (hapId2 \%in\% topClusters) {    haplotypeId1ofN[indiv,as.character(hapId2)] <- haplotypeId1ofN[indiv,as.character(hapId2)] + 1; }
    
  }
  colnames(haplotypeId1ofN ) <- paste("h_",c(topClusters),sep="")

  dataHapPh <- data.frame(haplotypeId1ofN, phenotypes)

  xnam <- paste(colnames(haplotypeId1ofN ), sep="")
  fmla <- as.formula(paste("phenotypes ~ ", paste(xnam, collapse= "+")))

  # fit the linear model
  # use the anova table: sum of squares and degrees of freedom
  anovaTest <- anova(lm(fmla,dataHapPh))

  ## test the hypothesis that there is no
  # f test

  noTopHaps <- length(topClusters)

  
  x <- (sum(anovaTest$"Sum Sq"[1:noTopHaps])/sum(anovaTest$"Df"[1:noTopHaps]))/anovaTest$"Mean Sq"[noTopHaps+1]

  
  # get a p-Value
  pValue <- 1 - pf(x, noTopHaps, anovaTest$"Df"[1+noTopHaps])

  logP <- -log(pValue)

  logP
  }
  else {
    print("The association test failed: clusters are too small.")
  }

  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
