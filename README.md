# AncHap
Code for analysis of identity-by-descent in genetic data from isolated populations.

Packages required:
snpStats
to install, from R, run   
	source("http://bioconductor.org/biocLite.R")
	biocLite("snpStats") 

copy anchap_2.2.tar.gz to the work folder
from the command line: R CMD install anchap_2.2.tar.gz 

open R session
library(anchap)

run demos on the toy data set
demo('toyData')
demo('chr2')

and run the attached example:
source('runme.R')

