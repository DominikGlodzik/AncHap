#include <iostream>
#include <cmath>
using namespace std;

// extracts IBD shared sequences, keeping track of the missing values

// compile with R CMD SHLIB ibdHapGen.cc

extern "C" {

  // id1 is the haplotype
  // id2 is the genotype
  int ibdThreshMiss (int *n, int *id1, int *id2,  int *lThresh, int* noSeqs, int* starts, int* ends,  int* chromends ) {

    // location of previous opposing homozygote
    int seqStart;
    int seqOn;

    // border SNP ids on respective chromosomes
    int firstSnp;
    int lastSnp;

    // count of non-missing loci in the IBD region
    int nonMissingLoci;

    for (int chrom = 1; chrom<=22; chrom++) {     
      
      if (chrom==1) {firstSnp = 1;}
      else {firstSnp = chromends[chrom - 2] + 1 ;}
      lastSnp = chromends[chrom - 1];
      seqStart = firstSnp;
      
      // whether we are currently traversing a sequence without opposing homozygotes
      seqOn = 0; 
      // number of loci not missing on either genotype, since the last opposing homozygote
      nonMissingLoci = 0;

      for (int pointer=firstSnp; pointer<=lastSnp; pointer++) {	   
	     
	if ((((*id1==1)&&(*id2==3))||((*id1==3)&&(*id2==1)))&&(seqOn>0)) {
	  // opposing homozygote
	  if (nonMissingLoci>=*lThresh) {
	    // exceeded a threshold length - found an ibd sequence      
	    *starts = seqStart;
	    *ends = pointer;
	    *noSeqs= *noSeqs + 1;	
	    
	    starts++; ends++;
	  }
	  seqOn = 0;
	  nonMissingLoci = 0;
	}      
	else {
	  // non-opposing homozygote

	  // check if both genotypes not missing
	  if ((*id1!=0)&&(*id2!=0)) {
	    if (seqOn==0) {seqStart = pointer+1; seqOn = 1;}
	    nonMissingLoci++;
	  }
	}

	id1++; id2++;
      }
      
      // check if there is an IBD sequence at the end of the region
      if ((seqOn>0)&&(nonMissingLoci>=*lThresh)) {
	// found an ibd sequence
	*starts = seqStart;
	*ends = lastSnp;
	starts++; ends++;
	*noSeqs= *noSeqs + 1;
      }
    } // end of the chromosome loop  
  } // ibdThresh
} // extern "C"
