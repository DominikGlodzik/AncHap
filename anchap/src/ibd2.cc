#include <iostream>
#include <cmath>
using namespace std;

// compile with R CMD SHLIB ibd2.cc

extern "C" {

  // id1 is the haplotype
  // id2 is the genotype
  int ibd2 (int *n, int *id1, int *id2,  int *lThresh, int* noSeqs, int* starts, int* ends,  int* chromends ) {

    // location of previous opposing homozygote
    int seqStart;
    int seqOn;

    // border SNP ids on respective chromosomes
    int firstSnp;
    int lastSnp;

    // count of non-missing loci in the IBD region
    int nonMissingLoci;

    // loop over chromosomes
    for (int chrom = 1; chrom<=22; chrom++) {     
      
      if (chrom==1) {firstSnp = 1;}
      else {firstSnp = chromends[chrom - 2] + 1 ;}
      lastSnp = chromends[chrom - 1];
      seqStart = firstSnp;
      
      seqOn = 0; 
      nonMissingLoci = 0;

      // loop over SNPs within one chromosome
      for (int pointer=firstSnp; pointer<lastSnp; pointer++) {	   
	     
	// non IBD2 - genotypes non missing, and different from each other

	if ((seqOn>0)&&(*id1!=0)&&(*id2!=0)&&(*id1!=*id2)) {
	  // non ibd matching
	  if ((pointer-seqStart)>=*lThresh) {
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
	  // start the sequence, if not started already
	  
	  // AA found on haplotype / genotype 1
	  if (*id1==1) {
	    if (seqOn==0) {seqStart = pointer; seqOn = 1;}
	  }
	// BB found on haplotype / genotype 2
	  else if (*id1==3) {
	    if (seqOn==0) {seqStart = pointer; seqOn = 1;}
	  }
	}

	id1++; id2++;
      }
      
      // finished the chromosome, check if there is an IBD sequence at the end of the region
      if ((seqOn>0)&&((lastSnp - seqStart)>=*lThresh)) {
	// found an ibd sequence
	*starts = seqStart;
	*ends = lastSnp;
	starts++; ends++;
	*noSeqs= *noSeqs + 1;
      }
    } // end of the chromosome loop  
  } // ibdThresh
} // extern "C"
