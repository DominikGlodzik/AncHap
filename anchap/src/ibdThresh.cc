#include <iostream>
#include <cmath>
using namespace std;

// compile with R CMD SHLIB ibdHapGen.cc

extern "C" {

  // id1 is the haplotype
  // id2 is the genotype
  int ibdThresh (int *n, int *id1, int *id2,  int *lThresh, int* noSeqs, int* starts, int* ends,  int* chromends ) {

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
      
      seqOn = 0; 
      nonMissingLoci = 0;


      for (int pointer=firstSnp; pointer<lastSnp; pointer++) {	   
	     
	if ((((*id1==1)&&(*id2==3))||((*id1==3)&&(*id2==1)))&&(seqOn>0)) {
	  // opposing homozygote
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
      
      // check if there is an IBD sequence at the end of the region
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
