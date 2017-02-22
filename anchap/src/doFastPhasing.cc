#include <iostream>
#include <cmath>
#include <stdio.h>
using namespace std;

// compile with R CMD SHLIB doFastPhasing.cc

extern "C" {

  int getIndex(int indId, int snp, int noSnps, int firstSnp) {
    // assuming the indices start at 0
    // ie. R's index 1 is 0
    //         index length(vec) is length(vec)-1
    return ((indId-1)*noSnps + (snp) - (firstSnp));
  }

  int getOppGene(int gene) {
    int value = 0;
    if (gene==1) {value=3;}
    else if (gene==3) {value=1;}
    else if (gene==4) {value=4;}
    return value;
  }

    int h1SeqsAtLocus = 0;
    int h2SeqsAtLocus = 0;
    int h1HetSeqsAtLocus = 0;
    int h2HetSeqsAtLocus = 0;


    int h1ones = 0;
    int h1threes = 0;
    int h2ones = 0;
    int h2threes = 0;
    int modeH1 = 0;
    int modeH2 = 0;

    int countIncH1 = 0;
    int countIncH2 = 0;

    int probandsGene = 0;
    int allele = 0;
    int indId = 0;
    
    int heteroLoci = 0;
    int heteroPhased = 0;

    int vecI = 0;

  // id1 is the haplotype
  // id2 is the genotype
  int doFastPhasing (int *noSnps, int* firstSnpChrom, int* lastSnpChrom, int* firstSnpAnalysed, int* lastSnpAnalysed, int *proband, int *noIndividiuals,  int *genes, int *noH1Seqs, double* h1Starts, double *h1Ends, int *h1Ids, int *noH2Seqs, double* h2Starts, double *h2Ends, int *h2Ids, int *noH1atLocus, int *noH2atLocus,  int *haplotype1, int *haplotype2, int *majority1, int *majority2, int *genotype, int *inconsH1, int *inconsH2, int *incons, int *pH) {

    //cout << "no Snps " << *noSnps << endl;
    //cout << "proband " << *proband << endl;
    //cout << "no individuals " <<  *noIndividiuals << endl;
    //cout << "input " << *firstSnpChrom << " "<<*lastSnpChrom << " " << *firstSnpAnalysed << " "<<*lastSnpAnalysed << endl;
  
    h1SeqsAtLocus = 0;
    h2SeqsAtLocus = 0;
    h1HetSeqsAtLocus = 0;
    h2HetSeqsAtLocus = 0;
   
    h1ones = 0;
    h1threes = 0;
    h2ones = 0;
    h2threes = 0;
    modeH1 = 0;
    modeH2 = 0;

    countIncH1 = 0;
    countIncH2 = 0;

    probandsGene = 0;
    allele = 0;
    indId = 0;
    
    heteroLoci = 0;
    heteroPhased = 0;

    vecI = 0;

    for (int indx=*firstSnpAnalysed; indx<=*lastSnpAnalysed; indx++) {      
      

      modeH1 = 0;
      modeH2 = 0;
      countIncH1 = 0;
      countIncH2 = 0;

      probandsGene = genes[getIndex(*proband, indx, *noSnps, *firstSnpChrom)];
      genotype[indx - *firstSnpAnalysed]=probandsGene;

      
      h1SeqsAtLocus = 0;
      h1HetSeqsAtLocus = 0;
      h1ones = 0;
      h1threes = 0;

        // H1: loop over all sequences at the locus
      for (int indx1=0; indx1<*noH1Seqs; indx1++) {
	if ((h1Starts[indx1]<=indx)&&(h1Ends[indx1]>indx)) 	 
	    { indId = h1Ids[indx1];
	      vecI = getIndex(indId, indx, *noSnps, *firstSnpChrom);
	      allele = genes[vecI];	      
	      
	      if (allele==1) {h1ones++; h1SeqsAtLocus++;}
	      else if (allele==3) {h1threes++; h1SeqsAtLocus++;}
	      else {h1HetSeqsAtLocus++; }
	      	    }	  
	}
      noH1atLocus[indx - *firstSnpAnalysed] = h1SeqsAtLocus + h1HetSeqsAtLocus;
      

      // H1: find the most frequent homozygous allele
      if ((h1ones+h1threes)>0) {
	if (h1ones>h1threes) { modeH1 =1; countIncH1 = h1threes;}
	else if (h1ones<h1threes) { modeH1 = 3; countIncH1 = h1ones;}
	else if (h1ones==h1threes) { modeH1 = 1; countIncH1 = h1threes;}
      } 

      majority1[indx - *firstSnpAnalysed] = modeH1;
      inconsH1[indx - *firstSnpAnalysed] = countIncH1;

      h2SeqsAtLocus = 0;
      h2HetSeqsAtLocus = 0;
      h2ones = 0;
      h2threes = 0;

      // H2: loop over all sequences at the locus
      for (int indx2=0; indx2<*noH2Seqs; indx2++) {
	if ((h2Starts[indx2]<=indx)&&(h2Ends[indx2]>indx)) 	 
	  { indId = h2Ids[indx2];
	    allele = genes[getIndex(indId, indx, *noSnps, *firstSnpChrom)];
	    //cout << "id "<< indId << " allele " <<allele << endl;	    
	    if (allele==1) {h2ones++; h2SeqsAtLocus++;}
	    else if (allele==3) {h2threes++; h2SeqsAtLocus++;}
	    else { h2HetSeqsAtLocus++;}
	      
	  }	 
      }
      noH2atLocus[indx - *firstSnpAnalysed] = h2SeqsAtLocus +  h2HetSeqsAtLocus;
      
      if ((h2ones+h2threes)>0) {
	if (h2ones>h2threes) { modeH2=1; countIncH2 = h2threes;}
	  else if (h2ones<h2threes) { modeH2=3; countIncH2 = h2ones;}
	  else if (h2ones==h2threes) { modeH2 = 1; countIncH2 = h2threes;}
      }
      majority2[indx - *firstSnpAnalysed] = modeH2;
      inconsH2[indx - *firstSnpAnalysed] = countIncH2;

      // sum all parent homozygotes
      pH[indx - *firstSnpAnalysed] =  h1threes + h2ones + h1ones + h2threes;



      if (probandsGene==1) { haplotype1[indx - *firstSnpAnalysed]=1; haplotype2[indx - *firstSnpAnalysed]=1; }
      else if (probandsGene==3) { haplotype1[indx - *firstSnpAnalysed]=3; haplotype2[indx - *firstSnpAnalysed]=3; }
      else {
	
	heteroLoci++;
	
	if ((h1SeqsAtLocus+h2SeqsAtLocus)>0) {
	  heteroPhased++;
 
	  
	  if (h1SeqsAtLocus>=h2SeqsAtLocus) {
	    haplotype1[indx - *firstSnpAnalysed] = modeH1;
	    haplotype2[indx - *firstSnpAnalysed] = getOppGene(modeH1);

	    //total count of inconsistencies
	    
	    if (modeH1==1) {
	      // mode H1 1
	      incons[indx - *firstSnpAnalysed] = h1threes + h2ones;
	    }
	    else {
	      // mode H1 3
	      incons[indx - *firstSnpAnalysed] = h1ones + h2threes;
	    }

	  }
	  else {
	    // H2: find the most frequent homozygous allele
	    haplotype2[indx - *firstSnpAnalysed] = modeH2;
	    haplotype1[indx - *firstSnpAnalysed] = getOppGene(modeH2);

	    //total count of inconsistencies
	    //count of minority alleles on H2 + count of H2 majority on H1
	    //incons[indx - *firstSnpAnalysed] = countIncH2;
	    if (modeH2==1) {
	      // mode H2 1
	      incons[indx - *firstSnpAnalysed] = h2threes + h1ones;
	    }
	    else {
	      // mode H2 3
	      incons[indx - *firstSnpAnalysed] = h2ones + h1threes;
	    }
	  }
	  //end else
	}
	//end else
      }
      // end else
    }
    // end for loop over all loci
    //cout << "hetero loci " <<  heteroLoci << endl;
    //cout << "hetero loci phased " <<  heteroPhased << endl;
  }

  
  
} 
// extern "C"
