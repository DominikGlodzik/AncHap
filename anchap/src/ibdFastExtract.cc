// does not care about the chromosomes


     #include <iostream>
     #include <cmath>

     extern "C" {
     
       int ibdFastExtract (int *n, int *id1, int *id2, int *lThresh, int* noSeqs, int* starts, int* ends ) {

	 
	 int pointer = 1;
	 int seqStart = 1;
	 
	 int snpDiff = -1;
	 int seqLength = -1;

	for (int pointer=1; pointer<=*n; pointer++) {	   
 
	  if (((*id1==1)&&(*id2==3))||((*id1==3)&&(*id2==1))) {
	      // opposing homozygote

              seqLength = pointer-seqStart;   	      
	      if (seqLength>=*lThresh) {
		// found a sequence whose length is above the threshold
		*starts = seqStart;
		*ends = pointer;
		starts++; ends++; *noSeqs= *noSeqs + 1;
	      }                     
	      seqStart = pointer + 1; 
	    }
	  id1++;
	  id2++;
	}
   
	seqLength = *n-seqStart; 
	if (seqLength>=*lThresh) {
	  // found a sequence whose length is above the threshold
	  *starts = seqStart;
	  *ends = *n+1;
	  starts++; ends++; *noSeqs= *noSeqs + 1;
	}    

      }
     
     } // extern "C"
