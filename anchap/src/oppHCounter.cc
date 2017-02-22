     
     #include <iostream>
     #include <cmath>

     extern "C" {
     
       int oppHCounter (int *n, int *id1, int *id2, int* count ) {

	 *count = 0;

	for (int pointer=1; pointer<=*n; pointer++) {	   
 
	  if (((*id1==1)&&(*id2==3))||((*id1==3)&&(*id2==1))) {
	    // opposing homozygote
	    *count = *count + 1;
	  }

	  id1++;
	  id2++;
	}
   
      }
     
     } // extern "C"
