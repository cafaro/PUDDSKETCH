/**
 * @file Quantiles.cc
 * @author CM
 */

#include "Quantiles.h"

// double QQ[] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};
// int nQ = sizeof(QQ)/sizeof(double);

extern double QQ[];
extern int nQ;



double quickselect(double *data, int len, int pos) {
  
  int i, ir, j, l, mid;
  double a, temp;

  l=0;
  ir=len-1;

  for(;;) {
    if (ir <= l+1) { 
      if (ir == l+1 && data[ir] < data[l]) {
	    SWAP(data[l],data[ir]);
      }
    return data[pos];
    }
    else 
    {
      mid=(l+ir) >> 1; 
      SWAP(data[mid],data[l+1]);
      if (data[l] > data[ir]) {
	    SWAP(data[l],data[ir]);
      }
      if (data[l+1] > data[ir]) {
	    SWAP(data[l+1],data[ir]);
      }
      if (data[l] > data[l+1]) {
	    SWAP(data[l],data[l+1]);
      }
      i=l+1; 
      j=ir;
      a=data[l+1];

      for (;;) { 
	    do i++; while (data[i] < a); 
	    do j--; while (data[j] > a); 
	    if (j < i) break; 
	    SWAP(data[i],data[j]);
      } 
      data[l+1]=data[j]; 
      data[j]=a;
      if (j >= pos) ir=j-1; 
      if (j <= pos) l=i;
    }
  }
}



double *getExactQuantiles(double *gdataset, long dlen, int *eLen, FILE *fp) {
  
  double *exact = (double *) malloc(nQ * sizeof(double));
  *eLen = nQ;
    

  if (fp) {
    fprintf(fp, "Exact Quantiles\n");
  }//fi
    
  int rank = 0;
  for (int c = 0; c < nQ; ++c) {
        
    rank = std::floor(QQ[c]*(dlen-1) + 1) - 1;
    exact[c] = quickselect(gdataset, dlen, rank);


    if (fp) {
      fprintf(fp, "q: %.3f(%d) \t Estimate: %.16f\n", QQ[c], rank, exact[c]);
    } else {
      fprintf(stderr,  "q: %.3f(%d) \t Estimate: %.16f\n", QQ[c], rank, exact[c]);
    }//fi fp

  }//for c

  if (fp) {
    fprintf(fp, "---------------------------------------------\n");
  }//fi

  return exact;
}
