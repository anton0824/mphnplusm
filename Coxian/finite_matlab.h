#ifndef _FINITE_
#define _FINITE_

#include <assert.h>
#include "generator.h"
#include "integral.h"
#include "gaussianReference.h"
#include "index.h"

class finiteElement {

 private:
  vector<vector<double> > grid;  
  int dim;
  int degree;

 public:
  finiteElement(){};
  finiteElement(const vector<vector<double> >& grid_, const int degree_){ 
    grid = grid_;
    dim = grid.size();
    degree = degree_;
    assert( degree == 4|| degree== 8 || degree==16 || degree == 32 || degree == 64
	    || degree == 128|| degree== 256 || degree== 512 || degree == 1024);
  }

  void getAb(generator&, gaussianReference&, double *, double *);
  void getAb(generator&, gaussianReference&, double *, double *, double *, double *);

};

void finiteElement::getAb(generator& g, 
			  gaussianReference & reference, 
			  double * A,
			  double * b ) { 
  
  vector<int> i_index1(dim),  r_index1(dim), i_index2(dim), r_index2(dim);
  vector<double> lo(dim), hi(dim);


  assert( dim ==2 ); // the following code is for 2d
  
  int numberOfInteriorPoints0=grid[0].size()-2;
  int numberOfInteriorPoints1=grid[1].size()-2;
  

  integral a(grid, degree);
  gridIndex gi(grid);
  int N= gi.getN();

  for (int ki0=1; ki0<=numberOfInteriorPoints0; ki0++)
    for (int kr0=0; kr0<2; kr0++) 
      for (int ki1=1; ki1<=numberOfInteriorPoints1; ki1++)
	for (int kr1=0; kr1<2; kr1++)
	  {
	    i_index1[0]=ki0;
	    i_index1[1]=ki1;
	    r_index1[0]=kr0;
	    r_index1[1]=kr1;

	    int k = gi.getIndex(i_index1, r_index1);

	    b[k] = a.inner_product_rhs(g, reference, i_index1, r_index1, lo, hi);

	    for (int ji0=max(ki0-1,1); 
		 ji0<=min(ki0+1,numberOfInteriorPoints0); ji0++) 
	      for (int ji1=max(ki1-1,1); 
		 ji1<=min(ki1+1,numberOfInteriorPoints1); ji1++) 
		for (int jr0=0; jr0<2; jr0++) 
		  for (int jr1=0; jr1<2; jr1++) 
		    {
		      
		      i_index2[0]=ji0;
		      i_index2[1]=ji1;
		      r_index2[0]=jr0;
		      r_index2[1]=jr1;
		      int j=gi.getIndex(i_index2, r_index2);

			  if (k<=j)
			  {
				  A[j*N+k]= a.inner_product(g, reference, i_index1, r_index1,
					       i_index2, r_index2, lo, hi);
				  if (k<j) A[k*N+j]=A[j*N+k];
			  }
		    }
	  }

  
}

void finiteElement::getAb(generator& g, 
			  gaussianReference & reference, 
			  double * Ar,
			  double * Ac,
			  double * Av,
			  double * b ) { 
  
  vector<int> i_index1(dim),  r_index1(dim), i_index2(dim), r_index2(dim);
  vector<double> lo(dim), hi(dim);


  assert( dim ==2 ); // the following code is for 2d
  
  int numberOfInteriorPoints0=grid[0].size()-2;
  int numberOfInteriorPoints1=grid[1].size()-2;
  int n=0;
  

  integral a(grid, degree);
  gridIndex gi(grid);
  int N= gi.getN();

  for (int ki0=1; ki0<=numberOfInteriorPoints0; ki0++)
    for (int kr0=0; kr0<2; kr0++) 
      for (int ki1=1; ki1<=numberOfInteriorPoints1; ki1++)
	for (int kr1=0; kr1<2; kr1++)
	  {
	    i_index1[0]=ki0;
	    i_index1[1]=ki1;
	    r_index1[0]=kr0;
	    r_index1[1]=kr1;

	    int k = gi.getIndex(i_index1, r_index1);

	    b[k] = a.inner_product_rhs(g, reference, i_index1, r_index1, lo, hi);

	    for (int ji0=max(ki0-1,1); 
		 ji0<=min(ki0+1,numberOfInteriorPoints0); ji0++) 
	      for (int ji1=max(ki1-1,1); 
		 ji1<=min(ki1+1,numberOfInteriorPoints1); ji1++) 
		for (int jr0=0; jr0<2; jr0++) 
		  for (int jr1=0; jr1<2; jr1++) 
		    {
		      
		      i_index2[0]=ji0;
		      i_index2[1]=ji1;
		      r_index2[0]=jr0;
		      r_index2[1]=jr1;
		      int j=gi.getIndex(i_index2, r_index2);

			  if (k<=j)
			  {
				  Ar[n]=j+1; 
				  Ac[n]=k+1;
				  Av[n]= a.inner_product(g, reference, i_index1, r_index1,
					       i_index2, r_index2, lo, hi);
				  n++;
				  if (k<j)
				  {
					  Ar[n]=k+1; 
					  Ac[n]=j+1; 
					  Av[n]=Av[n-1]; 
					  n++;
				  }
			  }
		    }
	  }  
}

#endif
