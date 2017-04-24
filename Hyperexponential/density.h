#ifndef _DENSITY_
#define _DENSITY_
#include "generator.h"
#include <assert.h>
#include "index.h"
#include <iostream>
#include "gaussianReference.h"

class density {
 private:
  vector<vector<double> > grid;
  int numberWithinGrid;
  vector<int> meshSize;
  int totalMeshSize;
 public:

  density() {};

  density(const vector<vector<double> >& grid_,
	  const int numberWithinGrid_) {
    grid = grid_;
    numberWithinGrid = numberWithinGrid_;
    assert( numberWithinGrid>=1);
    meshSize.resize(2);
    for (int i=0; i<2; i++)
      meshSize[i] = grid[i].size()+(grid[i].size()-1)*(numberWithinGrid-1);
    totalMeshSize = meshSize[0]*meshSize[1];
  }

  void compute_density(generator&, gaussianReference&, const double *, double * p); 

  int meshIndex(const int i0, const int i1, const int j, const int k ) {

    return (i0*numberWithinGrid)*meshSize[0] 
      + j*meshSize[0] 
      + i1* numberWithinGrid 
      + k;

  }
  int getTotalMeshSize() { return totalMeshSize;}

};


void density::compute_density(generator& g,
			      gaussianReference & reference,
			      const double *u,
			      double * p) {

  int n = numberWithinGrid;

  vector<int> i_index(2), r_index(2);
  vector<double> y(2);
  int ell;
  gridIndex gi(grid);
  unsigned i0;
  unsigned i1;


  for ( i0=1; i0< grid[0].size()-2; i0++)
    for ( i1=1; i1< grid[1].size()-2; i1++)
      {
	// interior
	for (int j=0; j<n; j++ )
	  for (int k=0; k<n; k++)
	    {
	      
	      y[0] = grid[0][i0] + j*(grid[0][i0+1]- grid[0][i0])/n;
	      y[1] = grid[1][i1] + k*(grid[1][i1+1]- grid[1][i1])/n;

	      ell = meshIndex(i0, i1, j, k);
	      p[0*totalMeshSize+ell]= y[0];
	      p[1*totalMeshSize+ell]= y[1];
	      p[2*totalMeshSize+ell]= reference(y);
	    
	      for (int m0=0; m0<2; m0++) 
		for (int m1=0; m1<2; m1++) 
		  {
		    i_index[0] = i0+m0;
		    i_index[1] = i1+m1;
		    for (int r0=0; r0<2; r0++)
		      for (int r1=0; r1<2; r1++)
			{

			  r_index[0]=r0;
			  r_index[1]=r1;
			  p[2*totalMeshSize+ell] -= u[gi.getIndex(i_index, r_index)]
			    *g.Lf(i_index, r_index, y)*reference(y); 
			} // r0, r1
		  } // m0, m1
	    } // j, k
      } // i0, i1


  for ( i0=1; i0< grid[0].size()-2; i0++)
    {
      i1 = 0; // bottom
      for (int j=0; j<n; j++ )
	for (int k=0; k<n; k++)
	  {
	      
	    y[0] = grid[0][i0] + j*(grid[0][i0+1]- grid[0][i0])/n;
	    y[1] = grid[1][i1] + k*(grid[1][i1+1]- grid[1][i1])/n;

	    ell = meshIndex(i0, i1, j, k);
	    p[0*totalMeshSize+ell]= y[0];
	    p[1*totalMeshSize+ell]= y[1];
	    p[2*totalMeshSize+ell]= reference(y);
	    
	    for (int m0=0; m0<2; m0++) 
		{
		  i_index[0] = i0+m0;
		  i_index[1] = i1+1;
		  for (int r0=0; r0<2; r0++)
		    for (int r1=0; r1<2; r1++)
		      {
			r_index[0]=r0;
			r_index[1]=r1;
			p[2*totalMeshSize+ell] -= u[gi.getIndex(i_index, r_index)]
			  *g.Lf(i_index, r_index, y)*reference(y); 
			} // r0, r1
		  } // m0, 

	    } // j, k
      } // i0



  for ( i0=1; i0< grid[0].size()-2; i0++)
    {
      i1 = grid[1].size()-2; // top
      for (int j=0; j<n; j++ )
	for (int k=0; k<=n; k++) // top
	  {
	      
	    y[0] = grid[0][i0] + j*(grid[0][i0+1]- grid[0][i0])/n;
	    y[1] = grid[1][i1] + k*(grid[1][i1+1]- grid[1][i1])/n;

	    ell = meshIndex(i0, i1, j, k);
	    p[0*totalMeshSize+ell]= y[0];
	    p[1*totalMeshSize+ell]= y[1];
	    p[2*totalMeshSize+ell]= reference(y);
	    
	    for (int m0=0; m0<2; m0++) 
		{
		  i_index[0] = i0+m0;
		  i_index[1] = i1;
		  for (int r0=0; r0<2; r0++)
		    for (int r1=0; r1<2; r1++)
		      {
			r_index[0]=r0;
			r_index[1]=r1;
			p[2*totalMeshSize+ell] -= u[gi.getIndex(i_index, r_index)]
			  *g.Lf(i_index, r_index, y)*reference(y); 
			} // r0, r1
		  } // m0, 
	    } // j, k
      } // i0

  for ( i1=1; i1< grid[1].size()-2; i1++)
    {
      i0=0; // left
      for (int j=0; j<n; j++ )
	for (int k=0; k<n; k++) 
	  {
	      
	    y[0] = grid[0][i0] + j*(grid[0][i0+1]- grid[0][i0])/n;
	    y[1] = grid[1][i1] + k*(grid[1][i1+1]- grid[1][i1])/n;

	    ell = meshIndex(i0, i1, j, k);
	    p[0*totalMeshSize+ell]= y[0];
	    p[1*totalMeshSize+ell]= y[1];
	    p[2*totalMeshSize+ell]= reference(y);
	    
	    for (int m1=0; m1<2; m1++) 
		  {
		    i_index[0] = i0+1; // left
		    i_index[1] = i1+m1;
		    for (int r0=0; r0<2; r0++)
		      for (int r1=0; r1<2; r1++)
			{

			  r_index[0]=r0;
			  r_index[1]=r1;
			  p[2*totalMeshSize+ell] -= u[gi.getIndex(i_index, r_index)]
			    *g.Lf(i_index, r_index, y)*reference(y); 
			} // r0, r1
		  } // m0, m1

	    } // j, k
      } // i0, i1
  


  for ( i1=1; i1< grid[1].size()-2; i1++)
    {
      i0=grid[0].size()-2; // right
      for (int j=0; j<=n; j++ ) // right
	for (int k=0; k<n; k++)
	  {
	      
	    y[0] = grid[0][i0] + j*(grid[0][i0+1]- grid[0][i0])/n;
	    y[1] = grid[1][i1] + k*(grid[1][i1+1]- grid[1][i1])/n;

	    ell = meshIndex(i0, i1, j, k);
	    p[0*totalMeshSize+ell]= y[0];
	    p[1*totalMeshSize+ell]= y[1];
	    p[2*totalMeshSize+ell]= reference(y);
	    
	    for (int m1=0; m1<2; m1++) 
		  {
		    i_index[0] = i0; // right
		    i_index[1] = i1+m1;
		    for (int r0=0; r0<2; r0++)
		      for (int r1=0; r1<2; r1++)
			{

			  r_index[0]=r0;
			  r_index[1]=r1;
			  p[2*totalMeshSize+ell] -= u[gi.getIndex(i_index, r_index)]
			    *g.Lf(i_index, r_index, y)*reference(y); 
			} // r0, r1
		  } // m0, m1
	    } // j, k
      } // i0, i1

  i0 = 0;
  i1 = 0; // southWest
  for (int j=0; j< n; j++ ) 
    for (int k=0; k<n; k++)
      {
	      
	y[0] = grid[0][i0] + j*(grid[0][i0+1]- grid[0][i0])/n;
	y[1] = grid[1][i1] + k*(grid[1][i1+1]- grid[1][i1])/n;

	ell = meshIndex(i0, i1, j, k);
	p[0*totalMeshSize+ell]= y[0];
	p[1*totalMeshSize+ell]= y[1];
	p[2*totalMeshSize+ell]= reference(y);
	    
	
	i_index[0] = i0+1; 
	i_index[1] = i1+1; // southwest
	for (int r0=0; r0<2; r0++)
	  for (int r1=0; r1<2; r1++)
	    {
	      r_index[0]=r0;
	      r_index[1]=r1;
	      p[2*totalMeshSize+ell] -= u[gi.getIndex(i_index, r_index)]
		*g.Lf(i_index, r_index, y)*reference(y); 
	    } // r0, r1
      } // j, k
  
  i0= 0;
  i1 =grid[1].size()-2; // northWest
  for (int j=0; j< n; j++ ) 
    for (int k=0; k<= n; k++) // notice the <= sign
      {
	      
	y[0] = grid[0][i0] + j*(grid[0][i0+1]- grid[0][i0])/n;
	y[1] = grid[1][i1] + k*(grid[1][i1+1]- grid[1][i1])/n;

	ell = meshIndex(i0, i1, j, k);
	p[0*totalMeshSize+ell]= y[0];
	p[1*totalMeshSize+ell]= y[1];
	p[2*totalMeshSize+ell]= reference(y);
	    
	
	i_index[0] = i0+1; 
	i_index[1] = i1; // northwest
	for (int r0=0; r0<2; r0++)
	  for (int r1=0; r1<2; r1++)
	    {
	      r_index[0]=r0;
	      r_index[1]=r1;
	      p[2*totalMeshSize+ell] -= u[gi.getIndex(i_index, r_index)]
		*g.Lf(i_index, r_index, y)*reference(y); 
	    } // r0, r1
      } // j, k
  

  i0= grid[0].size()-2;
  i1 =0; // southEast
  for (int j=0; j<= n; j++ )  // note the <=n sign
    for (int k=0; k<n; k++)
      {
	      
	y[0] = grid[0][i0] + j*(grid[0][i0+1]- grid[0][i0])/n;
	y[1] = grid[1][i1] + k*(grid[1][i1+1]- grid[1][i1])/n;

	ell = meshIndex(i0, i1, j, k);
	p[0*totalMeshSize+ell]= y[0];
	p[1*totalMeshSize+ell]= y[1];
	p[2*totalMeshSize+ell]= reference(y);
	    
	
	i_index[0] = i0; 
	i_index[1] = i1+1; // southwest
	for (int r0=0; r0<2; r0++)
	  for (int r1=0; r1<2; r1++)
	    {
	      r_index[0]=r0;
	      r_index[1]=r1;
	      p[2*totalMeshSize+ell] -= u[gi.getIndex(i_index, r_index)]
		*g.Lf(i_index, r_index, y)*reference(y); 
	    } // r0, r1
      } // j, k

  i0 = grid[0].size()-2;
  i1 = grid[0].size()-2; // NorthEast
  for (int j=0; j<= n; j++ )  // note the <= sign
    for (int k=0; k<=n; k++) // note the <= sign
      {
	      
	y[0] = grid[0][i0] + j*(grid[0][i0+1]- grid[0][i0])/n;
	y[1] = grid[1][i1] + k*(grid[1][i1+1]- grid[1][i1])/n;

	ell = meshIndex(i0, i1, j, k);
	p[0*totalMeshSize+ell]= y[0];
	p[1*totalMeshSize+ell]= y[1];
	p[2*totalMeshSize+ell]= reference(y);
	    
	
	i_index[0] = i0; 
	i_index[1] = i1; // NorthEast
	for (int r0=0; r0<2; r0++)
	  for (int r1=0; r1<2; r1++)
	    {
	      r_index[0]=r0;
	      r_index[1]=r1;
	      p[2*totalMeshSize+ell] -= u[gi.getIndex(i_index, r_index)]
		*g.Lf(i_index, r_index, y)*reference(y); 
	    } // r0, r1
      } // j, k
}


#endif
