#ifndef _INTEGRAL_
#define _INTEGRAL_
#include "gauss_legendre.h"
#include <vector>
#include <assert.h>
#include "generator.h"
using namespace std;

class integral {

 private:
  int dim;
  vector<vector<double> > grid;
  vector<double *> x, w; 
  vector<double> y0, y1, y2, y3; // dummy variable used in integration
  vector<double> C, D, Cx;
  vector<int> weightPoints;

 public:
  
  integral() {};

  integral(vector<vector<double> >& grid_, const int degree) {
    grid= grid_;
    dim = grid.size();
    assert(degree>=2 && degree <=1024 && degree%2 == 0);
    y0.resize(dim);
    y1.resize(dim);
    y2.resize(dim);
    y3.resize(dim);
    C.resize(dim);
    D.resize(dim);
    Cx.resize(dim);
    x.resize(dim);
    w.resize(dim);
    weightPoints.resize(dim);

    switch (degree) {

    case 4:
      for (int i=0; i<dim; i++) {
	x[i] = x4;
	w[i]=  w4;
      }
      break;

    case 8:
      for (int i=0; i<dim; i++) {
	x[i] = x8;
	w[i]=  w8;
      }
      break;

    case 16:
      for (int i=0; i<dim; i++) {
	x[i] = x16;
	w[i]=w16;
      }
      break;

    case 32:
      for (int i=0; i<dim; i++) {
	x[i] = x32;
	w[i]=w32;
      }
      break;

    case 64:
      for (int i=0; i<dim; i++) {
	x[i] = x64;
	w[i]=w64;
      }
      break;

    case 128:
      for (int i=0; i<dim; i++) {
	x[i] = x128;
	w[i]=w128;
      }
      break;

    case 256:
      for (int i=0; i<dim; i++) {
	x[i] = x256;
	w[i] = w256;
      }
      break;

    case 512:
      for (int i=0; i<dim; i++) {
	x[i] = x512;
	w[i] = w512;
      }
      break;
      
    case 1024:
      for (int i=0; i<dim; i++) {
	x[i] = x1024;
	w[i] = w1024;
      }
      break;

    default:
      assert(0);
    }
    for (int i=0; i<dim;i++)
      weightPoints[i]=degree/2;
  }
  


 double integral_2d(generator&, gaussianReference&, 
		    const vector<int>&, const vector<int>&,
		    const vector<int>&, 
		    const vector<int>&,
		    const vector<double>& lo,
		    const vector<double>& hi);

 double integral_rhs(generator&, gaussianReference&,
		    const vector<int>&, const vector<int>&,
		    const vector<double>& lo,
		    const vector<double>& hi);

 double inner_product(generator&, gaussianReference&,
		      const vector<int>&, const vector<int>&,
		      const vector<int>&,
		      const vector<int>&,
		      vector<double>& lo,
		      vector<double>& hi);

 double inner_product_rhs(generator&, gaussianReference&,
		      const vector<int>&, const vector<int>&,
		      vector<double>& lo,
		      vector<double>& hi);
			
};

double integral::inner_product(generator& g, 
			       gaussianReference& reference,
			       const vector<int>& i_index1, 
			       const vector<int>& r_index1, 
			       const vector<int>& i_index2, 
			       const vector<int>& r_index2, 
			       vector<double>& lo,
			       vector<double>& hi ) {

  assert (dim == 2); 
  double sum =0.0;


  if (i_index1[0]== i_index2[0] && i_index1[1]== i_index2[1]) {
    for (int k=0; k<2; k++) 
      for (int ell=0; ell<2; ell++) 
	{
	  lo[0] = grid[0][i_index1[0]-1+k];
	  hi[0] = grid[0][i_index1[0]+k];
	  lo[1] = grid[1][i_index1[1]-1+ell];
	  hi[1] = grid[1][i_index1[1]+ell];
	  sum += integral_2d(g, reference, i_index1, r_index1, i_index2, r_index2, lo, hi);
	}
    return sum;
  }

  if (i_index1[0]== i_index2[0] && i_index1[1] != i_index2[1]) {
    for (int k=0; k<2; k++) {
      lo[0] = grid[0][i_index1[0]-1+k];
      hi[0] = grid[0][i_index1[0]+k];
      lo[1] = max(grid[1][i_index1[1]-1], grid[1][i_index2[1]-1]);
      hi[1] = min(grid[1][i_index1[1]+1], grid[1][i_index2[1]+1]);
      sum += integral_2d(g, reference, i_index1, r_index1, i_index2, r_index2, lo, hi);
    }
    return sum;
  }

  if (i_index1[0] != i_index2[0] && i_index1[1] == i_index2[1]) {
    for (int k=0; k<2; k++) {
      lo[0] = max(grid[0][i_index1[0]-1], grid[0][i_index2[0]-1]);
      hi[0] = min(grid[0][i_index1[0]+1], grid[0][i_index2[0]+1]);
      lo[1] = grid[1][i_index1[1]-1+k];
      hi[1] = grid[1][i_index1[1]+k];
      sum += integral_2d(g, reference, i_index1, r_index1, i_index2, r_index2, lo, hi);
    }
    return sum;
  }
 

  if (i_index1[0] != i_index2[0] && i_index1[1] != i_index2[1]) {
      lo[0] = max(grid[0][i_index1[0]-1], grid[0][i_index2[0]-1]);
      hi[0] = min(grid[0][i_index1[0]+1], grid[0][i_index2[0]+1]);
      lo[1] = max(grid[1][i_index1[1]-1], grid[1][i_index2[1]-1]);
      hi[1] = min(grid[1][i_index1[1]+1], grid[1][i_index2[1]+1]);
      sum += integral_2d(g, reference, i_index1, r_index1, i_index2, r_index2, lo, hi);
  }
  return sum;


}


double integral::integral_2d(generator& g, 
			     gaussianReference & reference,
			     const vector<int>& i_index1, 
			     const vector<int>& r_index1, 
			     const vector<int>& i_index2, 
			     const vector<int>& r_index2,
			     const vector<double>& lo,
			     const vector<double>& hi) {

  
  for (int i=0; i<dim; i++) {
    C[i] = .5*(hi[i]-lo[i]);
    D[i] = .5*(hi[i]+lo[i]);
  }

  
  double sum =0;
  
  assert (dim == 2); 

  for (int i=0; i< weightPoints[0]; i++) {
    Cx[0]= C[0]*x[0][i];
    for (int j=0; j < weightPoints[1]; j++) {
      Cx[1]= C[1]*x[1][j];

      y0[0] = D[0]+Cx[0];
      y0[1] = D[1]+Cx[1];

      y1[0] = D[0]-Cx[0];
      y1[1] = D[1]+Cx[1];

      y2[0] = D[0]+Cx[0];
      y2[1] = D[1]-Cx[1];

      y3[0] = D[0]-Cx[0];
      y3[1] = D[1]-Cx[1];

      sum += w[0][i]*w[1][j]* (
				g.Lf(i_index1, r_index1, y0)*g.Lf(i_index2, r_index2, y0) *reference(y0) +
				g.Lf(i_index1, r_index1, y1)*g.Lf(i_index2, r_index2, y1) *reference(y1) +
				g.Lf(i_index1, r_index1, y2)*g.Lf(i_index2, r_index2, y2) *reference(y2) +
				g.Lf(i_index1, r_index1, y3)*g.Lf(i_index2, r_index2, y3) *reference(y3)
				); 

    }
  }
  return sum*C[0]*C[1];
    
}


double integral::inner_product_rhs(generator& g, 
			       gaussianReference & reference,
			       const vector<int>& i_index1, 
			       const vector<int>& r_index1, 
			       vector<double>& lo,
			       vector<double>& hi ) {

  assert (dim == 2); 
  double sum =0.0;

  
    for (int k=0; k<2; k++) 
      for (int ell=0; ell<2; ell++) 
	{
	  lo[0] = grid[0][i_index1[0]-1+k];
	  hi[0] = grid[0][i_index1[0]+k];
	  lo[1] = grid[1][i_index1[1]-1+ell];
	  hi[1] = grid[1][i_index1[1]+ell];
	  sum += integral_rhs(g, reference, i_index1, r_index1, lo, hi);
	}
    return sum;
}


double integral::integral_rhs(generator& g, 
			      gaussianReference & reference,
			      const vector<int>& i_index1, 
			      const vector<int>& r_index1, 
			      const vector<double>& lo,
			      const vector<double>& hi) {

  
  for (int i=0; i<dim; i++) {
    C[i] = .5*(hi[i]-lo[i]);
    D[i] = .5*(hi[i]+lo[i]);
  }

  
  double sum =0;
  
  assert (dim == 2); 

  for (int i=0; i< weightPoints[0]; i++) {
    Cx[0]= C[0]*x[0][i];
    for (int j=0; j < weightPoints[1]; j++) {
      Cx[1]= C[1]*x[1][j];

      y0[0] = D[0]+Cx[0];
      y0[1] = D[1]+Cx[1];

      y1[0] = D[0]-Cx[0];
      y1[1] = D[1]+Cx[1];

      y2[0] = D[0]+Cx[0];
      y2[1] = D[1]-Cx[1];

      y3[0] = D[0]-Cx[0];
      y3[1] = D[1]-Cx[1];

      sum += w[0][i]*w[1][j]* (
				g.Lf(i_index1, r_index1, y0) *reference(y0) +
				g.Lf(i_index1, r_index1, y1) *reference(y1) +
				g.Lf(i_index1, r_index1, y2) *reference(y2) +
				g.Lf(i_index1, r_index1, y3) *reference(y3)
				); 

    }
  }
  return sum*C[0]*C[1];
    
}



#endif
