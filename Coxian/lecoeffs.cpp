#include "mex.h"
#include "gaussianReference.h"
//#include "h2drift.h"
#include "finite_matlab.h"
#include <vector>
#include <assert.h>
using namespace std;


/*void fE_run( const vector<vector<double> > & Gamma,
	     //hyper2drift& mu,
	  const vector<vector<double> > & grid,
	     gaussianReference & reference, 
	  double * A,
	  double * b,
	  const int degree,
	  const int nserv) {


    assert( grid.size()== Gamma.size() );
    
    finiteElement  fE(grid, degree); 
    
    generator g(Gamma, mu, grid, nserv);
    fE.getAb(g, reference, A, b);
}
*/

void fE_run( const vector<vector<double> > & Gamma,
	    // hyper2drift& mu,
	  const vector<vector<double> > & grid,
	     gaussianReference & reference, 
	  double * Ar,
	  double * Ac,
	  double * Av,
	  double * b,
	  const int degree,
	  const double alpha,
	  const double lambda,
	  const vector<double> nu,
	  const double P12,
	  const double nserv) {

    assert( grid.size()== Gamma.size() );
    
    finiteElement  fE(grid, degree); 
    
    generator g(Gamma, grid, alpha, lambda, nu[0], nu[1], P12, nserv);
    fE.getAb(g, reference, Ar, Ac, Av, b);
}


/* Extract parameters from the MATLAB function call, 
 and run fE_run with them */
void mexFunction( int nlhs, mxArray *plhs[],

				 int nrhs, const mxArray *prhs[])

{

  int m, n;
  double *data;

  if(nrhs!=10)   mexErrMsgTxt("Ten inputs required!");
  if(nlhs!=4) 	mexErrMsgTxt("Four outputs required!");
 

  // getting Gamma from the first input
  m = mxGetM(prhs[0]);
  n = mxGetN(prhs[0]);
  assert( m==2 && n==2);
  vector<vector<double> >  Gamma(m, vector<double>(n));
  data = mxGetPr(prhs[0]);
  for (int i=0; i<m; i++) {
    for (int j=0; j<n; j++) {
      Gamma[i][j] = data[j*m+i];  // column mode 
    }
  }

  // getting grid from the second input
  m = mxGetM(prhs[1]);
  n = mxGetN(prhs[1]);
  assert( m==2);
  vector<vector<double> >  grid(m, vector<double>(n));
  data = mxGetPr(prhs[1]);
  for (int i=0; i<m; i++) {
    for (int j=0; j<n; j++) {
      grid[i][j] = data[j*m+i];  // column mode 
    }
  }

  if( mxGetM(prhs[2])*mxGetN(prhs[2]) !=1 ) mexErrMsgTxt("The last input must be a scalar 4, 8, 16, 32 etc <=1024.");
  int degree = (int) mxGetScalar(prhs[2]);
  
  //getting alpha, lambda, nu1,nu2,P12, nserv
  vector<double> nu(2);
  double alpha = mxGetScalar(prhs[3]);  
  double lambda = mxGetScalar(prhs[4]);  
  nu[0] = mxGetScalar(prhs[5]); 
  nu[1] = mxGetScalar(prhs[6]);
  double P12 = mxGetScalar(prhs[7]);
  double nserv = mxGetScalar(prhs[8]);

  

  //getting parameters for the reference density
  m = mxGetM(prhs[9]);
  n = mxGetN(prhs[9]);
  assert(m==3 && n==2);
  vector<double> ref_a(n);   
  vector<double> ref_b0(n);
  vector<double> ref_b1(n); 
  data = mxGetPr(prhs[9]);
  
  for (int i=0; i<n; i++)
  {
  	ref_a[i] = data[i*m];
  	ref_b0[i] = data[i*m+1];
  	ref_b1[i] = data[i*m+2];
  }

  gaussianReference reference(ref_a, ref_b0, ref_b1);

  gridIndex gi(grid);
  int N = gi.getN();
  
  assert(nlhs==4);
  plhs[0] = mxCreateDoubleMatrix(N*36, 1, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(N*36, 1, mxREAL);
  plhs[2] = mxCreateDoubleMatrix(N*36, 1, mxREAL);
  plhs[3] = mxCreateDoubleMatrix(N, 1, mxREAL);
  double * Ar = mxGetPr(plhs[0]);
  double * Ac = mxGetPr(plhs[1]);
  double * Av = mxGetPr(plhs[2]);
  double * b = mxGetPr(plhs[3]);
  fE_run(Gamma, grid, reference, Ar, Ac, Av, b, degree, alpha, lambda, nu, P12, nserv);


}

