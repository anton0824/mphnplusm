#include "mex.h"
#include "density.h"
//#include "h2drift.h"
#include <vector>
#include <assert.h>
using namespace std;


void mexFunction( int nlhs, mxArray *plhs[],

				 int nrhs, const mxArray *prhs[])

{

  int m, n;
  double *data;

  if(nrhs!=11)   mexErrMsgTxt("Eleven inputs required.");
  if(nlhs!=1) 	mexErrMsgTxt("One output required.");
 

  // getting Gamma from the first input
  m = mxGetM(prhs[0]);
  n = mxGetN(prhs[0]);
  assert( m==2 && n==2);
  vector<vector<double> >  Gamma(m, vector<double>(n));
  data = mxGetPr(prhs[0]);
  for (int i=0; i<m; i++) {
    for (int j=0; j<n; j++) {
      Gamma[i][j] = data[j*m+i];  // column mode 
      // mexPrintf(" Gamma[%d][%d] = %f ", i, j, Gamma[i][j]);
    }
    // mexPrintf("\n");
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

  gridIndex gi(grid);
  int N = gi.getN();

  // getting solution u from the third input
  m = mxGetM(prhs[2]);
  n = mxGetN(prhs[2]);
  assert( m== 1);
  assert( n == N);

  double * u= mxGetPr(prhs[2]);
  

  if( mxGetM(prhs[3])*mxGetN(prhs[3]) !=1 ) 
    mexErrMsgTxt("The number of meshes within a grid >=1 "); 
  int  numberWithinGrid = (int) mxGetScalar(prhs[3]);
  // mexPrintf(" The number of meshes within a grid = %d \n", numberWithinGrid);

  density d(grid, numberWithinGrid);
  int totalMeshSize = d.getTotalMeshSize();


  plhs[0] = mxCreateDoubleMatrix(totalMeshSize, 3, mxREAL);

  double * p = mxGetPr(plhs[0]);

  //getting alpha, lambda, srate, nserv
  vector<double> nu(2);
  double alpha = mxGetScalar(prhs[4]); 
  double lambda = mxGetScalar(prhs[5]);  
  nu[0] = mxGetScalar(prhs[6]); 
  nu[1] = mxGetScalar(prhs[7]);
  double P12 = mxGetScalar(prhs[8]); 
  double nserv = mxGetScalar(prhs[9]);
  
     
  //vector<double> nu(2);
  //nu[0] = 1/m1;
  //nu[1] = (1-pr)/(1-pr/nu[0]);
 
  //hyper2drift mu(alpha, beta, pr, nu, nserv);

  generator g(Gamma, grid, alpha, lambda, nu[0], nu[1], P12, nserv);
  
  //getting parameters for the reference density
  m = mxGetM(prhs[10]);
  n = mxGetN(prhs[10]);
  assert(m==3 && n==2);
  vector<double> ref_a(n);   
  vector<double> ref_b0(n);
  vector<double> ref_b1(n);
  data = mxGetPr(prhs[10]);
  
  for (int i=0; i<n; i++)
  {
  	ref_a[i] = data[i*m];
  	ref_b0[i] = data[i*m+1];
  	ref_b1[i] = data[i*m+2];
  }

 gaussianReference reference(ref_a, ref_b0, ref_b1);

 d.compute_density(g, reference, u, p);  

// for (int i=0; i< totalMeshSize; i++) 
//    mexPrintf("p[%d] = (%f, %f, %f)\n", i,
//	      p[i],
//	      p[totalMeshSize+i],
//	      p[2*totalMeshSize+i]);
}

