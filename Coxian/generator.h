#ifndef _GENERATOR_
#define _GENERATOR_
#include <vector>
#include <assert.h>
#include "hiDimElement.h"
//#include "h2drift.h"
#include <math.h> 
using namespace std;
typedef double (* functionType)(vector<double> );

class generator {
private:
    int dim;
	double alpha;
	double lambda;
	double rho;
	double nu1;
	double nu2;
	double P12;
	double nserv;
	double m;
	double util1;
	double util2;
    vector<vector<double> > Gamma;
    //hyper2drift mu;
    hiDimElement e;

public:

    generator() {};

    
    generator( const vector<vector<double> > & Gamma_,
               //const hyper2drift& mu_,
               const vector<vector<double> > & grid,
			   const double alpha_,
			   const double lambda_,
			   const double nu1_,
			   const double nu2_,
			   const double P12_,
			   const double nserv_) {

        Gamma = Gamma_;
		alpha = alpha_;
		lambda = lambda_;
		nu1 = nu1_;
		nu2 = nu2_;
		P12 = P12_;
		m = 1/(1/nu1 + P12/nu2);
		util1 = (1/nu1)/m;
		util2 = (P12/nu2)/m;
		nserv = nserv_;
		rho = lambda*m/nserv;
        e.setGrid(grid);
        assert( Gamma.size() == grid.size());
        dim = Gamma.size();
    }

	// This function evalues the generator of the diffusion as in (2.6) of Dai-He 2013
    double Lf(const vector<int>& i_index,
              const vector<int>& r_index,
              const vector<double>& y) {
				  
        double sum =0.0;
		double del = 1/sqrt(nserv);
		double aqinf = max(lambda -nserv/m, 0.0); //alpha * equilibrium queue length
		double rhominone = min(rho, 1.0); //The smaller of rho and one
		// The terms corresponding to the drift of the diffusion process
		sum += (del*(lambda - nu1*util1*nserv) - alpha*max(y[0]+y[1], 0.0) - nu1*(y[0] - max(y[0]+y[1], 0.0)))* e.d_elem(i_index, r_index, y, 0); 
		sum += (del*(P12*nu1*util1*nserv -nu2*util2*nserv) -nu2*y[1] + P12*nu1*(y[0] - max(y[0]+y[1], 0.0)))* e.d_elem(i_index, r_index, y, 1);
		
				
		// One may choose to use either a constant, or state-dependent diffusion coefficient for the diffusion approximation
		// Uncomment each block of code as necessary (only one of the blocks should be uncommented at a time)
		// Don't forget to run the mex command in MATLAB to recompile the c++ code after making modifications.
				
		// constant diffusion coefficient
        // sum += .5* (del*del)*(lambda+ rhominone*nu1*util1*nserv + aqinf) * e.dd_elem(i_index, r_index, y, 0, 0);//a_11
        // sum += .5* (-(del*del)*P12*rhominone*nserv*nu1*util1)* e.dd_elem(i_index, r_index, y, 0, 1);//a_12
        // sum += .5* (-(del*del)*P12*rhominone*nserv*nu1*util1)* e.dd_elem(i_index, r_index, y, 1, 0);//a_21
        // sum += .5* (del*del)*rhominone*nserv*(nu2*util2+P12*nu1*util1)* e.dd_elem(i_index, r_index, y, 1, 1);//a_22 
        
		// state dependent diffusion coefficient
        sum += .5* ((del*del)*(lambda+ nu1*util1*nserv + aqinf) + (del*alpha*max(y[0]+y[1], 0.0) - del*del*aqinf + del*nu1*(y[0] - max(y[0]+y[1], 0.0))))* e.dd_elem(i_index, r_index, y, 0, 0);//a_11
        sum += .5* (-(del*del)*P12*nserv*nu1*util1 - (del*P12*nu1*(y[0] - max(y[0]+y[1], 0.0))))* e.dd_elem(i_index, r_index, y, 0, 1);//a_12
        sum += .5* (-(del*del)*P12*nserv*nu1*util1 - (del*P12*nu1*(y[0] - max(y[0]+y[1], 0.0))))* e.dd_elem(i_index, r_index, y, 1, 0);//a_21
        sum += .5* ((del*del)*nserv*(nu2*util2+P12*nu1*util1) + ((del*del)*nserv*(nu2-nu2*util2)+ del*nu2*min(y[1]+del*nserv*(util2-1),0.0) + del*P12*nu1*(y[0] - max(y[0]+y[1], 0.0))))* e.dd_elem(i_index, r_index, y, 1, 1);//a_22
        
		return sum;
    }
    
};

#endif


