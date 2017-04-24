#ifndef _GENERATOR_
#define _GENERATOR_
#include <vector>
#include <assert.h>
#include "hiDimElement.h"
#include "h2drift.h"
#include <math.h> 
using namespace std;
typedef double (* functionType)(vector<double> );

class generator {
private:
    int dim;
	int nserv;
    vector<vector<double> > Gamma;
    hyper2drift mu;
    hiDimElement e;

public:

    generator() {};

    
    generator( const vector<vector<double> > & Gamma_,
               const hyper2drift& mu_,
               const vector<vector<double> > & grid,
			   const int nserv_) {

        Gamma = Gamma_;
        mu = mu_;
		nserv = nserv_;
        e.setGrid(grid);
        assert( Gamma.size() == grid.size());
        dim = Gamma.size();
    }

    double Lf(const vector<int>& i_index,
              const vector<int>& r_index,
              const vector<double>& y) {

        double sum =0.0;
		double del = 1/sqrt(nserv);
        for (int i=0; i<dim; i++){
			sum += mu(i,y)* e.d_elem(i_index, r_index, y, i);
		}
		// state dependent diffusion coefficient
        // for (int i=0; i<dim; i++){
                // sum += .5*(Gamma[i][i] - mu.mu_v2(i,y)/sqrt(double(nserv)))* e.dd_elem(i_index, r_index, y, i, i);
		// }
		// constant diffusion coefficient
        for (int i=0; i<dim; i++)
            for (int j=0; j<dim; j++)
                sum += .5* Gamma[i][j]* e.dd_elem(i_index, r_index, y, i, j); 
        return sum;
    }
    
};

#endif


