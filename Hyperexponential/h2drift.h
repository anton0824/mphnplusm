#ifndef _H2DRIFT_
#define _H2DRIFT_
#include <vector>
#include <math.h>
#include "drift.h"
using namespace std;

class hyper2drift : drift {
 private:
  double alpha, beta, p;
  vector<double> nu;
  double mean;
  int nserv;

 public: 
  hyper2drift(){};
  hyper2drift(const double alpha_, const double beta_, const double p_,
	      const vector<double>& nu_, const int nserv_) {  
    alpha = alpha_;
    beta=beta_;
    p= p_;
	nserv = nserv_;
    nu= nu_;
    assert ( 0<= p && p<=1);
    assert (nu.size()==2);
    assert (nu[0]>0 && nu[1]>0);
    mean= p/nu[0] + (1-p)/nu[1];
    assert( mean == 1.0 );
  }
  
  double operator()(const int i, const vector<double> & y){
    assert (i==0 || i==1);
    return (i == 0)?  -beta*p/mean - nu[0]* y[0] +
      p*(nu[0]-alpha)*max(y[0]+y[1],0.0):
      -beta*(1-p)/mean - nu[1]* y[1] +
      (1-p)*(nu[1]-alpha)*max(y[0]+y[1],0.0);
  }
  double mu_tilde(const int i, const vector<double> & y) {

        assert (i==0 || i==1);
    return (i == 0)?  -(alpha-nu[0])*p*sqrt(double(nserv))/nu[0] - alpha*y[0]  +
	 (alpha-nu[0])*max((y[0] + p*sqrt(double(nserv))/nu[0] - p*max(y[0]+y[1],0.0)),0.0):
	 -(alpha-nu[1])*(1-p)*sqrt(double(nserv))/nu[1] - alpha*y[1]  +
	 (alpha-nu[1])*max((y[1] + (1-p)*sqrt(double(nserv))/nu[1] - (1-p)*max(y[0]+y[1],0.0)),0.0);
	 
    }
  double mu_v2(const int i, const vector<double> & y) {

        assert (i==0 || i==1);
    return (i == 0)?  - nu[0]* y[0] +
      p*(nu[0]-alpha)*max(y[0]+y[1],0.0):
      - nu[1]* y[1] +
      (1-p)*(nu[1]-alpha)*max(y[0]+y[1],0.0);
	 
    }

};
#endif
