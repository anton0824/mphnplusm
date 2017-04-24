#ifndef _GaussianReference_H
#define _GaussianReference_H
#include <vector>
#include <cmath>
#include <assert.h>
using namespace std;

class gaussianReference {
 private:
  vector<double> a, b0, b1;
  double c0, c1;
   
 public:
  gaussianReference() {};
  gaussianReference(const vector<double>& a_, const vector<double> & b0_, 
                    const vector<double>& b1_) 
    {
      assert( a_.size() == 2 && b0_.size() == 2 && b1_.size() == 2);
      a = a_;
      b0 = b0_;
      c0 = a[1]*b0[1]*b0[1]-a[0]*b0[0]*b0[0];
      b1 = b1_;
      c1 = a[1]*b1[1]*b1[1]-b1[0];
    }
  
  double operator()(const vector<double>& y) {
    return marginal0(y[0])*marginal1(y[1]);
  }
  
  double marginal0(const double x) {
    return (x>0)? exp(-a[0]*(x-b0[0])*(x-b0[0])-c0):
      // minus c to make sure the 
      // density is continous.
      exp(- a[1]*(x-b0[1])*(x-b0[1]));
  }
  
  double marginal1(const double x) {
    return (x>0)? exp(-(x-b1[0])-c1):
      // minus c to make sure the 
      // density is continous.
      exp(- a[1]*(x-b1[1])*(x-b1[1]));
  }
};
#endif