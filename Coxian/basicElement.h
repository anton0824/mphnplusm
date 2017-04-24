#ifndef _BASIC_ELEMENT_
#define _BASIC_ELEMENT_
#include <assert.h>

class basicElement {

public:
  
  static double phi(const double x) { 
    assert( x>=-1 && x<=1);
    return ( x> 0) ? (x-1)*(x-1)*(2*x+1): (-x-1)*(-x-1)*(-2*x+1);
  }
  static double psi(const double x) {
    assert( x>=-1 && x<=1);
    return ( x>0)?  x*(x-1)*(x-1) : x*(-x-1)*(-x-1);
  }

  static double d_phi(const double x) {
    assert( x>=-1 && x<=1);
    return (x>0)?  6*x*(x-1): -6*x*(x+1);
  }

  static double d_psi(const double x) {
    assert( x>=-1 && x<=1);
    return (x>0)?  (x-1)*(3*x-1): (x+1)*(3*x+1);
  }

  static double dd_phi(const double x) {
    assert( x>=-1 && x<=1);
    return (x>0)?  12*x-6: -12*x-6;
  }

  static double dd_psi(const double x) {
    assert( x>=-1 && x<=1);
    return (x>=0)?  6*x-4: 6*x+4;
  }

};

#endif
