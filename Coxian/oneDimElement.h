#ifndef _ONE_DIM_ELEMENT_
#define _ONE_DIM_ELEMENT_
#include <vector>
#include "basicElement.h"
using namespace std;

class oneDimElement {
private:
    vector<double> x, h;
  
public:

    oneDimElement() {};

    oneDimElement(const vector<double>& x_) {
      setGrid(x_);
    }


  void setGrid(const vector<double>& x_) { 
    x = x_; 
    int gridSize=x.size();
    h.resize(gridSize);
    for (int  i=0; i< gridSize-1; i++)
      h[i]= x[i+1]- x[i];
    h[gridSize-1]=1.0; // this line should not be here;
    // It is a lazy way to take care of the end of the 
    // grid point; otherwise, line x below has a problem when y is the 
    // right most end point of the grid 
  }
  
  double elem(const int k, const int r, const double y) const {
      int gridSize=x.size();
      assert( k>=0 && k< gridSize);
      assert( r==0 || r==1);
      assert( y>= x[0] && y<= x[gridSize-1]);
      double temp= y-x[k];
      if ( r == 0)
	return (temp <0) ?  basicElement::phi(temp/h[k-1])
	  //line x
	  : basicElement::phi(temp/h[k]);
      else 
	return (temp <0) ?  h[k-1]*basicElement::psi(temp/h[k-1])
	  : h[k]*basicElement::psi(temp/h[k]); 
  }
  

  double d_elem(const unsigned k, const int r, const double y) const {
    assert( k>=0 && k< x.size());
    double temp= y-x[k];
    if ( r == 0)
      return (temp <0) ?  basicElement::d_phi(temp/h[k-1])/h[k-1]
	//line x
	: basicElement::d_phi(temp/h[k])/h[k];
    else 
      return (temp <0) ?  basicElement::d_psi(temp/h[k-1])
	: basicElement::d_psi(temp/h[k]); 
  }
  


  double dd_elem(const unsigned k, const int r, const double y) const {
    assert( k>=0 && k< x.size());
    double temp= y-x[k];
    if ( r == 0)
      return (temp <0) ?  basicElement::dd_phi(temp/h[k-1])/h[k-1]/h[k-1]
	//line x
	: basicElement::dd_phi(temp/h[k])/h[k]/h[k];
    else 
      return (temp <0) ?  basicElement::dd_psi(temp/h[k-1])/h[k-1]
	: basicElement::dd_psi(temp/h[k])/h[k]; 
  }
  

  
/*   // for testing purpose  */
/*   void fixElement(const int k_) { k = k_;} */

/*   double operator()(const double y) const { */
/*     return psi(k, y); */
/*   } */

};

#endif
