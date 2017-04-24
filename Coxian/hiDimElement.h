#ifndef _HI_DIM_ELEMENT_
#define _HI_DIM_ELEMENT_
#include "oneDimElement.h"
class hiDimElement {
  int dim;
  vector<oneDimElement> marginals;
  
public:

    hiDimElement() {};
    
    hiDimElement(const vector<vector<double> > & grid)
    {
        setGrid(grid);
    }

    void setGrid(const vector<vector<double> > & grid) {
        dim = grid.size();
        marginals.resize(dim);
        for (int i=0; i< dim; i++)
            marginals[i].setGrid(grid[i]);
    }        

  double elem(const vector<int> & i_index,
	      const vector<int> & r_index,
	      const vector<double>& y)  { 
    double product = 1.0;
    for (int i=0; i< dim; i++) 
      product *=   marginals[i].elem(i_index[i], r_index[i], y[i]);
    return product;
  }


  double d_elem(const vector<int> & i_index,
		const vector<int> & r_index,
		const vector<double>& y,
		const int k)  { 
    assert( k>= 0 && k< dim);
    
    double product = 1.0;
    for (int i=0; i<k; i++) 
      product *=  marginals[i].elem(i_index[i], r_index[i], y[i]);
    product *=  marginals[k].d_elem(i_index[k], r_index[k], y[k]);
    for (int i =k+1; i<dim; i++) 
      product *=  marginals[i].elem(i_index[i], r_index[i], y[i]);
    return product;
  }
  

  double dd_elem(const vector<int> & i_index,
		 const vector<int> & r_index,
		 const vector<double>& y,
		 const int k_,
		 const int ell_)  { 
    int k, ell;
    if (k_> ell_) {
      k= ell_;
      ell = k_;
    } else
        {
            k=k_;
            ell = ell_;
        }
    // make sure k <= ell;
    assert( k>= 0 && k< dim);
    assert( ell>= 0 && ell< dim);

    
    double product = 1.0;


    for (int i=0; i<k; i++) 
      product *=  marginals[i].elem(i_index[i], r_index[i], y[i]);

    if (k==ell)
      product *=  marginals[k].dd_elem(i_index[k], r_index[k], y[k]);
    else 
      {
	product *=  marginals[k].d_elem(i_index[k], r_index[k], y[k]);
	for (int i= k+1; i<ell; i++)
	  product *=  marginals[i].elem(i_index[i], r_index[i], y[i]);
	product *=  marginals[ell].d_elem(i_index[ell], r_index[ell], y[ell]);
      }

    for (int i =ell+1; i<dim; i++) 
      product *=  marginals[i].elem(i_index[i], r_index[i], y[i]);
    return product;
  }

};

#endif
