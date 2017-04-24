#ifndef _INDEX_
#define _INDEX_
#include <vector>
#include <assert.h>

using namespace std;

class gridIndex{

 private:
  vector<vector<double> > grid;
  int dim;
  int N;
  vector<int> gridLowLimits;
  vector<int> gridUpperLimits;
  

 public:
  gridIndex(){};
  gridIndex( const vector<vector<double> > & grid_) {
    grid=grid_;
    dim = grid.size();
    gridLowLimits.resize(dim);
    gridUpperLimits.resize(dim);
    N = 1;
    for (int i=0; i<dim; i++) {
      gridLowLimits[i] = 1;
      gridUpperLimits[i] = grid[i].size()-2;
      assert(gridLowLimits[i] <= gridUpperLimits[i]);
      N *= grid[i].size()-2;  // number of interior grid points;
    }
    N *=4;
  }


  int getN() { return N;}
  int getIndex(const vector<int> i_index, const vector<int> r_index) {

    for (int i=0; i<dim; i++) {
      assert( i_index[i] >= gridLowLimits[i] && i_index[i]<= gridUpperLimits[i]);
      assert( r_index[i] == 0 || r_index[i] ==1); 
    }
      
    return
      4*( (i_index[1]-1)*gridUpperLimits[0] + i_index[0]-1) + 2* r_index[1]+ r_index[0];

  }

};


#endif
