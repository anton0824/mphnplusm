#ifndef _DRIFT_
#define _DRIFT_

#include <vector>
class drift{
 private:

 public: 

 virtual double operator()(const int i, const vector<double> & y) = 0;

};

#endif
