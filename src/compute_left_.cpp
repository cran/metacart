#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
NumericVector compute_left_(NumericVector x1,NumericVector x2, 
                           NumericVector x3, NumericVector xuni){
  // x1 is the effect size in the unsplit leaves
  // x2 is the sampling variance in the unsplit leaves
  // x3 is the leaves memberships
  // xuni is the unique leaves 
  int i;
  double sumQoleaf = 0;
  double sumColeaf = 0;
  for (i = 0; i < xuni.length(); i++) {
    double tempWY2 = 0;
    double tempWY = 0;
    double tempW = 0;
    double tempW2 = 0;
    int ii;
    for (ii = 0; ii < x3.length(); ii++) {
      if (x3[ii] == xuni[i]) {
        tempWY2 = tempWY2 + pow(x1[ii], 2)/x2[ii];
        tempWY = tempWY + x1[ii]/x2[ii];
        tempW = tempW + 1/x2[ii];
        tempW2 = tempW2 + pow(x2[ii], -2);
      }
    }
    sumQoleaf = sumQoleaf + tempWY2 - pow(tempWY,2)/tempW;
    sumColeaf = sumColeaf + tempW2/tempW;
  }
  NumericVector res;
  res.push_back(sumQoleaf);
  res.push_back(sumColeaf);
  return res;
}




