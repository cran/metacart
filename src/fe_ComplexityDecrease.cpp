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
NumericVector complexity_decrease_(int inNode, NumericVector delQ,
                            IntegerVector pnode) {
  IntegerVector inNodeV;
  inNodeV.push_back(inNode);
  IntegerVector temp = match(inNodeV, pnode);
  // base case
  if (any(is_na(temp))) {
    NumericVector resBase;
    resBase.push_back(0);
    resBase.push_back(0);
    return resBase;
  }
  else {
    int inx = temp[0];
    double res1 = delQ[inx - 1];
    NumericVector childLeft = complexity_decrease_(inNode*2, delQ, pnode);
    NumericVector childRight = complexity_decrease_(inNode*2 + 1, delQ, pnode);
    double temp1 = res1 + childLeft[0] + childRight[0];
    double temp2 = 1 + childLeft[1] + childRight[1];
    NumericVector res;
    res.push_back(temp1);
    res.push_back(temp2);
    return res;
  }
  
}
// double complexity_decrease_(int inNode, NumericVector delQ,
//                                    IntegerVector pnode) {
//   IntegerVector inNodeV;
//   inNodeV.push_back(inNode);
//   IntegerVector temp = match(inNodeV, pnode);
//   // base case
//   if (any(is_na(temp))) {
//     return 0;
//   }
//   else {
//     int inx = temp[0];
//     double res = delQ[inx-1];
//     double childLeft = complexity_decrease_(inNode*2, delQ, pnode);
//     double childRight = complexity_decrease_(inNode*2 + 1, delQ, pnode);
//     return res + childLeft + childRight;
//     
//   }
//   
// }

// int test_match(IntegerVector inNode, IntegerVector pnode) {
//   int inx = match(inNode, pnode)[0];
//   if (any(is_na(match(inNode, pnode)))) return 0;
//   else return inx;
// }

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

