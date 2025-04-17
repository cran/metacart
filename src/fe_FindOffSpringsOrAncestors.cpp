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
IntegerVector find_offsprings_(int nodeID, IntegerVector allNodes) {
  IntegerVector nodeIDv;
  nodeIDv.push_back(nodeID);
  if (all(is_na(match(nodeIDv, allNodes)))) return R_NilValue;
  IntegerVector left = find_offsprings_(nodeID*2, allNodes);
  IntegerVector right = find_offsprings_(nodeID*2+1, allNodes);
  IntegerVector res;
  res.push_back(nodeID);
  for (int i = 0; i < left.length(); i++) {
    res.push_back(left[i]);
  }
  for (int j = 0; j < right.length(); j++) {
    res.push_back(right[j]);
  }
  return res;
}

// [[Rcpp::export]]
IntegerVector find_children_vec(IntegerVector nodeIDv, IntegerVector allNodes){
  IntegerVector res = clone(nodeIDv);
  for (int i = 0; i < nodeIDv.length(); i++) {
    IntegerVector temp;
    temp.push_back(nodeIDv[i]*2);
    temp.push_back(nodeIDv[i]*2 + 1);
    res = union_(res, temp);

  }
  //res = setdiff(res, nodeIDv); // remove themselves
  return res;
}

// [[Rcpp::export]]
IntegerVector find_ancestor_(int nodeID) {
  IntegerVector res;
  // res.push_back(nodeID);
  while (true) {
    nodeID = nodeID/2;
    if (nodeID <= 0) break;
    res.push_back(nodeID);
  }
  return res;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

