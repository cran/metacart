#' Permutation test
#'
#' Perform permutation test for an RE-tree
#' @param mf the data object of the RE-tree
#' @param nsplit the number of splits in the RE-tree
#' @param P the number of permuted data sets
#' @param sss boolean indicating whether the SSS strategy is used or not.
#' @param lookahead a boolean argument indicating whether to apply the "look-ahead" strategy when fitting the tree
#' @param minbucket the minimum number of observations in any terminal <leaf> node.
#' @param minsplit the minimum number of observations that must exist in a node in order for a split to be attempted.
#' @param cp complexity parameter as in rpart.
#' @param alpha.endcut parameter used in the splitting algorithm to avoid the endcut preference problem.
#' @param a parameter used in the sss to determine the slope of the logistic function that replaces the indicator function.
#' @param multi.start boolean indicating whether multiple starts must be used
#' @param n.starts number of multiple starts
#' @return a vector of Q-between computed from the permuted data sets
#' @keywords internal
permuteRE <- function(mf, nsplit, P = 999, sss, lookahead, 
                        minbucket = 3, minsplit = 6, cp = 0.0001,
                        alpha.endcut = 0.02, a = 50,
                        multi.start=T, n.starts=3){
  Qb.p <- numeric(P)

  nsplit.init <- nsplit
  if(nsplit == 1 & lookahead == TRUE) {
    nsplit <- 2
  }
  for (p in 1:P){
    inx.p <- sample(1:nrow(mf))
    mf.p <- mf
    #mf.p <- mf
    cols.mods <-  match(labels(terms(mf)), colnames(mf))
    mf.p[, cols.mods] <- mf[inx.p, cols.mods]
    if (sss == TRUE) {
      tree.p <- REmrt_SSS(mf.p, maxL = nsplit, minbucket = minbucket, 
                          minsplit = minsplit, cp = cp, lookahead = lookahead,
                          alpha.endcut = alpha.endcut, a = a,
                          multi.start = multi.start, n.starts = n.starts)
    } else {
      tree.p <- REmrt_GS_(mf.p, maxL = nsplit,  minbucket = minbucket, 
                          minsplit = minsplit, cp = cp, lookahead = lookahead)
    }

    if(nsplit.init > 1){
      Qb.p[p] <- tree.p$tree$Qb[nsplit+1]
    }
    else { 
      Qb.p[p] <- tree.p$tree$Qb[2]
    }
    
  }
  Qb.p
}