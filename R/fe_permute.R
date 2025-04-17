#' A function returns the Q-between from the tree with given size
#' @param tree A initial tree fitted by rpart, needs to an rpart object.
#' @param nsplit the required number of splits
#' @param mods the moderators found by the tree.
#' @param minbucket the minimum number of observations in any terminal <leaf> node.
#' @param minsplit the minimum number of observations that must exist in a node in order for a split to be attempted.
#' @param ... Additional arguments passed to prune.rpart().
#' @return The pruned tree
#' @importFrom rpart prune
#' @importFrom methods is
#' @keywords internal
Q_selected_size_GS <- function(tree, nsplit, mods, minbucket, minsplit,...){
  # 
  # if the tree with required size is not directly available
  # first prune the tree to the smaller size
  # then further split the tree until reaching the required size
  tsplit <- max(tree$cptable[tree$cptable[,2] <= nsplit,2])
  prunedTree <- prune(tree,tree$cptable[tree$cptable[,2] == tsplit,1])
  nodes <- as.numeric(rownames(prunedTree$frame)[prunedTree$where])
  Qb.pruned <-prunedTree$frame$dev[1] - 
    sum(prunedTree$frame$dev[prunedTree$frame$var == "<leaf>"])
  if (tsplit < nsplit) {
    Qincrease <- NULL
    for (i in 1:(nsplit - tsplit)) {
      Dev <- 0
      candidate.pleaf <- names(table(nodes))[table(nodes) >= minsplit]
      if(length(candidate.pleaf) == 0) break
      for (j in as.numeric(candidate.pleaf)) {
        inx.node <- nodes == j
        yi <- tree$y[inx.node]
        vi <- 1/tree$wt[inx.node]
        g <- sum(yi/vi)/sum(1/vi)
        Q <- sum((yi-g)^2/vi)
        for (k in 1:ncol(mods)){
          new.nodes <- nodes
          xk <- mods[inx.node, k]
          if (is.numeric(xk)) {
            xk.ord <- xk
          } else {
            xk.rank <- rank(tapply(yi, xk, mean))
            xk.ord <- xk.rank[as.character(xk)]
          }
          temp <- fe_cutoff_GS(yi, vi, xk.ord, minbucket = minbucket)
          if (is.null(temp)) {
            Dev.new <- 0
          } else {
            Dev.new <-  Q - sum(temp[2:3])
          }
          if (Dev.new > Dev) {
            Dev <- Dev.new
            pleaf <- j
            chosenx <- xk.ord
            xstar <- temp[1]
          }
        }
      }
      nodes[nodes == pleaf] <- ifelse(chosenx < xstar, pleaf*2, pleaf*2+1)
      Qincrease <- c(Qincrease, Dev)
    }
    #print(c(tsplit, nsplit)) to check
    #print(Qincrease) to check
    Qb <- Qb.pruned + sum(Qincrease)
  } else {
    Qb <- Qb.pruned
  }
  return(Qb)
}

Q_selected_size_sss <- function(initial.tree, cptable, 
                                nsplit, minsplit,
                                y, vi, mods,
                                a, alpha.endcut,
                                multi.start, n.starts){
  # 
  # if the tree with required size is not directly available
  # first prune the tree to the smaller size
  # then further split the tree until reaching the required size
  tsplit <- max(cptable[cptable[,2] <= nsplit,2])
  prunedtree <- fe_prune(initial.tree, cptable[cptable[,2] == tsplit,1])
  if (tsplit == 0) {
    nodes <- rep(1, nrow(mods)) 
  } else {
    nodes <- numeric(nrow(mods))
    for (k in as.character(prunedtree[[4]][,1])){
      nodes[prunedtree[[2]][[k]]] <- as.numeric(k)
      }
    }
 
  if(nrow(prunedtree[[1]]) > 0) {
    Qb.pruned <- prunedtree[[1]]$Q[1] - sum(prunedtree[[4]][,3])
  } else{
    Qb.pruned <- 0
  }
  if (tsplit < nsplit) {
    Qincrease <- NULL
    for (i in 1:(nsplit - tsplit)) {
      Dev <- 0
      candidate.pleaf <- names(table(nodes))[table(nodes) >= minsplit]
      if(length(candidate.pleaf) == 0) break
      for (j in as.numeric(candidate.pleaf)) {
        inx.node <- nodes == j
        g <- sum(y[inx.node]/vi[inx.node])/sum(1/vi[inx.node])
        Q <- sum((y[inx.node]-g)^2/vi[inx.node])
        for (k in 1:ncol(mods)){
          new.nodes <- nodes
          xk <- mods[inx.node, k]
          if (is.numeric(xk)) {
            xk.ord <- xk
          } else {
            xk.rank <- rank(tapply(y, xk, mean))
            xk.ord <- xk.rank[as.character(xk)]
          }
          temp <- fe_cutoff_SSS(y[inx.node], vi[inx.node], xk.ord,
                                a = a, alpha.endcut = alpha.endcut,
                                multi.start = multi.start, n.starts = n.starts)
          if (is.null(temp)) {
            Dev.new <- 0
          } else {
            Dev.new <-  Q - sum(temp[2:3])
          }
          if (Dev.new > Dev) {
            Dev <- Dev.new
            pleaf <- j
            chosenx <- xk.ord
            xstar <- temp[1]
          }
        }
      }
      nodes[nodes == pleaf] <- ifelse(chosenx < xstar, pleaf*2, pleaf*2+1)
      Qincrease <- c(Qincrease, Dev)
    }
    #print(c(tsplit, nsplit)) to check
    #print(Qincrease) to check
    Qb <- Qb.pruned + sum(Qincrease)
  } else {
    Qb <- Qb.pruned
  }
  return(Qb)
}

#' A function returns the Q-between from the tree with given size
#' @param mf data transformed to fit the tree.
#' @param Call The matched call.
#' @param nsplit the required number of splits.
#' @param P the number of permuted data sets.
#' @param sss boolean indicating whether the SSS strategy is used or not.
#' @param lookahead a boolean argument indicating whether to apply the "look-ahead" strategy when fitting the tree
#' @param minbucket the minimum number of observations in any terminal <leaf> node.
#' @param minsplit the minimum number of observations that must exist in a node in order for a split to be attempted.
#' @param cp complexity parameter as in rpart.
#' @param maxdepth set the maximum depth of any node of the final tree, with the root node counted as depth 0.
#' @param alpha.endcut parameter used in the splitting algorithm to avoid the endcut preference problem.
#' @param a parameter used in the sss to determine the slope of the logistic function that replaces the indicator function.
#' @param multi.start boolean indicating whether multiple starts must be used
#' @param n.starts number of multiple starts
#' @param ... Additional arguments passed to prune.rpart().
#' @return The pruned tree
#' @importFrom rpart rpart
#' @keywords internal
permuteFE <- function(mf, Call, nsplit, P = 100, sss, lookahead, 
                       minbucket = 3, minsplit = 6, cp = 0.0001,
                       maxdepth = 10, alpha.endcut = 0.02, a = 50,
                       multi.start=T, n.starts=3,...){
  # a function that runs permutation test for FEmrt
  # and returns the computed Q-between for each permuted data set
  Qb.p <- numeric(P)
  if (sss == FALSE) {
    for (p in 1:P){
      inx.p <- sample(1:nrow(mf))
      mf.p <- mf
      cols.mods <-  match(labels(terms(mf)), colnames(mf))
      mf.p[, cols.mods] <- mf[inx.p, cols.mods]
      mf.p$wts <- 1/mf.p$`(vi)`
      tree.p <- rpart(as.formula(Call), weights = mf.p$wts, data = mf.p,
                      control = rpart.control(maxdepth = maxdepth, minbucket = minbucket, minsplit = minsplit, cp = cp),
                      x = TRUE)
      mods <- model.frame(delete.response(tree.p$terms), mf.p)
      Qb.p[p] <- Q_selected_size_GS(tree.p, nsplit, mods = mods,
                                    minbucket = minbucket, minsplit = minsplit)
      
    }
  } else{
    for (p in 1:P){
      inx.p <- sample(1:nrow(mf))
      mf.p <- mf
      cols.mods <-  match(labels(terms(mf)), colnames(mf))
      mf.p[, cols.mods] <- mf[inx.p, cols.mods]
      initial.tree <- FEpartitionSSS.wrapper(mf.p, minsplit, maxdepth, cp,
                                             a, alpha.endcut, multi.start, n.starts)
      cps <-  sort(unique(initial.tree[[1]]$cp))
      split.left <- sapply(fe_prune_path(initial.tree[[1]], cps), length)
      nsplits <- cumsum(split.left[length(split.left):1])
      CP <- c(cps[length(cps):1], cp)
      cptable <- cbind(CP, nsplits)
      mods <- mf.p[, cols.mods]
      y <- model.response(mf.p)
      vi <- c(t(mf.p["(vi)"]))
      Qb.p[p] <- Q_selected_size_sss(initial.tree = initial.tree, cptable = cptable, 
                                     nsplit = nsplit, minsplit = minsplit,
                                     y = y, vi = vi, mods = mods,
                                     a = a, alpha.endcut = alpha.endcut,
                                     multi.start = multi.start, 
                                     n.starts = n.starts)
      
    }
    
  }
  
  Qb.p
}


