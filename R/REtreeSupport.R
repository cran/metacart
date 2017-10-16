#' A function to compute Q-between and residual heterogeneity
#'
#' @param yi the effect sizes
#' @param vi the sampling variances
#' @param mods the subgrouping moderator
#' @return Q-between and the residual heterogeneity under mixed effects model assumption
#' @keywords internal
rebetQ<- function(yi, vi, mods){
  wts = 1/vi
  wy = wts*yi
  wy2 = wts * yi^2
  Q <- tapply(wy2, mods, sum) - tapply(wy, mods, function(x) (sum(x))^2)/tapply(wts, mods, sum)
  df <- tapply(wy, mods, length)-1
  C <- tapply(wts, mods, sum) - tapply(wts, mods, function(x) sum(x^2))/ tapply(wts, mods, sum)
  tau2 <- (sum(Q) - sum(df))/sum(C)
  tau2 <- max(0, tau2)
  wstar = 1/(vi+tau2)
  wystar = wstar*yi
  wy2star = wstar*yi^2
  Qstar <- tapply(wy2star, mods, sum) - tapply(wystar, mods, function(x) (sum(x))^2)/tapply(wstar, mods, sum)
  Qstar.total <- sum(wy2star) - (sum(wystar))^2/sum(wstar)
  Qbet <- Qstar.total - sum(Qstar)
  if (is.na(Qbet)) {
    Qbet <- Inf
  }
  return(c(Qbet, tau2))

}


#' A function to subgroup newdata according to the fitted tree
#'
#' @param x the fitted RE meta-tree results
#' @param newdata the new data set
#' @return a matrix consists of node lables
#' @importFrom stats model.frame terms
#' @keywords internal
REmrt.prednode <- function(x, newdata){
  tt <- terms(x$data)
  ms <- model.frame(delete.response(tt), newdata)
  oms <- model.frame(delete.response(tt), x$data)
  tree <- x[["tree"]]
  # if (any(sapply(ms, class) != sapply(oms, class)))
  #   stop("The type of the variables do not match")
  if(nrow(tree) < 2) {
    pred.node <- rep(1, nrow(ms))
  } else {
    tnode <- rep(1, nrow(ms))
    nodes <- tnode
    for (i in 1:(nrow(tree) - 1)){
      tinx <- which(tnode == tree[i+1, "pleaf"])
      tempm <- ms[tree[i+1, "mod"]]
      if(sapply(tempm, is.numeric) == TRUE) {
        tnode[tinx] <- ifelse(tempm[tinx,1] <= x[["cpt"]][[i]], 2*i, 2*i+1)
      } else {
        tnode[tinx] <- ifelse(tempm[tinx,1] %in% oms[,tree[i+1, "mod"]],
                              ifelse(tempm[tinx,1] %in% x[["cpt"]][[i]], 2*i, 2*i+1),
                              NA)
      }
      nodes <- cbind(nodes, tnode)
    }
    pred.node <- nodes
  }
  pred.node
}


#' A function to grow the tree
#'
#' @param mf the data.frame to grow the tree
#' @param maxL the maximum number of splits
#' @param minisplit the minimal number of studies in a parent node to be split
#' @param delQ the stopping rule for decrease of between-subgroups Q. Any split that does not decrease the between-subgroups Q is not attempted.
#' @return a list including a tree, the split points, the data, and the nodes after each split
#' @importFrom stats terms model.response
#' @keywords internal
REmrt.fit1<- function(mf, maxL, minsplit = 2, delQ = 0.001){
  y <- model.response(mf)
  vi <- c(t(mf["(vi)"]))
  mods.names <-  labels(terms(mf))
  mods <- mf[mods.names]
  nmod <- ncol(mods)
  cpt <- list()
  nodemark <- data.frame(node = rep(1, nrow(mf)))
  res <- data.frame(Qb = rebetQ(y, vi, mods = nodemark)[1],
                    tau2 = rebetQ(y, vi, mods = nodemark)[2],
                    split = NA, mod = NA, pleaf = NA)
  delta.Q <- Inf
  for (i in 1:maxL){
    Dev<- -Inf
    cnode <- nodemark[ ,i]
    # Check if the parent leaves are all smaller than the minsplit
    # and only split nodes with subjects more than the minsplit
    len.node <- tapply(vi, cnode, length)
    nodes <- names(len.node) [len.node >= minsplit]
    if (length(nodes) == 0) {
      break
    } else {
      for (j in 1:length(nodes)){
        pleaf.inx <- cnode == as.numeric(nodes[j])
        for (k in 1:nmod){
          if(sapply(mods[mods.names[k]], is.numeric) == FALSE) {
            chosenmod <- mods[pleaf.inx, k]
            mod.order <- rank(tapply(y[pleaf.inx],chosenmod,mean))
            cmod.ordinal <- mod.order[as.character(chosenmod)]
            cpoints <- sort(mod.order)
            if (length(cpoints) >= 2) {
              for (g in 1:(length(cpoints)-1)) {
                cnode.test <- cnode
                cnode.test[pleaf.inx] <- ifelse( cmod.ordinal <= cpoints[g], 2*i, 2*i+1)
                temp <- rebetQ(y, vi, mods = as.factor(cnode.test))
                if (temp[1] > Dev) {
                  Dev <- temp[1]
                  msplit <- paste(mods.names[k], "=", paste(names(mod.order[mod.order <= cpoints[g]]), collapse = "/"), collapse = " ")
                  tres <- data.frame(Qb = temp[1], tau2 = temp[2],
                                     split = msplit, mod = mods.names[k], pleaf = as.numeric(nodes[j]))
                  tcpt <- names(mod.order[mod.order <= cpoints[g]])
                  tnode <- cnode.test
                }
              }
            }
          } else{
            chosenmod <- mods[pleaf.inx, k]
            cpoints <- sort(unique(chosenmod))
            if (length(cpoints) >= 2) {
              for (g in 1:(length(cpoints)-1)) {
                cnode.test <- cnode
                cnode.test[pleaf.inx] <- ifelse( chosenmod <= cpoints[g], 2*i, 2*i+1)
                temp <- rebetQ(y, vi, mods = as.factor(cnode.test))
                if (temp[1] > Dev) {
                  Dev <- temp[1]
                  tcpt <- cpoints[g]
                  msplit <- paste(mods.names[k], "<=", tcpt, collapse = " ")
                  tres <- data.frame(Qb = temp[1], tau2 = temp[2], split = msplit, mod = mods.names[k], pleaf = as.numeric(nodes[j]))
                  tnode <- cnode.test
                }
              }
            }
          }

        }
      }

      delta.Q <- abs(tres$Qb[1] - res$Qb[i])
      if (delta.Q < delQ) break
      new.node <- tnode
      nodemark <- cbind(nodemark, new.node)
      res <- rbind(res, tres)
      cpt[[i]] <- tcpt

    }



  }
  list(tree = res, node.split = nodemark, cpt = cpt, data = mf)

}

#' A function to fit the tree for cross-validation
#'
#' @param mf the data.frame to grow the tree
#' @param maxL the maximum number of splits
#' @return a list including a tree, the split points, the data, and the nodes after each split
#' @keywords internal
#' @importFrom stats terms
REmrt.fit0<- function(mf, maxL){
  y <- model.response(mf)
  vi <- c(t(mf["(vi)"]))
  mods.names <-  labels(terms(mf))
  mods <- mf[mods.names]
  nmod <- ncol(mods)
  cpt <- list()
  nodemark <- data.frame(node = rep(1, nrow(mf)))
  res <- data.frame(Qb = rebetQ(y, vi, mods = nodemark)[1],
                    tau2 = rebetQ(y, vi, mods = nodemark)[2],
                    split = NA, mod = NA, pleaf = NA)
  for (i in 1:maxL){
    Dev<- -Inf
    cnode <- nodemark[ ,i]
    len.node <- tapply(vi, cnode, length)
    nodes <- names(len.node) [len.node >= 2]
    if (length(nodes) == 0) break
    tres <- NULL

    for (j in 1:length(nodes)){
      pleaf.inx <- cnode == as.numeric(nodes[j])
      for (k in 1:nmod){
        if(sapply(mods[mods.names[k]], is.numeric) == FALSE) {
          chosenmod <- mods[pleaf.inx, k]
          mod.order <- rank(tapply(y[pleaf.inx],chosenmod,mean))
          cmod.ordinal <- mod.order[as.character(chosenmod)]
          cpoints <- sort(mod.order)
          if (length(cpoints) >= 2) {
            for (g in 1:(length(cpoints)-1)) {
              cnode.test <- cnode
              cnode.test[pleaf.inx] <- ifelse( cmod.ordinal <= cpoints[g], 2*i, 2*i+1)
              temp <- rebetQ(y, vi, mods = as.factor(cnode.test))
              if (temp[1] > Dev) {
                Dev <- temp[1]
                msplit <- paste(mods.names[k], "=", paste(names(mod.order[mod.order <= cpoints[g]]), collapse = "/"), collapse = " ")
                tres <- data.frame(Qb = temp[1], tau2 = temp[2],
                                   split = msplit, mod = mods.names[k], pleaf = as.numeric(nodes[j]))
                tcpt <- names(mod.order[mod.order <= cpoints[g]])
                tnode <- cnode.test
              }
            }
          }
        } else{
          chosenmod <- mods[pleaf.inx, k]
          cpoints <- sort(unique(chosenmod))
          if (length(cpoints) >= 2) {
            for (g in 1:(length(cpoints)-1)) {
              cnode.test <- cnode
              cnode.test[pleaf.inx] <- ifelse( chosenmod <= cpoints[g], 2*i, 2*i+1)
              temp <- rebetQ(y, vi, mods = as.factor(cnode.test))
              if (temp[1] > Dev) {
                Dev <- temp[1]
                tcpt <- cpoints[g]
                msplit <- paste(mods.names[k], "<=", tcpt, collapse = " ")
                tres <- data.frame(Qb = temp[1], tau2 = temp[2], split = msplit, mod = mods.names[k], pleaf = as.numeric(nodes[j]))
                tnode <- cnode.test
              }
            }
          }
        }

      }
    }
    if (is.null(tres)) {
      break
    } else {
      new.node <- tnode
      nodemark <- cbind(nodemark, new.node)
      res <- rbind(res, tres)
      cpt[[i]] <- tcpt
    }

  }
  list(tree = res, node.split = nodemark, cpt = cpt, data = mf)

}

#' A function to perform cross-validation for RE meta-tree
#'
#' @param mf the data.frame to grow the tree
#' @param maxL the maximum number of splits
#' @param n.fold the number of folds to perform the cross-validation
#' @return a cp table
#' @importFrom stats model.response
#' @keywords internal
REmrt.xvalid <- function(mf, maxL, n.fold = 10){
  N <- nrow(mf)
  if (n.fold > N | n.fold <=1 ) stop("n.fold is not acceptable")
  if (maxL > N) {
    warning("The maximum number of split is too large")
    maxL <- N-1
  }
  pred <- matrix(NA, nrow = N, ncol = maxL+1)
  inx <- sample(1:N)
  inx.xvalid <- c(round(seq(from = 1, by= N/n.fold, length.out = n.fold)), N+1)
  for (i in 1:n.fold){
    inx.test <- inx[inx.xvalid[i]:(inx.xvalid[i+1]-1)]
    test <- mf[inx.test, ]
    train <- mf[-inx.test, ]
    fit.train <- REmrt.fit0(train, maxL = maxL)
    yi.train <- model.response(fit.train$data)
    vi.train <- c(t(fit.train$data["(vi)"]))
    nsplt <- nrow(fit.train$tree)
    train.wy <- sapply(1:nsplt, function(x) yi.train/(vi.train+fit.train$tree[,2][x]))
    train.wts <-  sapply(1:nsplt, function(x) 1/(vi.train+fit.train$tree[,2][x]))
    train.pred <- sapply(1:nsplt, function(x) tapply(train.wy[, x], fit.train$node.split[, x], sum)/tapply(train.wts[, x], fit.train$node.split[, x], sum))
    test.nodes <- REmrt.prednode(fit.train, test)
    if (is.null(dim(test.nodes))) {
      pred[inx.test, 1:nsplt] <- train.pred[1]
    } else {
      test.pred <- lapply(1:nsplt, function(x) train.pred[[x]][as.character(test.nodes[ ,x])])
      test.pred.rmna <- sapply(1:nsplt, function(x){test.pred[[x]][is.na(test.pred[[x]])] <- sum(train.wy[, x])/sum(train.wts[, x]);test.pred[[x]]})
      pred[inx.test, 1:nsplt] <- test.pred.rmna
    }

  }
  pred <- pred[ ,!is.na(colSums(pred))]
  y <- model.response(mf)
  if (!is.null(dim(pred))) {
    x.error <- apply(pred,2, function(x) sum((y-x)^2)/sum((y-mean(y))^2))
    sdx.error <- apply(pred, 2, function(x)  sqrt(sum(((y-x)^2-mean((y-x)^2))^2))/sum((y-mean(y))^2))
  } else {
    x.error <- sum((y-pred)^2)/sum((y-mean(y))^2)
    sdx.error <- sqrt(sum(((y-pred)^2-mean((y-pred)^2))^2))/sum((y-mean(y))^2)
  }
  cbind(x.error, sdx.error)

}



