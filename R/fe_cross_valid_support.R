FExvalid <- function(mf, initial.tree, n.fold, minsplit = 6, maxdepth = 4, cp = 0.01,
                     a = 50,
                     alpha.endcut = 0.02,
                     multi.start=T, n.starts=3,
                     lookahead = F) {
  N <- nrow(mf)
  cps <-  sort(unique(initial.tree[[1]]$cp))
  if (n.fold > N | n.fold <=1 ) stop("n.fold is not acceptable")
  inx <- sample(1:N)
  inx.xvalid <- as.numeric(cut(1:N, n.fold))
  predY <- matrix(nrow = N, ncol = length(cps) + 1)
  if (lookahead == TRUE) func <- fe_SSS_LH_fit else func <- FEpartitionSSS.wrapper
  for (i in 1:n.fold){
    inx.test <- inx[inx.xvalid == i]
    test <- mf[inx.test, ]
    train <- mf[-inx.test, ]
    fit.train <- do.call(func, list(mf = train, minsplit = minsplit, 
                                    maxdepth = maxdepth, cp = cp,
                                    a = a, alpha.endcut = alpha.endcut, 
                                    multi.start = multi.start, 
                                    n.starts = n.starts)) 
    if (is.null(fit.train[[1]])){
      predY[inx.test, ] <- fit.train[[4]][1, 4]
    } else {
      allLeaves <- c(fit.train[[1]]$pnode, fit.train[[4]][,1])
      allWmeans <- c(fit.train[[1]]$d, fit.train[[4]][,4])
      prunePath <- fe_prune_path(fit.train[[1]], cps)
      testnode <- fe_predict_node(fit.train[[1]], fit.train[[3]], rep(1,nrow(test)), test, train)
      testY <- fe_predictY(fit.train[[1]], fit.train[[3]], prunePath, test, train,
                                  allLeaves, allWmeans)
      testY[is.na(testY)] <- allWmeans[which(allLeaves == 1)]
      testY[, 1] <- allWmeans[which(allLeaves == 1)] # force no split for the largest complexity parameter
      
      predY[inx.test, ] <- testY
    }
  } 
  y <- model.response(mf)
  wt <- 1/mf$`(vi)`/sum(1/mf$`(vi)`)
  split.left <- sapply(fe_prune_path(initial.tree[[1]], cps), length)
  decreasedQ <- sapply(fe_prune_path(initial.tree[[1]], cps),
                       function(x) sum(initial.tree[[1]]$delQ[match(x, initial.tree[[1]]$pnode)]))
  rel.error <- (initial.tree[[1]]$Q[1] - cumsum(decreasedQ[length(decreasedQ):1]))/initial.tree[[1]]$Q[1] 
  nsplit <- cumsum(split.left[length(split.left):1])
  CP <- c(cps[length(cps):1], cp)
  x.error <- apply(predY, 2, function(x) sum(wt*(y-x)^2)/initial.tree[[1]]$Q[1]*sum(1/mf$`(vi)`))
  sdx.error <- apply(predY, 2, function(x) sqrt(sum(wt*((y-x)^2-sum(wt*(y-x)^2))^2))/initial.tree[[1]]$Q[1]*sum(1/mf$`(vi)`)/sqrt(N))
  cbind(CP, nsplit, rel.error, x.error, sdx.error)
}
#FExvalid(mf, FEpartitionSSS.wrapper(mf, cp = 0.01), n.fold = 10)
#FExvalid(mf, fe_SSS_LH_fit(mf, cp = 0.01), n.fold = 10, lookahead = T)
