FEmrt_SSS <- function(mf, cp = 0.01, maxdepth = 10L, minsplit = 6L,
                      n.fold = 10L, lookahead = FALSE, 
                      a = 50,
                      alpha.endcut = 0.02,
                      multi.start=T, n.starts=3, c.pruning = 0.5, ...){
  if (alpha.endcut < 0 | alpha.endcut > 1) {
    stop("alpha.endcut should range between 0 and 1")
  }
  if (!is.logical(lookahead)) {
    stop("lookahead should be TRUE or FALSE")
  }
  if(!lookahead) {
    initial.tree <- FEpartitionSSS.wrapper(mf, minsplit, maxdepth, cp,
                                           a, alpha.endcut, multi.start, n.starts)
    cptable <- FExvalid(mf = mf, initial.tree = initial.tree, 
                        n.fold = n.fold, minsplit = minsplit, 
                        maxdepth = maxdepth, cp = cp,
                        a = a, alpha.endcut = alpha.endcut, multi.start = multi.start,
                        n.starts = n.starts, lookahead = FALSE)
  } else {
    initial.tree <- fe_SSS_LH_fit(mf, minsplit = minsplit, maxdepth = maxdepth,
                                  cp = cp, a = a, alpha.endcut = alpha.endcut,
                                  multi.start = multi.start, n.starts = n.starts)
    cptable <- FExvalid(mf = mf, initial.tree = initial.tree, 
                        n.fold = n.fold, minsplit = minsplit, 
                        maxdepth = maxdepth, cp = cp,
                        a = a, alpha.endcut = alpha.endcut, multi.start = multi.start,
                        n.starts = n.starts, lookahead = TRUE)
    
  }
  mindex <- which.min(cptable[,4]) 
  cp.minse <- cptable[mindex,4] + c.pruning*cptable[mindex,5]  
  cp.row <- min(which(cptable[,4]<= cp.minse))  
  cp.take <- cptable[cp.row, 1] 
  prunedtree <- fe_prune(initial.tree, cp.take)
  frame <- prunedtree[[1]]
  if (nrow(frame) > 0) {
    terminalNodes <- prunedtree[[4]]
    n <- terminalNodes[, 2]
    names(n) <- terminalNodes[, 1]
    Qw <- terminalNodes[, 3]
    g <- terminalNodes[, 4]
    names(g) <- terminalNodes[, 1]
    Qb <- frame$Q[1] - sum(Qw)
    df <- nrow(frame)
    pval.Qb <- pchisq(Qb, df, lower.tail = FALSE)
    se <- sqrt(1/terminalNodes[ ,5])
    zval <- g/se
    pval <- pnorm(abs(zval),lower.tail=FALSE)*2
    ci.lb <- g - qnorm(0.975)*se
    ci.ub <- g + qnorm(0.975)*se
    moderators <- unique(frame$var)
    node.membership <- numeric(nrow(mf))
    for (k in as.character(prunedtree[[4]][,1])){
      node.membership[prunedtree[[2]][[k]]] <- as.numeric(k)
    }
    cpt <- prunedtree[[3]]
    res <- list(tree =frame, n = n, moderators = moderators , Qb = Qb, df = df, pval.Qb = pval.Qb,
                Qw = Qw, g = g, se = se, zval =zval, pval = pval, ci.lb = ci.lb,
                ci.ub = ci.ub, cptable = cptable, node.membership = node.membership,
                initial.tree = initial.tree[[1]], cpt = cpt,
                data = mf)
    
  } else{
    warning("no moderator effect was detected")
    y <- model.response(mf)
    v <- c(t(mf["(vi)"]))
    n <- length(y)
    g <- sum(y/v)/sum(1/v)
    Q <- sum((y-g)^2/v)
    df <- n - 1
    pval.Q <- pchisq(Q, df, lower.tail = FALSE)
    se <- 1/sqrt(sum(1/v))
    zval <- g/se
    pval <- pnorm(abs(zval), lower.tail=FALSE)*2
    ci.lb <- g - qnorm(0.975)*se
    ci.ub <- g + qnorm(0.975)*se
    res <- list(n = n ,  Q = Q,
                df = df, pval.Q = pval.Q, g = g, se = se, zval = zval,
                pval = pval, ci.lb = ci.lb, ci.ub = ci.ub, cptable = cptable, 
                initial.tree = initial.tree[[1]], data = mf)
  }
  res
}

