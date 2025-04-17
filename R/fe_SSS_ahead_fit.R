fe_SSS_LH_first_two <- function(y, vi, mods, assump = "fe", 
                                a, alpha.endcut,
                                multi.start, n.starts){
  Dev = Inf
  for (k1 in 1:ncol(mods)){
    xk1 <- mods[, k1]
    sigma1 <- sd(xk1); mu1 <- mean(xk1)
    xk1 <- scale(xk1)
    for (k2 in 1:ncol(mods)){
      if (k1 == k2) {
        temp <- find_two_cutoffs_on_one_var(y, vi, xk1, a, assump,
                                            alpha.endcut, multi.start, n.starts)
        temp[3:4] <- temp[3:4]*sigma1 + mu1
      } else {
        xk2 <- mods[, k2]
        sigma2 <- sd(xk2); mu2 <- mean(xk2)
        xk2 <- scale(xk2)
        temp <- find_two_cutoffs(y, vi, xk1, xk2, a, assump,
                                 alpha.endcut, multi.start, n.starts)
        temp[3] <- temp[3]*sigma1 + mu1
        temp[4] <- temp[4]*sigma2 + mu2
      }
      if(temp[2] < Dev) {
        Dev <- temp[2]
        pleaf <- as.integer(temp[1])
        cpt1 <- temp[3]
        cpt2 <- temp[4]
        mod1 <- k1 
        mod2 <- k2
      }
    }
    
  }
  if(is.infinite(Dev)) stop("no possible split points for look-ahead")
  c(pleaf, mod1, mod2, cpt1, cpt2)
} 

fe_SSS_LH_continue <- function(y, vi, mods,
                               pleaf, mod1, mod2, cpt1, cpt2,
                               minsplit, maxdepth, cp,
                               a, alpha.endcut, multi.start, n.starts) {
  mods.names <- colnames(mods)
  pnode <- c(1, pleaf)
  msplit1 <- paste(mods.names[mod1], "<", cpt1, collapse = " ")
  msplit2 <- paste(mods.names[mod2], "<", cpt2, collapse = " ")
  # msplit1 <- paste(mods.names[mod1], "<", round(cpt1,2), collapse = " ")
  # msplit2 <- paste(mods.names[mod2], "<", round(cpt2,2), collapse = " ")
  childL <- which(mods[, mod1] < cpt1)
  childR <- which(mods[, mod1] >= cpt1)
  d1 <- sum(y/vi)/sum(1/vi)
  d2 <- sum(y[childL]/vi[childL])/sum(1/vi[childL])
  d3 <- sum(y[childR]/vi[childR])/sum(1/vi[childR])
  Q1 <- sum((y - d1)^2/vi)
  Q2 <- sum((y[childL] - d2)^2/vi[childL])
  Q3 <- sum((y[childR] - d3)^2/vi[childR])
  delQ1 <- Q1- Q2 -Q3
  if(pleaf == 2) {
    grandL <- childL[which(mods[childL, mod2] < cpt2)]
    grandM <- childL[which(mods[childL, mod2] >= cpt2)]
    grandR <- childR
    #  sort(c(grandL, grandM, grandR)) # to check
    d4 <-  sum(y[grandL]/vi[grandL])/sum(1/vi[grandL])
    d5 <- sum(y[grandM]/vi[grandM])/sum(1/vi[grandM])
    Q4 <- sum((y[grandL] - d4)^2/vi[grandL])
    Q5 <- sum((y[grandM] - d5)^2/vi[grandM])
    delQ2 <- Q2 - Q4 - Q5
    delQ <- c(delQ1, delQ2)
    d <- c(d1, d2)
    wt <- c(sum(1/vi), sum(1/vi[childL]))
    var <- c(mods.names[mod1], mods.names[mod2])
    split <- c(msplit1, msplit2)
    No. <- c(length(y), length(childL))
    left.No. <- c(length(childL), length(grandL))
    right.No. <- c(length(childR), length(grandM))
    Q <- c(Q1, Q2)
    leafL <- 4; leafM <- 5; leafR <- 3
    inx0 <- list(childL, grandL, grandM, grandR)
    inx0 <- setNames(inx0, c(pleaf, leafL, leafM, leafR))
    stats0 <- cbind(c(leafL, leafM, leafR),
                    c(length(grandL),length(grandM), length(grandR)),
                    c(Q4, Q5, Q3),
                    c(d4, d5, d3),
                    c(sum(1/vi[grandL]), sum(1/vi[grandM]), sum(1/vi[grandR])))
  } else {
    grandL <- childL
    grandM <- childR[which(mods[childR, mod2] < cpt2)]
    grandR <- childR [which(mods[childR, mod2] >= cpt2)]
    # sort(c(grandL, grandM, grandR)) # to check
    d6 <-  sum(y[grandM]/vi[grandM])/sum(1/vi[grandM])
    d7 <- sum(y[grandR]/vi[grandR])/sum(1/vi[grandR])
    Q6 <- sum((y[grandM] - d6)^2/vi[grandM])
    Q7 <- sum((y[grandR] - d7)^2/vi[grandR])
    delQ2 <- Q3 - Q6 - Q7
    delQ <- c(delQ1, delQ2)
    d <- c(d1, d3)
    wt <- c(sum(1/vi), sum(1/vi[childR]))
    var <- c(mods.names[mod1], mods.names[mod2])
    split <- c(msplit1, msplit2)
    No. <- c(length(y), length(childR))
    left.No. <- c(length(childL), length(grandM))
    right.No. <- c(length(childR), length(grandR))
    Q <- c(Q1, Q3)
    leafL <- 2; leafM <- 6; leafR <- 7
    inx0 <- list(childR, grandL, grandM, grandR)
    inx0 <- setNames(inx0, c(pleaf, leafL, leafM, leafR))
    stats0 <- cbind(c(leafL, leafM, leafR),
                    c(length(grandL),length(grandM), length(grandR)),
                    c(Q2, Q6, Q7),
                    c(d2, d6, d7),
                    c(sum(1/vi[grandL]), sum(1/vi[grandM]), sum(1/vi[grandR])))
  }
  frame0 <- data.frame(pnode, No., left.No., right.No., Q,
                       delQ, d, wt, var, split, stringsAsFactors = FALSE)
  cpt <- list()
  cpt[[1]] <- cpt1
  cpt[[2]] <- cpt2
  
  resL <- FEpartitionSSS(y[grandL], vi[grandL], mods[grandL, ], 
                         pnode = leafL,
                         inx.within.node = grandL,
                         split.candidate = mods.names[apply(mods[grandL, ], 2, function(x) length(unique(x))) >= 2],  
                         minsplit, maxdepth, cp = 0,
                         a, alpha.endcut, multi.start, n.starts)
  resM <- FEpartitionSSS(y[grandM], vi[grandM], mods[grandM, ], 
                         pnode = leafM,
                         inx.within.node = grandM,
                         split.candidate = mods.names[apply(mods[grandM, ], 2, function(x) length(unique(x))) >= 2],  
                         minsplit, maxdepth, cp = 0,
                         a, alpha.endcut, multi.start, n.starts)
  resR <- FEpartitionSSS(y[grandR], vi[grandR], mods[grandR, ], 
                         pnode = leafR,
                         inx.within.node = grandR,
                         split.candidate = mods.names[apply(mods[grandR, ], 2, function(x) length(unique(x))) >= 2],  
                         minsplit, maxdepth, cp = 0,
                         a, alpha.endcut, multi.start, n.starts)
  res <- list()
  res[[1]] <- rbind(frame0, resL[[1]], resM[[1]], resR[[1]])
  res[[2]] <- do.call(c, list(inx0, resL[[2]], resM[[2]], resR[[2]]))
  res[[3]] <- do.call(c, list(cpt, resL[[3]], resM[[3]], resR[[3]]))
  res[[4]] <- do.call(rbind, list(stats0, resL[[4]], resM[[4]], resR[[4]]))
  res.path <- fe_find_pruning_path(res)
  cps <- unique(res.path[[1]]$cp)
  if (any(cps <= cp)){
    cp.lb <- max(cps[cps <= cp])
  } else {
    cp.lb <- 0
  }
  res <- fe_prune(res.path, cp.lb)
  res
  # sort(unlist(res[[2]])) #check
}


fe_SSS_LH_fit <- function(mf, minsplit = 6, maxdepth = 4, cp,
                                     a = 50,
                                     alpha.endcut = 0.02,
                                     multi.start=T, n.starts=3){
  y <- model.response(mf)
  vi <- c(t(mf["(vi)"]))
  mods.names <-  labels(terms(mf))
  mods <- mf[mods.names]
  first.two.splits <- fe_SSS_LH_first_two(y, vi, mods, assump = "fe", 
                                            a, alpha.endcut,
                                            multi.start, n.starts)
  pleaf <- as.integer(first.two.splits[1])
  mod1 <- as.integer(first.two.splits[2])
  mod2 <- as.integer(first.two.splits[3])
  cpt1 <- first.two.splits[4]
  cpt2 <- first.two.splits[5]
  fe_SSS_LH_continue(y, vi, mods,
                     pleaf, mod1, mod2, cpt1, cpt2,
                     minsplit, maxdepth, cp,
                     a, alpha.endcut, multi.start, n.starts)
  
}




    