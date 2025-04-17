SSS_fe_Qb <- function(y, vi, xk, cpt, a){
  # y, vi, xk are the values within the parent leaf and have the same length! 
  i.r <- Sright(cpt, xk, a)
  #i.r <- Iright(cpt, xk) #for testing
  i.l <- 1 - i.r
  Q.l <- comp_Q_within(y, vi, i.l)
  Q.r <- comp_Q_within(y, vi, i.r)
  Q.l + Q.r
}


fe_cutoff_GS <- function(yi, vi, xk, minbucket){
  # given the moderator, the effecti size and the sampling variance
  # find the best split point
  # To do: add minbucekt and minsplit option
  n <- length(yi)
  xk.order <- order(xk)
  xk <- xk[xk.order]
  yi <- yi[xk.order]
  vi <- vi[xk.order]
  c.split <- (xk[-1] - xk[-n]) != 0
  if (minbucket > 1) {
    c.split[c(1:(minbucket-1), (n-minbucket+1):(n-1))] <- FALSE
  }
  if (all(c.split == FALSE)) {
    return(NULL)
  } else {
    wy <- yi / vi
    wy2 <- yi^2 / vi
    wts <- 1/vi
    wts2 <- wts^2
    cwy <- cumsum(wy[-n])
    cwy2 <- cumsum(wy2[-n])
    cwts <- cumsum(wts[-n])
    cwts2 <- cumsum(wts2[-n])
    Ql <- cwy2 - cwy^2/cwts
    Qr <- sum(wy2) - cwy2 - (sum(wy) - cwy)^2/(sum(wts) - cwts)
    Qrl <- Ql + Qr
    inx.star <- which(Qrl ==  min(Qrl[c.split]))
    res <- c((xk[inx.star]+xk[inx.star+1])/2, Ql[inx.star], Qr[inx.star])
    names(res) <- c("c.star", "Ql", "Qr")
    return(res)
  }
  
  wy <- yi / vi
  wy2 <- yi^2 / vi
  wts <- 1/vi
  wts2 <- wts^2
  cwy <- cumsum(wy[-n])
  cwy2 <- cumsum(wy2[-n])
  cwts <- cumsum(wts[-n])
  cwts2 <- cumsum(wts2[-n])
  Ql <- cwy2 - cwy^2/cwts
  Qr <- sum(wy2) - cwy2 - (sum(wy) - cwy)^2/(sum(wts) - cwts)
  Qrl <- Ql + Qr
  inx.star <- which(Qrl ==  min(Qrl))
  res <- c((xk[inx.star]+xk[inx.star+1])/2, Ql[inx.star], Qr[inx.star])
  names(res) <- c("c.star", "Ql", "Qr")
  return(res)
  
}

fe_cutoff_SSS <- function(y, vi, xk,
                          a = 50,
                          alpha.endcut = 0.02,
                          multi.start=T, n.starts=3){
  # xk is the value of moderator within the pleaf!
  n <- length(xk)
  sigma <- sd(xk); mu <- mean(xk)
  xk <- scale(xk)  # IMPORTANT TO STANDARDIZE x IN ORDER TO APPLY A CONSTANT a
  # FINDING THE SEARCH RANGE TO AVOID ENDCUT PREFERENCE PROBLEM
  LB <- quantile(xk, probs = alpha.endcut); UB <- quantile(xk, probs = 1-alpha.endcut); 
  if (multi.start==T) {
    B <- seq(LB, UB, length.out = n.starts)
    Q.min <- Inf
    for (b in 2:n.starts) {
      OPT <- optimize(SSS_fe_Qb, lower=B[b-1], upper=B[b], maximum=FALSE, 
                      y = y, vi = vi, xk = xk, a = a)
      if (OPT$objective < Q.min) {
        Q.min <- OPT$objective
        cstar <- OPT$minimum
      }
    }
  } else {
    cstar <- optimize(SSS_fe_Qb, lower=LB, upper=UB, maximum=F, 
                      a = a, y = y, xk = xk, vi = vi)$minimum
  }
  i.r <- Iright(cstar, xk)
  i.l <- 1 - i.r
  Q.l <- comp_Q_within(y, vi, i.l)
  Q.r <- comp_Q_within(y, vi, i.r)
  if(any(is.na(c(Q.l, Q.r)))) {
    return(NULL)
  } else {
    d.l <- sum(y*i.l/vi)/sum(i.l/vi)
    d.r <- sum(y*i.r/vi)/sum(i.r/vi)
    se.l <- sqrt(1/sum(i.l/vi))
    se.r <- sqrt(1/sum(i.r/vi))
    cstar <- cstar*sigma + mu	# TRANSFORM BACK
    return(c(cstar, Q.l, Q.r, d.l, d.r, se.l, se.r))
  }
}


