#' A function to find the split point
#'
#' @param g the effect size
#' @param vi the sampling variance
#' @param x the splitting moderator 
#' @param minbucket the minimum number of the studies in a terminal node
#' @param inx.s indicates whether a study belongs to the candidate parent leaf
#' @param cnode the terminal nodes that the studies belong to in the current tree
#' @useDynLib metacart
#' @return a vector including the split point, Q, and tau2
#' @keywords internal
re.cutoff_cpp <- function(g, vi, x, inx.s, cnode, minbucket) {
  n <- sum(inx.s)
  xinx.order <- order(x)
  x.sort = x[xinx.order]
  c.split <- (x.sort[-1] - x.sort[-n]) != 0
  if (minbucket > 1) {
    c.split[c(1:(minbucket - 1), (n - minbucket + 1):(n - 1))] <- FALSE
  }
  if (all(c.split == FALSE)) {
    return(NULL)
  }
  else {
    g.sort <- g[inx.s][xinx.order]
    vi.sort <- vi[inx.s][xinx.order]
    inx.unsplit <- inx.s == FALSE
    g.unsplit <- g[inx.unsplit]
    vi.unsplit <- vi[inx.unsplit]
    cnode.unsplit <- cnode[inx.unsplit]
    tau2 <- .compute_tau_(g.unsplit, vi.unsplit, cnode.unsplit, 
                          unique(cnode.unsplit), g.sort, vi.sort)
    tau2 <- pmax(0, tau2)
    qb.star <- .compute_re_Q_(g.unsplit, vi.unsplit, cnode.unsplit, 
                              tau2, unique(cnode.unsplit), g.sort, vi.sort)
    inx.star <- which(qb.star == max(qb.star[c.split]))
    cstar <- x.sort[inx.star]
    res <- c(cstar, qb.star[inx.star], tau2[inx.star])
    res
  }
  
} 

re_cutoff_SSS <- function(y, vi, xk, cnode, pleaf, 
                          a = 50,
                          alpha.endcut = 0.02,
                          multi.start=T, n.starts=3){
  # xk is the value of moderator within the pleaf!
  n <- length(xk)
  # FINDING THE SEARCH RANGE TO AVOID ENDCUT PREFERENCE PROBLEM
  LB <- quantile(xk, probs = alpha.endcut); UB <- quantile(xk, probs = 1-alpha.endcut); 
  if (multi.start==T) {
    (B <- seq(LB, UB, length.out=n.starts))
    Q.min <- Inf
    for (b in 2:n.starts) {
      OPT <- optimize(SSS_re_Qb, lower=B[b-1], upper=B[b], maximum=FALSE, 
                      y = y, vi = vi, xk = xk, cnode = cnode, 
                      pleaf = pleaf, a = a)
      if (OPT$objective < Q.min) {
        Q.min <- OPT$objective
        cstar <- OPT$minimum
      }
    }
  } else {
    cstar <- optimize(SSS_re_Qb, lower=LB, upper=UB, maximum=F, 
                      a = a, y = y, xk = xk, vi = vi,  cnode = cnode, pleaf = pleaf)$minimum
  }
  
  return(cstar)
}

SSS_re_Qb <- function(y, vi, xk, cnode, pleaf, cpt, a){
  df <- length(y) - length(unique(cnode)) - 1
  inx.unsplit <- cnode!= pleaf
  inx.split <- cnode == pleaf
  y_pleaf <- y[inx.split]
  vi_pleaf <- vi[inx.split]
  y_unsplit <- y[inx.unsplit]
  vi_unsplit <- vi[inx.unsplit]
  cnode_unsplit <- cnode[inx.unsplit]
  i.r <- Sright(cpt, xk, a)
  #i.r <- Iright(cpt, xk) for testing
  i.l <- 1 - i.r
  QC.unsplit <- compute_left_(y_unsplit, vi_unsplit, 
                              cnode_unsplit, unique(cnode_unsplit))
  Q.l <- comp_Q_within(y_pleaf, vi_pleaf, i.l)
  Q.r <- comp_Q_within(y_pleaf, vi_pleaf, i.r)
  C.l <- comp_C(vi_pleaf, i.l)
  C.r <- comp_C(vi_pleaf, i.r)
  tau2 <- comp_tau2(vi, Q.l + Q.r + QC.unsplit[1],
                    C.l + C.r + QC.unsplit[2], df)
  vi.star <- vi + tau2  
  vi.star_pleaf <- vi.star[inx.split]
  vi.star_unsplit <- vi.star[inx.unsplit]
  Qstar.unsplit <- compute_left_(y_unsplit, vi.star_unsplit, cnode_unsplit, unique(cnode_unsplit))[1]
  Qstar.l <- comp_Q_within(y_pleaf, vi.star_pleaf, i.l)
  Qstar.r <- comp_Q_within(y_pleaf, vi.star_pleaf, i.r)
  Q.between <- sum(y^2/vi.star) - (sum(y/vi.star))^2/sum(1/vi.star) - Qstar.l - Qstar.r -Qstar.unsplit
  #c(Q.between, tau2)
  -Q.between
}



