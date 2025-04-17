find_two_cutoffs <- function(y, vi, xk1, xk2, a, assump,
                             alpha.endcut = 0.02,
                             multi.start=T, n.starts=3){
  fn1 <- function(x) {
    compute_Qb_two_splits(y, vi, xk1, xk2, x[1], x[2], a,
                          splitfrom = "left", assump)
  }
  fn2 <- function(x) {
    compute_Qb_two_splits(y, vi, xk1, xk2, x[1], x[2], a,
                          splitfrom = "right", assump)
  }
  LB1 <- quantile(xk1, probs = alpha.endcut); UB1 <- quantile(xk1, probs =1-alpha.endcut); 
  LB2 <- quantile(xk2, probs = alpha.endcut); UB2 <- quantile(xk2, probs =1-alpha.endcut); 
  if(multi.start == T) {
    Q.min <- 1e15
    for (b in 1:n.starts) {
      rstart1 <- runif(1, LB1, UB1)
      rstart2 <- runif(1, LB2, UB2)
      #print(rbind(c(rstart1, LB1, UB1), c(rstart2, LB2, UB2)))
      res1 <- optim(c(rstart1, rstart2), fn1, gr = NULL, method = "L-BFGS-B",
                    lower = c(LB1,LB2), upper =c(UB1,UB2))
      res2 <- optim(c(rstart1, rstart2), fn2, gr = NULL, method = "L-BFGS-B",
                    lower = c(LB1,LB2), upper =c(UB1,UB2))
      if (res1$value < res2$value) {
        temp <- c(2, res1$value, res1$par)
      } else {
        temp <- c(3, res2$value, res2$par)
      }
      if (temp[2] < Q.min) {
        Q.min <- temp[2]
        res <- temp
      }
    }
    return(res)
    
  } else {
    rstart1 <- runif(1, LB1, UB1)
    rstart2 <- runif(1, LB2, UB2)
    res1 <- optim(c(rstart1, rstart2), fn1, gr = NULL, method = "L-BFGS-B",
                  lower = c(LB1,LB2), upper =c(UB1,UB2))
    res2 <- optim(c(rstart1, rstart2), fn2, gr = NULL, method = "L-BFGS-B",
                  lower = c(LB1,LB2), upper =c(UB1,UB2))
    if (res1$value < res2$value) {
      return(c(2, res1$value, res1$par))
    } else {
      return(c(3, res2$value, res2$par))
    }
  }
}


find_two_cutoffs_on_one_var <- function(y, vi, xk1, a, assump,
                                        alpha.endcut = 0.02,
                                        multi.start=T, n.starts=3){
  LB1 <- quantile(xk1, probs = alpha.endcut); UB1 <- quantile(xk1, probs =1-alpha.endcut)
  if (UB1 - LB1 < alpha.endcut) {
    return(c(NA, Inf ,NA))
  } else {
    fn2 <- function(x) {
      compute_Qb_two_splits(y, vi, xk1, xk1, x[1], x[2], a,
                            splitfrom = "right", assump)
    } 
    if(multi.start == T) {
      Q.min <- 1e15
      for (b in 1:n.starts) {
        rstart1 <- runif(1, LB1, UB1 - alpha.endcut)
        rstart2 <- runif(1, rstart1 + alpha.endcut, UB1)
        res2 <- constrOptim(c(rstart1, rstart2), f = fn2, grad = NULL,
                            ui = cbind(c(-1,1,0,-1,0),c(1,0,1,0,-1)),
                            ci = c(alpha.endcut, LB1, LB1, -UB1, -UB1))
        temp <- c(3, res2$value, res2$par)
        if (temp[2] < Q.min) {
          Q.min <- temp[2]
          res <- temp
        }
      }
      return(res)
      
    } else {
      rstart1 <- runif(1, LB1, UB1 - alpha.endcut)
      rstart2 <- runif(1, rstart1 + alpha.endcut, UB1)
      res2 <- constrOptim(c(rstart1, rstart2), f = fn2, grad = NULL,
                          ui = cbind(c(-1,1,0,-1,0),c(1,0,1,0,-1)),
                          ci = c(alpha.endcut, LB1, LB1, -UB1, -UB1))
      return(c(3, res2$value, res2$par))
      
    }
  }
  
}

compute_Qb_two_splits <- function(y, vi, xk1, xk2, cpt1, cpt2, a,
                                  splitfrom, assump){
  df <- length(y) - 3
  i1.r <-Sright(cpt1, xk1, a)
  i2.r <- Sright(cpt2, xk2, a)
  i1.l <- 1-i1.r
  i2.l <- 1-i2.r
  if (splitfrom == "right") {
    i.l <- i1.l
    i.m <- i1.r*i2.l
    i.r <- i1.r*i2.r
  } else {
    i.l <- i1.l*i2.l
    i.m <- i1.l*i2.r
    i.r <- i1.r
  }
  
  Q.l <- comp_Q_within(y, vi, i.l)
  Q.m <- comp_Q_within(y, vi, i.m)
  Q.r <- comp_Q_within(y, vi, i.r)
  if (assump == "re"){
    C.l <- comp_C(vi, i.l)
    C.m <- comp_C(vi, i.m)
    C.r <- comp_C(vi, i.r)
    tau2 <-  (sum(c(Q.l, Q.m, Q.r), na.rm = T) - df)/(sum(1/vi) - sum(c(C.l, C.m, C.r), na.rm = T))
    tau2 <- max(tau2,0)
    vi.star <- vi + tau2
    Qstar.l <- comp_Q_within(y, vi.star, i.l)
    Qstar.m <- comp_Q_within(y, vi.star, i.m)
    Qstar.r <- comp_Q_within(y, vi.star, i.r)
    Q.between <- sum(y^2/vi.star) - (sum(y/vi.star))^2/sum(1/vi.star) - sum(c(Qstar.l, Qstar.r, Qstar.m), na.rm = T)
    
  } else { # checked
    Q.between <- sum(y^2/vi) - (sum(y/vi))^2/sum(1/vi) - sum(c(Q.l, Q.m, Q.r), na.rm = T)
  }
  -Q.between
}