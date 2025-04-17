FEpartitionSSS<- function(y, vi, mods, pnode, inx.within.node,
                          split.candidate,  minsplit, maxdepth, cp,
                          a, alpha.endcut, multi.start, n.starts) {
  # y: the effect size
  # vi: the sample variance 
  # mods: the data.frame of moderators
  # pnode: an integer indicating the label of the root node, usually equal to one
  # inx.within.node: rows within the root node
  # split.candidate: the candidate for splitting moderators
  # 
  
  a = a;
  alpha.endcut = alpha.endcut
  multi.start = multi.start
  n.starts = n.starts
  
  # base case
  if (length(y) <  minsplit | length(split.candidate) == 0 | pnode >= 2^maxdepth) {
    terminalNode <- c(pnode, length(y), sum(y^2/vi) - (sum(y/vi))^2/sum(1/vi),
                      sum(y/vi)/sum(1/vi), sum(1/vi))
    return(list(NULL, NULL, NULL, terminalNode))
  } else {
    # print(pnode)
    # print(y)
    # print(split.candidate)
    # print(nodes)
    Q <- sum(y^2/vi) - (sum(y/vi))^2/sum(1/vi)
    d <- sum(y/vi)/sum(1/vi)
    wt <- sum(1/vi)
    Dev <- -Inf
    for (k in split.candidate){
      xk <- mods[, k]
      if (is.numeric(xk)) {
        xk.ord <- xk
      } else {
        xk.rank <- rank(tapply(y, xk, mean))
        xk.ord <- xk.rank[as.character(xk)]
      }
      temp <- fe_cutoff_SSS(y, vi, xk.ord, a = a ,
                            alpha.endcut = alpha.endcut,
                            multi.start = multi.start, 
                            n.starts = n.starts)
      if (is.null(temp)) {
        Dev.new <- -Inf
      } else {
        Dev.new <- Q - sum(temp[2:3])
      }
      if (Dev.new > Dev) {
        Dev <- Dev.new
        if (is.numeric(xk)) {
          c.star <- temp[1]
          msplit <- paste(k, "<", c.star, collapse = " ")
          # msplit_round <- paste(k, "<", round(c.star, 2), collapse = " ")
        } else {
          c.star <- names(xk.rank[xk.rank <= temp[1]])
          msplit <- paste(k, "=", paste(c.star, collapse = "/"), collapse = " ")
          
        }
        mod <- k
        inx.left <- xk.ord < temp[1]
        inx.right <- xk.ord >= temp[1]
      }
    }
    if (Dev < cp) {
      terminalNode <- c(pnode, length(y), sum(y^2/vi) - (sum(y/vi))^2/sum(1/vi),
                        sum(y/vi)/sum(1/vi), sum(1/vi))
      return(list(NULL, NULL, NULL, terminalNode))
    } else {
      inx.within.node.left <- inx.within.node[inx.left]
      inx.within.node.right <- inx.within.node[inx.right]
      nodes <- do.call(c,list(setNames(list(inx.within.node.left), pnode*2),
                              setNames(list(inx.within.node.right), pnode*2 + 1)))
      # print(apply(mods[inx.left, ], 2, function(x) length(unique(x))) >= 2)
      # print(apply(mods[inx.right, ], 2, function(x) length(unique(x))) >= 2)
      if (length(split.candidate) == 1) {
        new.split.candidate.left <-  new.split.candidate.right <-  split.candidate
      } else {
        new.split.candidate.left <- split.candidate[apply(mods[inx.left, split.candidate], 2, function(x) length(unique(x))) >= 2]
        new.split.candidate.right <- split.candidate[apply(mods[inx.right, split.candidate], 2, function(x) length(unique(x))) >= 2]
      }
      res.left <- FEpartitionSSS(y[inx.left], vi[inx.left], mods[inx.left, ],
                               pnode*2, inx.within.node.left,
                               new.split.candidate.left,  minsplit, maxdepth, cp,
                               a, alpha.endcut, multi.start, n.starts)
      res.right <- FEpartitionSSS(y[inx.right], vi[inx.right], mods[inx.right, ],
                                pnode*2 + 1, inx.within.node.right, 
                                new.split.candidate.right,  minsplit, maxdepth, cp,
                                a, alpha.endcut, multi.start, n.starts)
      list(rbind(data.frame(pnode = pnode, No. = length(y),
                            left.No. = sum(inx.left), right.No. = sum(inx.right),
                            Q = Q, delQ = Dev, d = d, wt = wt, var = mod, split = msplit,
                            stringsAsFactors = FALSE),
                 res.left[[1]],
                 res.right[[1]]), 
           do.call(c, list(nodes, res.left[[2]], res.right[[2]])),
           do.call(c, list(list(c.star), res.left[[3]], res.right[[3]])),
           rbind(res.left[[4]], res.right[[4]]))
      
    }
  }
}

FEpartitionSSS.wrapper <- function(mf, minsplit = 6, maxdepth = 4, cp,
                                   a = 50,
                                   alpha.endcut = 0.02,
                                   multi.start=T, n.starts=3){
  # Returns:
  # list1: the tree frame
  # list2: the distributed rows in terminal nodes
  # list3: the split points
  # list4: 1st col terminal nodes, 2nd col # of studies,
  #        3rd col effect size, 4th col Q, 5th col weights
  y <- model.response(mf)
  vi <- c(t(mf["(vi)"]))
  split.candidate <- labels(terms(mf))
  mods <- mf[split.candidate]
  inx.within.node <- 1:nrow(mf)
  split.candidate <- split.candidate[apply(mods, 2, function(x) length(unique(x))) >= 2]
  totalQ <- sum(y^2/vi) - (sum(y/vi))^2/sum(1/vi)
  res <- FEpartitionSSS(y, vi, mods, pnode = 1, inx.within.node = inx.within.node,
                        split.candidate = split.candidate,  minsplit = minsplit, 
                        maxdepth = maxdepth, cp = 0,
                        a = a, alpha.endcut = alpha.endcut, 
                        multi.start = multi.start, n.starts = n.starts)
  if(is.null(res[[1]])) {res} else {
    res.path <- fe_find_pruning_path(res)
    cps <- unique(res.path[[1]]$cp)
    if (any(cps <= cp)){
      cp.lb <- max(cps[cps <= cp])
    } else {
      cp.lb <- 0
    }
    
    res <- fe_prune(res.path, cp.lb) 
    #res[[1]]$cp[res[[1]]$cp == cp.lb] <- cp
    res    
  }
  
}


