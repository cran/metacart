fe_find_pruning_path <- function(tree) {
  treeFrame <- tree[[1]]
  PrunedFrame <- treeFrame
  cp <- numeric(0)
  optimalSubtree <- list()
  treeFrame$cp <- NA
  while(nrow(PrunedFrame) >= 1) {
    complexity_decrease_per_node <- sapply(PrunedFrame$pnode,function(x) {
      temp <- complexity_decrease_(x, PrunedFrame$delQ/PrunedFrame$Q[1], PrunedFrame$pnode) 
      temp[1]/temp[2]
    })
    cp <- min(complexity_decrease_per_node)
    inx.weakest <- which(complexity_decrease_per_node == cp)
    nodeToBePruned <- find_offsprings_(PrunedFrame$pnode[inx.weakest],
                                       PrunedFrame$pnode) 
    toBeRemoved <- match(nodeToBePruned, PrunedFrame$pnode)
    PrunedFrame <- PrunedFrame[-toBeRemoved , ]
    cpadd <- match(nodeToBePruned, treeFrame$pnode)
    treeFrame$cp[cpadd] <- cp
  }
  if (any(is.na(treeFrame$cp))) warning("results contain missing values in complexity parameters")
  tree[[1]] <- treeFrame
  tree
}

fe_prune <- function(tree, cp){
  #stop('something is wrong with the reverse of Q and d!')
  treeFrame <- tree[[1]]
  indices <- tree[[2]]
  cpt <- tree[[3]]
  terminalNodes <- tree[[4]]
  node <- c(treeFrame$pnode, terminalNodes[,1])
  No. <- c(treeFrame$No., terminalNodes[,2])
  Q <- c(treeFrame$Q, terminalNodes[,3])
  d <- c(treeFrame$d, terminalNodes[,4])
  wts <- c(treeFrame$wt, terminalNodes[ ,5])
  inx1 <- which(treeFrame$cp > cp)
  nodes.kp <- treeFrame$pnode[inx1]
  treeFrame <- treeFrame[inx1, ]
  cpt <- cpt[inx1]
  allnodes <- as.numeric(names(indices))
  terNodes <- setdiff(find_children_vec(nodes.kp, c(1,allnodes)), nodes.kp)
  indices <- indices[as.character(c(terNodes, nodes.kp[nodes.kp!=1]))]
  inx4 <- match(terNodes, node)
  terminalNodes <- cbind(node[inx4], No.[inx4], Q[inx4], d[inx4], wts[inx4])
  list(treeFrame, indices, cpt, terminalNodes)
}

fe_prune_path <- function(treeFrame, cps){
  # given the treeFrame and sorted cp values
  # return the splits that have been removed for each cp
  res <- list()
  for (i in 1:length(cps)){
    inx.rm <- which(treeFrame$cp <= cps[i])
    res[[i]] <- treeFrame$pnode[inx.rm]
    if (length(inx.rm) > 0) {
      treeFrame <- treeFrame[- inx.rm, ]
    }
  }
  res[[i+1]] <- treeFrame$pnode
  res
}
