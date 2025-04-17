fe_predict_node <- function(treeFrame, cpt, currentNode, newdata, olddata){
  while (any(unique(currentNode) %in% treeFrame$pnode)) {
    for (i in intersect(unique(currentNode), treeFrame$pnode)) {
      inx <- currentNode == i 
      inx[is.na(inx)] <- FALSE # NA memberships are kept unchanged
      inx.split <- which(treeFrame$pnode == i)
      varname <- treeFrame$var[inx.split]
      xi <- newdata[ ,as.character(varname)][inx]
      if (is.numeric(xi)) {
        currentNode[inx] <- ifelse(xi < cpt[[inx.split]], i*2, i*2 + 1)
      } else {
        currentNode[inx] <- ifelse(xi %in% olddata[ ,varname], ifelse(xi %in% cpt[[inx.split]], i*2, i*2 + 1) , NA)
      }
    }
  }
  currentNode
}


fe_predictY<- function(treeFrame, cpt, prunePath, newdata, olddata,
                               allLeaves, allWmeans){
  res <- matrix(ncol = length(prunePath), nrow = nrow(newdata))
  currentNode <- rep(1, nrow(newdata))
  for (j in length(prunePath):1) {
    inx.temp <- match(prunePath[[j]], treeFrame$pnode)
    if (length(inx.temp) > 0) {
      currentFrame <- treeFrame[inx.temp, ]
      currentCPT <- cpt[inx.temp]
      currentNode <- fe_predict_node(currentFrame, currentCPT, currentNode, newdata, olddata)
    }
    currentY <- allWmeans[match(currentNode, allLeaves)]
    res[,length(prunePath)+1-j] <- currentY
  }
  res
}
