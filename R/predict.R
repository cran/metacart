#' Predictions from a fitted metacart object
#'
#' Returns a data frame of predicted effect sizes and moderators from a fitted metacart object
#'
#' @param object fitted model object of class "FEmrt".
#' @param newdata data frame containing the values at which predictions are required.
#' @param ... Arguments that pass to other methods.
#' @return  A data frame containing the predicted effect size, the moderators, and the corresponding node lables in the fitted tree.
#' @importFrom stats as.formula delete.response model.frame model.response terms predict
#' @export

predict.FEmrt <- function(object, newdata, ...){
  if (!inherits(object, "FEmrt"))
    warning("calling predict.FEmrt(<fake-FEmrt-object>) ...")
  if (length(object$n) < 2) {
    warning("No tree was detected, all effect sizes are predicted as overall effect size")
    data.frame(pred.y = rep(object$g, nrow(newdata)))
  } else {
    mf <- as.formula(object$call)
    tt <- terms(mf)
    ms <- model.frame(delete.response(tt), newdata)
    tree <- object$tree
    pred.efk <- predict(tree, newdata, type = "vector", ...)
    inx <- match(pred.efk, predict(tree, type="vector"))
    pred.node <- tree$where[inx]
    data.frame(pred.y = pred.efk, ms, node = pred.node)
  }

}

#' Predictions from a fitted metacart object
#'
#' Returns a data frame of predicted effect sizes and moderators from a fitted metacart object
#'
#' @param object fitted model object of class "REmrt".
#' @param newdata data frame containing the values at which predictions are required.
#' @param ... Arguments that pass to other methods.
#' @return  A data frame containing the predicted effect size, the moderators, and the corresponding node lables in the fitted tree.
#' @importFrom stats as.formula delete.response model.frame model.response terms predict
#' @export
predict.REmrt <- function(object, newdata, ...){
  tt <- terms(object$data)
  ms <- model.frame(delete.response(tt), newdata)
  oms <- model.frame(delete.response(tt), object$data)
  tree <- object[["tree"]]
   if (any(sapply(ms, class) != sapply(oms, class)))
     stop("The type of the variables do not match")
  if(is.null(tree)) {
    pred.node <- rep(1, nrow(ms))
    pred.y <- object[["g"]]
  } else {
    tnode <- rep(1, nrow(ms))
    nodes <- tnode
    for (i in 1:(nrow(tree) - 1)){
      tinx <- which(tnode == tree[i+1, "pleaf"])
      tempm <- ms[tree[i+1, "mod"]]
      if(sapply(tempm, is.numeric) == TRUE) {
        tnode[tinx] <- ifelse(tempm[tinx,1] <= object[["cpt"]][[i]], 2*i, 2*i+1)
      } else {
        tnode[tinx] <- ifelse(tempm[tinx,1] %in% oms[,tree[i+1, "mod"]],
                              ifelse(tempm[tinx,1] %in% object[["cpt"]][[i]], 2*i, 2*i+1),
                              NA)
      }
      nodes <- cbind(nodes, tnode)
    }
    pred.node <- nodes[ ,ncol(nodes)]
    pred.y <- object[["g"]][as.character(pred.node)]
  }
  cbind(newdata, pred.node, pred.y)
}
