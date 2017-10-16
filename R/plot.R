#' Visualisation of a FE meta-tree
#'
#' Plot function for a \code{FEmrt} object. The plot shows the result of \code{FEmrt}.
#' The plot function uses the plot method from the package \pkg{rpart.plot} of Stephen Milborrow
#' (2016).
#'
#' For categorical variables we recommend to use short names for levels to avoid overlapping
#' labels at split points.
#' For each split of the tree it holds that if the condition (in the parent node) is satisfied,
#' then the study goes to the left node and if not it goes to the right node.
#' @method plot FEmrt
#' @param x A FEmrt object.
#' @param ... Arguments that pass to prp()
#' @importFrom rpart.plot prp
#' @export
plot.FEmrt <- function(x, ...){
  if (length(x$n) < 2) {stop("no tree was detected")}
  else {prp(x$tree, ...)}
}
