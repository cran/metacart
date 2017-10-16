#' A function to update node
#'
#' @param nodes a data frame
#' @param newNode a new row
#' @param name a string
#' @return a data.frame updated the nodes
#' @keywords internal
updateNodes <- function(nodes, newNode, name = "leaf") {
  rows <- nodes[name] == newNode[1, name]
  nodes <- rbind(nodes[!rows, ], newNode)
  nodes
}

#' A function to deal with symbols
#'
#' @param input a string
#' @return converted string
#' @keywords internal
encodeHtml <- function(input) {
  # reference: https://en.wikipedia.org/wiki/List_of_XML_and_HTML_character_entity_references
  dict <- data.frame(
    c('&', '\u0026'),
    c('<=', '\u2264'),
    c('>=', '\u2265'),
    c('<', '\u003C'),
    c('>', '\u003E')
  )
  tmp <- input
  for (i in 1:ncol(dict)) {
    tmp <- gsub(dict[1, i], dict[2, i], tmp)
    Encoding(tmp) <- "UTF-8"
  }
  tmp
}

#' A function to draw an oval
#'
#' @param plotobj the obj to be plot
#' @param x x
#' @param y y
#' @param c c
#' @param x.scale x.scale
#' @param y.scale y.scale
#' @return converted string
#' @importFrom ggplot2 geom_polygon aes
#' @keywords internal
oval_draw <- function(plotobj, x, y, c, x.scale = 1, y.scale = 1, ...){
  t <- seq(-1 * pi, 1 * pi, length = 100)
  df <- data.frame(
    x = x.scale * sin(t) + x,
    y = y.scale * cos(t) / c + y
  )
  plotobj <- plotobj +
    geom_polygon(data = df,
                 aes(x, y),
                 fill = "lightgrey",
                 color = "black")
}
