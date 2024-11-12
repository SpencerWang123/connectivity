
#' gererate_G_BR
#'
#' @description
#' generate igraph object and adjacency matrix from BR_data
#'
#' @param trade br_y or br_q
#' @param time_wanted time needed
#' @return igraph object G and adjacency matrix
#' @export
#'
gererate_G_BR=function(trade,time_wanted)
{
  g=igraph::make_graph(as.character(matrix(t(trade[,c(2,6)]))))
  E(g)$weight=trade[,grep(time_wanted,colnames(trade))]
  mat=as.matrix(as_adjacency_matrix(g,attr = "weight"))
  return(list(g,mat))
}
