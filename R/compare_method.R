
#' compare_method
#'
#' @description
#' compare community detection algorithm with the given partition of nodes using "infomap","fastgreedy","walktrap","louvain" and "leidun"
#' The comparation is performed on undirected and weighted network
#'
#' @param adj_mat adjacency matrix of the target undirected and weighted network
#' @param truth partition of the node compared with
#' @return NMI and ARI of the comparation result with "infomap","fastgreedy","walktrap","louvain" and "leidun"
#' @export
#'
compare_method=function(adj_mat,truth)
{
  G=graph.adjacency(adj_mat,mode = "undirected")
  cd_infomap=cluster_infomap(G)
  cd_leiden=cluster_leiden(G,objective_function = "modularity")
  cd_fastgreedy=cluster_fast_greedy(G)
  cd_louvain=cluster_louvain(G)
  cd_walktrap=cluster_walktrap(G)
  #cd_optimal=cluster_optimal(G)
  info1=compare(truth,cd_infomap,method = "adjusted.rand")
  info2=compare(truth,cd_infomap,method = "nmi")
  fast1=compare(truth,cd_fastgreedy,method = "adjusted.rand")
  fast2=compare(truth,cd_fastgreedy,method = "nmi")
  wt1=compare(truth,cd_walktrap,method = "adjusted.rand")
  wt2=compare(truth,cd_walktrap,method = "nmi")
  lv1=compare(truth,cd_louvain,method = "adjusted.rand")
  lv2=compare(truth,cd_louvain,method = "nmi")
  ld1=compare(truth,cd_leiden,method = "adjusted.rand")
  ld2=compare(truth,cd_leiden,method = "nmi")
  #opt1=compare(truth,cd_optimal,method = "adjusted.rand")
  #opt2=compare(truth,cd_optimal,method = "nmi")
  results=matrix(c(info1,info2,fast1,fast2,wt1,wt2,lv1,lv2,ld1,ld2),nrow = 2)
  rownames(results)=c("ARI","NMI")
  colnames(results)=c("infomap","fastgreedy","walktrap","louvain","leidun")
  return(results)
}
