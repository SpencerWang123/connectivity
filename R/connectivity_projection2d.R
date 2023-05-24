
#' connectivity_projection2d
#' @description
#' function to combine 2 connect seq to construct a 2D-projection of each node which retains the community structure.
#' in addition, we perform a dip test based on the projection to test whether there are significant cluster structure in the graph
#'
#' @param adjmat adjacency matrix of a network which should have identical colnames and rownames
#' @param tuning_in weight of the in-connect when evaluating the connectivity. float in [0,1]
#' @param tuning_out weight of the out-connect when evaluating the connectivity. float in [0,1]
#' @param alpha smooth parameter which reflects how much the connectivity rely on the temp result. float in [0,1]
#' @param seed random seed used to initialize the first node in Node seq
#' @param memory memory length, which decides how many recent nodes are taken into consideration when calculating the connectivity. int
#'
#' @return
#' node_seq: 2 node sequences generated
#' connectivity_seq: 2 connectivity sequences generated
#'
#' coordinate: 2D coordinate generated
#' dip_result: result of the dip test on the 2D projection. A significant result means significant community structures
#'
#' @export
#'


connectivity_projection2d=function(adj_mat,tuning_in=0.5,tuning_out=0.5,alpha=0.3,seeds=10,memory=60)
{
  N=nrow(adj_mat)
   #two seqs denote as "a & b"
  a=connectivity_seq(adj_mat=adj_mat,memory = memory,sample_memory = N,length_seq = N,seeds = seeds)
  b=connectivity_seq(adj_mat=adj_mat,memory = memory,sample_memory = N,length_seq = N,seeds = seeds+1)
  aa=(1:N)[order(a$node_seq)]
  bb=(1:N)[order(b$node_seq)]
  library(diptest)
  coordinate=cbind(aa,bb)
  test_seq=as.vector(dist(coordinate))
  dip_result=dip.test(test_seq)
  result=list(connectivity_1=a,connectivity_1=b,coordinate=coordinate,dip_result=dip_result)
  return(result)
}
