#' connectivity_adj_mat
#'
#' @description
#' target function reveal the hierarchy structure
#'
#' @param adj_mat adjacency matrix of the target directed and weighted network
#' @param self_loop whether the NULL model consider self loop
#' @return target function matrix like modularity matrix
#' @export
#'
connectivity_adj_mat=function(adj_mat,self_loop=T)
{
  if(self_loop)
  {
  N=ncol(adj_mat)
  out_degree=apply(adj_mat,1,sum)
  in_degree=apply(adj_mat,2,sum)
  adj_connectivity=adj_mat/out_degree-t((in_degree%x%t(rep(1,N)))/sum(in_degree))
  #无自环网络
  #diag(adj_connectivity)=0
  return(adj_connectivity/N)
  }
  else{
    N=ncol(adj_mat)
    out_degree=apply(adj_mat,1,sum)
    in_degree=apply(adj_mat,2,sum)
    adj_connectivity=adj_mat/out_degree-t((in_degree%x%t(rep(1,N)))/(sum(in_degree)-in_degree))
    #无自环网络
    diag(adj_connectivity)=0
    return(adj_connectivity/N)
  }
}
