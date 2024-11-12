random_adj_mat=function(adj_mat)
{
  #adj_mat=adj_mat_sbm120
  library(igraph)
  out_deg=apply(adj_mat,1,sum)
  in_deg=apply(adj_mat,2,sum)
  RG=sample_degseq(out.deg = out_deg,in.deg = in_deg)
  RG_adjmat=as.matrix(as_adjacency_matrix(RG))
  colnames(RG_adjmat)=colnames(adj_mat)
  rownames(RG_adjmat)=rownames(adj_mat)
  return(RG_adjmat)
}
