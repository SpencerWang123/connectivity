
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


connectivity_projection2d=function(adj_mat,tuning_in=0.5,tuning_out=0.5,alpha=0.3,seeds=10,memory=60,cut=0.5,dim_4=FALSE,fisher=F)
{
  N=nrow(adj_mat)
   #two seqs denote as "a & b"
  a=connectivity_seq(adj_mat=adj_mat,memory = memory,sample_memory = N,length_seq = N,seeds = seeds,tuning_in=tuning_in,tuning_out=tuning_out)
  b=connectivity_seq(adj_mat=adj_mat,memory = memory,sample_memory = N,length_seq = N,seeds = seeds+1,tuning_in=tuning_in,tuning_out=tuning_out)
  aa=(1:N)[order(a$node_seq)]
  bb=(1:N)[order(b$node_seq)]
  coordinate=cbind(aa,bb)
  #
  N=length(aa)
  cut_index=seq(0,N,N*cut)
  x_level=cut(aa,cut_index)
  y_level=cut(bb,cut_index)
  test_chisq=chisq.test(x_level,y_level)
  test_fisher=NULL
  if(fisher)
  {
  test_fisher=fisher.test(x_level,y_level)$p.value
  }
  #test permute
  aP=connectivity_seq(adj_mat=adj_mat,memory = memory,sample_memory = N,length_seq = N,seeds = seeds+3,tuning_in=tuning_in,tuning_out=tuning_out)
  bP=connectivity_seq(adj_mat=adj_mat,memory = memory,sample_memory = N,length_seq = N,seeds = seeds+4,tuning_in=tuning_in,tuning_out=tuning_out)
  aa_p=(1:N)[order(aP$node_seq)]
  bb_p=(1:N)[order(bP$node_seq)]
  x_level_p=cut(aa_p,cut_index)
  y_level_p=cut(bb_p,cut_index)
  x_level_pp=as.factor(paste(x_level,x_level_p))
  y_level_pp=as.factor(paste(y_level,y_level_p))
  test_chisq_p3=chisq.test(x_level_pp,y_level)
  test_chisq_p4=chisq.test(x_level_pp,y_level_pp)

##
  all_result=c(test_chisq$p.value,test_fisher,test_chisq_p3$p.value,test_chisq_p4$p.value)
  result=list(connectivity_1=a,connectivity_2=b,coordinate=coordinate,p.value_all=all_result,test_fisher=test_fisher
              ,test_chisq=test_chisq,test_chisq_p3=test_chisq_p3,test_chisq_p4=test_chisq_p4)
  return(result)
}
