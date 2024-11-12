#' parallel community affiliation test for nodes
#'
#' @param adj_mat adjacency matrix of the graph which should be featured with colnames as well as rownames
#' @param N default to be number of nodes
#' @param aplha  first order sliding smoothing parameter of the depedency sequence, default to be 0.7
#' @param num_cores parallel computation cores used,default to be 8
#' @param m number of sequences sampled from the network
#' @param memory sample memory, default to be 50
#' @param tuning_in importance of in-edges, default to be 1 which only make differences when the graph is directed
#' @param tuning_out importance of out-edges, default to be 1 which only make differences when the graph is directed
#'
#' @return test_stat:test statistics,p_value:test p values,test_record:all sampling results

test_connectivity=function(adj_mat,N=NULL,aplha=0.7,num_cores=8,m=100,memory=50,tuning_in=1,tuning_out=0)
{
  if(is.null(N)){N=nrow(adj_mat)}
  test1=wavelet_test_long_parallel(adj_mat,memory = memory,length_seq = N,sample_memory = N,m=m,tuning_in=tuning_in,tuning_out=tuning_out,alpha = alpha,num_cores = num_cores)
  test10=wavelet_test_long_rg_parallel(adj_mat,memory = memory,length_seq = N,sample_memory = N,m=m,tuning_in=tuning_in,tuning_out=tuning_out,alpha = alpha,num_cores = num_cores)
  #对齐两次模拟结果的节点，防止有的节点在某些采样中从未出现
  mutual_index1=colnames(test1$local_sum_mean) %in% colnames(test10$local_sum_mean)
  mutual_index0=colnames(test10$local_sum_mean) %in% colnames(test1$local_sum_mean)

  mean_stat=(test1$local_sum_mean-apply(test1$local_sum_mean,1,mean))[,mutual_index1]-(test10$local_sum_mean-apply(test10$local_sum_mean,1,mean))[,mutual_index0]
  var_stat=((((N-2)/N)^2*test1$local_sum_sd^2)/test1$n_seq)[,mutual_index1]+(1/N^2)*apply(test1$local_sum_sd^2/test1$n_seq,1,sum)+((((N-2)/N)^2*test10$local_sum_sd^2)/test10$n_seq)[,mutual_index0]+(1/N^2)*apply(test10$local_sum_sd^2/test10$n_seq,1,sum)
  test_stat=mean_stat/sqrt(var_stat)
  test_p=1-pnorm(test_stat)
  return(list(test_stat=test_stat,p_value=test_p,test_record=list(test1=test1,test10=test10)))
}
