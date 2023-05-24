
#' connectivity_cut_fda
#' @description
#' function to generate node seq and the corresponding connectivity seq on the given network.
#' this algorithm will fit the connectivity seq with a continuous functional object called connectivity function.
#' In this algorithm the connectivity function will be cut at where its first derivative equals to 0 while the second derivative is positive.
#' the cut of the connectivity function will return a partition of node seq accordingly which served as community detection result.
#'
#' @param adjmat adjacency matrix of a network which should have identical colnames and rownames
#' @param tuning_in weight of the in-connect when evaluating the connectivity. float in [0,1]
#' @param tuning_out weight of the out-connect when evaluating the connectivity. float in [0,1]
#' @param alpha smooth parameter which reflects how much the connectivity rely on the temp result. float in [0,1]
#' @param seed random seed used to initialize the first node in Node seq
#' @param memory memory length, which decides how many recent nodes are taken into consideration when calculating the connectivity. int
#' @param miniest_cluster_scale the minimum scale of the community could be recognized
#' @param knot_gap knot_gap used to generate spline basis in fitting the connectivity function
#' @param basis_order basis order used to generate spline basis in fittinng the connectivity function
#' @return
#' node_seq: node_seq
#'
#' cluster: nodes and their community label in order of the truth community partition
#'
#' connect_record_smooth: smoothed connectivity seq according to parameter alpha
#'
#' @export
#'
connectivity_cut_fda=function(adj_mat,tuning_in=0.5,tuning_out=0.5,alpha=0.3,seeds=10,miniest_cluster_scale=10,knot_gap=30,basis_order=4,memory=60)##更新了平滑措施
{
  set.seed(seeds)
  N=ncol(adj_mat);out_degree=apply(adj_mat,1,sum);in_degree=apply(adj_mat,2,sum)
  #空表
  connect_record=0;
  temp_smooth=0;connect_record_smooth=0
  connect_maximun_point=0
  connect_record_memory=0;
  #孤立点问题
  isolated_node=which(out_degree+in_degree==0);
  occupied_node=c(isolated_node)
  #初始化点
  temp_node=sample(c(1:N)[!c(1:N)%in%isolated_node],1);occupied_node=c(occupied_node,temp_node)

  #循环开始生成社群
  for (i in 1:(N-1-length(isolated_node)))
  {
    #出度比例
    memory_out_degree=matrix(adj_mat[,colnames(adj_mat)%in%tail(temp_node,memory)],nrow = N)
    temp_connect_out_ratio=apply(memory_out_degree,1,sum)/(out_degree+1)
    #入度比例
    memory_in_degree=matrix(t(adj_mat[rownames(adj_mat)%in%tail(temp_node,memory),]),nrow = N)
    temp_connect_in_ratio=apply(memory_in_degree,1 ,sum )/(in_degree+1)
    #汇总
    temp_connect=(tuning_in*temp_connect_in_ratio+tuning_out*temp_connect_out_ratio);names(temp_connect)=c(1:N)
    #提取现有节点之外最大的，并平滑
    index_beixuan=which(!names(temp_connect)%in%occupied_node)
    max_connect=max(temp_connect[index_beixuan],na.rm = T)
    high_rank=index_beixuan[index_beixuan%in%names(which.max(temp_connect[index_beixuan]))][1]
    temp_smooth=(1-alpha)*temp_smooth+alpha*max_connect
    #根据记忆，记录最大连接度节点
    temp_node=c(temp_node,high_rank);
    occupied_node=c(occupied_node,high_rank)
    #记录连接度（原始和未平滑）
    connect_record_smooth=c(connect_record_smooth,temp_smooth)
  }

  #函数型数据分割部分
  fda_cluster=c(0,fda_clu(connect_record_smooth,miniest_cluster_scale=10,knot_gap=30,basis_order=4))
  clu_scale_fd=c(diff(fda_cluster),rep(1,length(isolated_node)))
  #孤立点单独成为社群
  temp_node=c(temp_node,isolated_node)
  membership_df=cbind(temp_node,rep(1:length(clu_scale_fd),clu_scale_fd))
  cluster=membership_df[order(membership_df[,1]),];colnames(cluster)=c("node_index","cluster_fd")
  results=c(list(temp_node),list(cluster),list(connect_record_smooth))
  names(results)=c("node_seq","membership","connect_record_smooth")
  return(results)
}
