#' connectivity_cut_thres
#' @description
#' function to generate node seq and the corresponding connectivity seq on the given network.
#' In this algorithm the connectivity seq will be cut according to the given threshold.
#' the cut of the connectivity will return a partition of node seq accordingly which served as community detection result.
#'
#' @param adjmat adjacency matrix of a network which should have identical colnames and rownames
#' @param tuning_in weight of the in-connect when evaluating the connectivity. float in [0,1]
#' @param tuning_out weight of the out-connect when evaluating the connectivity. float in [0,1]
#' @param alpha smooth parameter which reflects how much the connectivity rely on the temp result. float in [0,1]
#' @param seed random seed used to initialize the first node in Node seq
#' @param memory memory length, which decides how many recent nodes are taken into consideration when calculating the connectivity. int

#'
#' @return
#' membership: node_seq and the community label of each node in order of the connectivity seq
#'
#' connect_record: connectivity seq
#'
#' cluster: nodes and their community label in order of the truth community partition
#'
#' @export
#'
connectivity_cut_thres=function(adj_mat,tuning_in=0.5,tuning_out=0.5,threshold=0.75,alpha=1,seeds=1,memory)##更新了平滑措施
{
  #
  N=ncol(adj_mat)
  set.seed(seeds)
  #空表
  connect_record=0;
  membership=list()
  connect_record_smooth=0##平滑记录空表
  connect_maximun_point=0##当前社团最大连接度-来自平滑结果
  temp_smooth=0
  #网络信息统计
  out_degree=apply(adj_mat,1 ,sum )
  in_degree=apply(adj_mat,2 ,sum )
  #孤立点问题
  isolated_node=which(out_degree+in_degree==0)
  occupied_node=c(isolated_node)
  temp_node=sample(c(1:N)[!c(1:N)%in%isolated_node],1)
  occupied_node=c(occupied_node,temp_node)
  for (i in 1:(N-1-length(isolated_node)))
  {
    ##出度比例
    memory_out_degree=matrix(adj_mat[,colnames(adj_mat)%in%tail(temp_node,memory)],nrow = N)
    temp_connect_out_ratio=apply(memory_out_degree,1 ,sum )/(out_degree+1)
    ##入度比例
    memory_in_degree=matrix(t(adj_mat[rownames(adj_mat)%in%tail(temp_node,memory),]),nrow = N)
    temp_connect_in_ratio=apply(memory_in_degree,1 ,sum )/(in_degree+1)
    ##备选节点连接程度
    temp_connect=(tuning_in*temp_connect_in_ratio+tuning_out*temp_connect_out_ratio);names(temp_connect)=c(1:N)
    #提取现有节点之外最大的（多个满足条件数值取其中之一）
    index_beixuan=which(!names(temp_connect)%in%occupied_node)
    max_connect=max(temp_connect[index_beixuan],na.rm = T)
    high_rank=index_beixuan[index_beixuan%in%names(which.max(temp_connect[index_beixuan]))][1]
    #计算比值是否低于阈值，低于阈值的时社团归零，然后网络排除已有社团，记录下当前membership
    temp_smooth=(1-alpha)*temp_smooth+alpha*max_connect
    thres=temp_smooth/(connect_maximun_point+0.01)#防止出现NA
    if(i>2 & thres<threshold)
    {
      connect_record_smooth=0;connect_maximun_point=0
      membership=c(membership,list(temp_node));
      if(length(which(!rownames(adj_mat)%in%occupied_node))>1){
        temp_node=sample(which(!rownames(adj_mat)%in%occupied_node),1)
      }else{temp_node=which(!rownames(adj_mat)%in%occupied_node)}
      occupied_node=c(occupied_node,temp_node)
      connect_record=c(connect_record,0)
      temp_smooth=0
    }
    else
    {
      temp_node=c(temp_node,high_rank);
      occupied_node=c(occupied_node,high_rank)
      connect_record=c(connect_record,max_connect)
      temp_smooth=(1-alpha)*temp_smooth+alpha*max_connect##更新平滑的record
      connect_record_smooth=c(connect_record_smooth,temp_smooth)###
      connect_maximun_point=max(connect_record_smooth)
    }
  }
  #进度条3
  #end_time <- Sys.time();run_time <- end_time - star_time
  #print(run_time)
  #close(pb)
  #
  membership=c(membership,list(temp_node));
  #孤立点单独成为社群
  membership=c(membership,as.list(isolated_node))
  #整理类簇分类编号
  clu_scale=unlist(lapply(membership,length))
  membership_df=cbind(unlist(membership),rep(1:length(clu_scale),clu_scale))
  cluster=membership_df[order(membership_df[,1]),];colnames(cluster)=c("node_index","cluster")
  results=c(list(membership_df),list(connect_record),list(cluster))
  names(results)=c("membership","connect_record","cluster")
  return(results)
}
