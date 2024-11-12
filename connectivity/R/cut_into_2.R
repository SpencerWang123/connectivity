#' cut_into_2
#'
#' @description
#' cut the network into 2 sub-graph while taking the refine process into consideration.
#'
#' @param adj_mat_all adjacency matrix of the target directed and weighted network
#' @param adj_mat local adjacency matrix of the target directed and weighted network
#' @param target_node node list
#' @param min_scale min scale to be cut
#' @param temp_seq_overview results returned by the connectivity_seq function
#' @return two node lists and the corresponding adjacency matrix
#' @export
#'
cut_into_2=function(adj_mat_all,adj_mat,target_node,min_scale=3,temp_seq_overview)
{
  k=min_scale
  N=nrow(adj_mat)
  Nall=nrow(adj_mat_all)
  #尝试着全局切割
  #temp_seq_overview=connectivity_seq(adj_mat_all,length_seq = Nall,memory = min(Nall,60),sample_memory = Nall)
  temp_node_seq=temp_seq_overview$node_seq[temp_seq_overview$node_seq%in%target_node]
  temp_connectivity=temp_seq_overview$connect_record_smooth[temp_seq_overview$node_seq%in%target_node]
  temp_seq=list(node_seq=temp_node_seq,connect_record_smooth=temp_connectivity)
  #
  adj_rely_all=connectivity_adj_mat(adj_mat_all)##这里我们计算adjrely的时候为了保存全局信息，用的是全局的adjrely的子矩阵
  adj_rely=adj_rely_all[target_node,target_node]
  adj_rely_local=connectivity_adj_mat(adj_mat)
  ############################开始根据连接度序列二分并优化二分结果################################
  #从此开始refine-yinweirefine会调整顺序
  #重排矩阵的时候出了问题，识别的位置是数字标识的，这里给的要求是节点编号，要转化
  index_1_order=data.frame(node=as.numeric(colnames(adj_mat)),rank=1:N);index_2_order=data.frame(node=as.numeric(temp_seq$node_seq),rank2=1:N)
  merge_order=merge(index_2_order,index_1_order,by="node")
  reorder_rank=merge_order$rank[order(merge_order$rank2)]
  adj_rely_adj=adj_rely[reorder_rank,reorder_rank]
  adj_rely_adj_local=adj_rely_local[reorder_rank,reorder_rank]


  index_mat_l=rep(1,N)%*%t(rep(1,N));index_mat_l[lower.tri(index_mat_l)]=0
  index_mat_u=rep(1,N)%*%t(rep(1,N));index_mat_u[upper.tri(index_mat_u)]=0
  cut_back=c(diag(index_mat_l%*%adj_rely_adj%*%index_mat_u),0)#右侧全部，右侧N-1,右侧N-2
  cut_forth=c(0,diag(index_mat_u%*%adj_rely_adj%*%index_mat_l))#左侧全无，一个，2个
  cut_all=(cut_back+cut_forth)[(1+k):(N+1-k)];
  #切的时候用local来切，算的时候用全局来算
  cut_back_local=c(diag(index_mat_l%*%adj_rely_adj_local%*%index_mat_u),0)
  cut_forth_local=c(0,diag(index_mat_u%*%adj_rely_adj_local%*%index_mat_l))
  cut_all_local=cut_back_local+cut_forth_local;
  #全局不切分计算
  score_uncut=cut_all[1]#选择第k个比较，我们对最小社团规模有要求
  max_score=cut_all[which.max(cut_all)]
  if(score_uncut>=max_score)
  {return(list(node_seq_1_refine=target_node,node_seq_2_refine=NULL,
               adj_mat_1=adj_mat,
               adj_mat_2=NULL,best_score=score_uncut))}


  else{

    temp_score=0
    node_seq_1=temp_seq$node_seq[1:(which.max(cut_all)+k)];N1=length(node_seq_1)
    node_seq_2=temp_seq$node_seq[(which.max(cut_all)+k+1):N];N2=length(node_seq_2)
    node_seq_refine=temp_seq$node_seq
    node_seq_1_refine=node_seq_1
    node_seq_2_refine=node_seq_2
    #refine_process
    switch_index_1=c(rep(1,N1),rep(0,N2));switch_index_2=rep(1,N)-switch_index_1
    #社团1换位
    switch_num=1
    while(switch_num>0 & max_score>temp_score & min(c(length(node_seq_1_refine),length(node_seq_2_refine)))>1)
    {
      temp_score=max_score
      mat_0=matrix(rep(0,N1*N),ncol=N1);mat_1=matrix(rep(1,N1*N),ncol=N1)
      switch_mat11=mat_0;switch_mat11[1:N1,]=1;diag(switch_mat11[1:N1,])=0;switch_mat12=mat_1-switch_mat11
      switch_result_1=diag(t(switch_mat11)%*%adj_rely_adj%*%switch_mat11)+diag(t(switch_mat12)%*%adj_rely_adj%*%switch_mat12)-c(temp_score)
      #print(switch_result_1[which(switch_result_1>0)])
      #社团2换位
      mat_0=matrix(rep(0,N2*N),ncol=N2);mat_1=matrix(rep(1,N2*N),ncol=N2)
      switch_mat22=mat_0;switch_mat22[(N1+1):N,]=1;diag(switch_mat22[(N1+1):N,])=0;switch_mat21=mat_1-switch_mat22
      switch_result_2=diag(t(switch_mat21)%*%adj_rely_adj%*%switch_mat21)+diag(t(switch_mat22)%*%adj_rely_adj%*%switch_mat22)-c(temp_score)

      switch_num=sum(c(switch_result_1>0,switch_result_2>0))
      #输出换位结果
      new_index_1=as.logical(c(1-(switch_result_1>0),(switch_result_2>0)+0));
      new_index_2=as.logical(c((switch_result_1>0)+0,1-(switch_result_2>0)))

      #更新最优切分结果
      max_score=new_index_1%*%adj_rely_adj%*%(new_index_1)+new_index_2%*%adj_rely_adj%*%(new_index_2)
      if(max_score>temp_score)#如果分数更高更新分个序列，调整矩阵顺序，用于下一次调换，反之不更新
      {
        node_seq_1_refine=node_seq_refine[new_index_1];N1=length(node_seq_1_refine)
        node_seq_2_refine=node_seq_refine[new_index_2];N2=length(node_seq_2_refine)
        node_seq_refine=c(node_seq_1_refine,node_seq_2_refine)
        print(max_score)
        #调整矩阵顺序
        index_1_order=data.frame(node=as.numeric(colnames(adj_rely_adj)),rank=1:N);index_2_order=data.frame(node=as.numeric(node_seq_refine),rank2=1:N)
        merge_order=merge(index_2_order,index_1_order,by="node")
        reorder_rank=merge_order$rank[order(merge_order$rank2)]
        adj_rely_adj=adj_rely_adj[reorder_rank,reorder_rank]
        #
        print(switch_num)
      }

    }


    return(list(node_seq_1_refine=node_seq_1_refine,node_seq_2_refine=node_seq_2_refine,
                adj_mat_1=adj_mat_all[node_seq_1_refine,node_seq_1_refine],
                adj_mat_2=adj_mat_all[node_seq_2_refine,node_seq_2_refine],best_score=max_score
    ))
  }

}

