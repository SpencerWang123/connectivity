#' connectivity_test
#' @description
#' function to generate node seq and the corresponding connectivity seq on the given network
#'
#' @param adjmat adjacency matrix of a network which should have identical colnames and rownames
#' @param tuning_in weight of the in-connect when evaluating the connectivity. float in [0,1]
#' @param tuning_out weight of the out-connect when evaluating the connectivity. float in [0,1]
#' @param alpha smooth parameter which reflects how much the connectivity rely on the temp result. float in [0,1]
#' @param seed random seed used to initialize the first node in Node seq
#' @param memory memory length, which decides how many recent nodes are taken into consideration when calculating the connectivity. int
#' @param sample_memory denote how recent of those nodes in node seq won't be the candidate of the new node. int
#' @param jump possibility of choosing a random node regardless of the connectivity when generating a new node in node seq. float in [0,1]
#' @param length_seq length of the connectivity seq and node seq
#'
#' @return a list contenting node seq and connectivity seq
#' @export
#'
connectivity_test=function(adj_mat,tuning_in=0.5,tuning_out=0.5,alpha=1,seeds=10,memory,sample_memory,jump=0,length_seq)#
{
  set.seed(seeds)
  N=ncol(adj_mat);out_degree=apply(adj_mat,1,sum);in_degree=apply(adj_mat,2,sum)
  #空表
  connect_record=0;
  temp_smooth=0;connect_record_smooth=0
  connect_maximun_point=0
  connect_record_memory=0;
  temp_node_list=NULL
  dip_seq=rep(1,memory)
  dip_stat_seq=rep(1,memory)
  #孤立点问题
  isolated_node=which(out_degree+in_degree==0);
  occupied_node=c(isolated_node)
  #初始化点
  temp_node=sample(c(1:N)[!c(1:N)%in%isolated_node],1);

  #循环开始生成社群
  for (i in 1:(length_seq-1))
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
    index_beixuan=which(!names(temp_connect)%in%tail(temp_node,sample_memory))
    max_connect=max(temp_connect[index_beixuan],na.rm = T)
    if(i>=memory)
    {
      temp_connectivity=temp_connect#[index_beixuan]
      temp_dip=diptest::dip.test(temp_connectivity)
      dip_seq=c(dip_seq,temp_dip$p.value)
      dip_stat_seq=c(dip_stat_seq,temp_dip$statistic)

    }
    high_rank_all=index_beixuan[index_beixuan%in%names(which.max(temp_connect[index_beixuan]))]
    high_rank=high_rank_all[sample(length(high_rank_all),1)]
    #从这里设置jump的步骤-截断temp_node，塑造新的temp_node
    if(rbinom(1,1,jump)>0)
    {
      temp_node_list=c(temp_node_list,temp_node)
      temp_node=sample(c(1:N)[!c(1:N)%in%isolated_node],1)
      connect_record_smooth=c(connect_record_smooth,0)
    }
    else
    {
      temp_smooth=(1-alpha)*temp_smooth+alpha*max_connect
      #根据记忆，记录最大连接度节点
      temp_node=c(temp_node,high_rank);
      #记录连接度（原始和未平滑）
      connect_record_smooth=c(connect_record_smooth,temp_smooth)
    }
  }
  temp_node_list=c(temp_node_list,temp_node)
  node_seq_all=unlist(temp_node_list)
  dip_significant=sum((dip_seq[(memory+1):length(dip_seq)]<=0.05))
  results=c(list(node_seq_all),list(connect_record_smooth),list(dip_seq),list(dip_significant),list(dip_stat_seq))

  names(results)=c("node_seq","connect_record_smooth","dip_test","dip_significant","dip_statistic")
  return(results)
}

