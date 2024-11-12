#' cd_based_on_score
#'
#' @description
#' iteratively divide the network into 2 section--communities, until the target function won't increase. Need function connectivity_adj_mat and cut_into_2 in advance
#'
#' @param adj_mat adjacency matrix of the target undirected and weighted network
#' @return partition list and community labels
#' @export
#'
cd_based_on_score=function(adj_mat,seeds=1)
{
  N=nrow(adj_mat)
  adj_rely=connectivity_adj_mat(adj_mat)
  #
  adj_mat_cutable=list(adj_mat)
  adj_list_uncutable=list()
  #
  node_seq_cutable=list(1:N)
  node_seq_uncutable=list()
  temp_cut=1;temp_score=0
  temp_seq_overview=connectivity_seq(adj_mat,memory = 60,length_seq = N,sample_memory = N,seeds = seeds)

  #adj_list不断更新，不能分的合能分的分开放
  while (length(adj_mat_cutable)>0)
  {
    adj_mat_cutable_beifen=list()
    node_seq_cutable_beifen=list()
    adj_list_uncutable_beifen=adj_list_uncutable
    node_seq_uncutable_beifen=node_seq_uncutable

    for (i in 1:length(adj_mat_cutable))
    {
      temp_cut_seq=node_seq_cutable[i][[1]]
      temp_adj_mat=adj_mat_cutable[i][[1]]
      temp_cut=cut_into_2(adj_mat,temp_adj_mat,target_node =temp_cut_seq,temp_seq_overview = temp_seq_overview)
      #############################
      #泽丽默认是那个176节点的算法了
      #############################
      #分割完了合并分割结果判断是否有提升
      node_seq_1=temp_cut$node_seq_1_refine
      node_seq_2=temp_cut$node_seq_2_refine
      temp_partition=c(node_seq_uncutable,node_seq_cutable[-i],list(node_seq_1),list(node_seq_2))
      new_score=score_calculate(adj_rely,temp_partition)
      if(new_score>temp_score & min(length(node_seq_1),length(node_seq_2))>2)#如果可分割，那么从可分割list移除可分割的社团，然后替换成其分割的两个新的
        #不需要用分数来判断，前面已经帮你判断了！只需要看node-seq-2-refine长度是否为零就可以
      {
        adj_mat_cutable_beifen=c(adj_mat_cutable_beifen,list(temp_cut$adj_mat_1),list(temp_cut$adj_mat_2))
        node_seq_cutable_beifen=c(node_seq_cutable_beifen,list(node_seq_1),list(node_seq_2))
        #print(new_score)
      }
      else{
        adj_list_uncutable_beifen=c(adj_list_uncutable_beifen,list(temp_adj_mat))
        node_seq_uncutable_beifen=c(node_seq_uncutable_beifen,list(temp_cut_seq))
      }
      #print(community_lfr500[node_seq_1]);
      #print(community_lfr500[node_seq_2])
    }
    adj_mat_cutable=adj_mat_cutable_beifen
    adj_list_uncutable=adj_list_uncutable_beifen
    node_seq_cutable=node_seq_cutable_beifen
    node_seq_uncutable=node_seq_uncutable_beifen
    temp_score=score_calculate(adj_rely,c(node_seq_cutable,node_seq_uncutable))

  }
  #梳理社团
  community_scale=sapply(node_seq_uncutable,length)
  temp_community_df=cbind(unlist(node_seq_uncutable),rep(1:length(node_seq_uncutable),community_scale))
  community_df=temp_community_df[order(temp_community_df[,1]),][,2]
  return(c(partition_list=list(node_seq_uncutable),community_df=list(community_df)))
}
