
#' plot_data_function
#'
#' @description
#' generating data for visualization according to the given partition of the nodes in the connectivity analysis
#'
#' @param temp_seq node seq needed
#' @param temp_connect_record connectivity seq needed
#' @param temp_community_df data.frame composed of node labels and corresponding community labels
#' @return data.frame composed of 4 columns which contents node, connect_record, community and color
#' @export
#'
plot_data_function=function(temp_seq,temp_connect_record,temp_community_df,seed=1,col_given=NULL)
{
  #下面按照真实颜色画图！
  library(RColorBrewer)
  #可视化真实社团
  colnames(temp_community_df)=c("node","community")
  temp_community=temp_community_df[,2]
  community_num=length(unique(temp_community))
  set.seed(seed)
  if(length(col_given)==0){col_list=rainbow(community_num)
  col_available=cbind(sample(1:community_num,community_num),col_list);colnames(col_available)=c("community","color")
  }
  else
  {
  col_list=col_given
  col_available=cbind(1:community_num,col_list);colnames(col_available)=c("community","color")
  }
  true_community=as.data.frame(merge(temp_community_df,col_available))
  plot_data=as.data.frame(cbind(temp_seq,temp_connect_record));colnames(plot_data)=c("node","connect_record")
  plot_data_all=plyr::join(plot_data,true_community)
  #对节点进行重新排序按照原本顺序
  return(plot_data_all)
}
