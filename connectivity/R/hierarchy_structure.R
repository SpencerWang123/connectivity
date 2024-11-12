#' hierarchy_community
#'
#' @description
#' hierarchy bottom-up algorithm using target function
#'
#' @param adj_rely targrt function matrix
#' @return combine seq and combine score
#' @export
#'
hierarchay_community=function(adj_rely)
{
  combine_record=NULL
  combine_score_all=NULL
  #测试##
  #adj_rely=connectivity_adj_mat(adj_test)
  colnames(adj_rely)=1:ncol(adj_rely);rownames(adj_rely)=1:nrow(adj_rely)
  #######
  while(nrow(adj_rely)>1)
  {
    adj_rely_judge=adj_rely;diag(adj_rely_judge)=-999
    combine_index=which(adj_rely_judge==max(adj_rely_judge),arr.ind = T)
    #合并行列-根据矩阵含义，i行j列表是的i对于j的依赖，所以我们吧i融入j中，然后我们更换行列名

    for(i in 1)#nrow(combine_index))
    {
      index_1=rownames(adj_rely)[combine_index[i,1]];index_2=rownames(adj_rely)[combine_index[i,2]]
      #更新行列
      adj_rely[index_2,]=adj_rely[index_1,]+adj_rely[index_2,];adj_rely[,index_2]=adj_rely[,index_1]+adj_rely[,index_2]
      K=nrow(adj_rely)
      adj_rely=adj_rely[-which(rownames(adj_rely)==index_1),];
      #处理最后一部分为两个社团的时候的情形
      if(K==2)
      {
        adj_rely=t(adj_rely);adj_rely=adj_rely[,-which(colnames(adj_rely)==index_1)]
        combine_record=rbind(combine_record,c(index_1,index_2));combine_score_all=c(combine_score_all,adj_rely)
        adj_rely=matrix(adj_rely)
      }
      else
      {
        adj_rely=adj_rely[,-which(colnames(adj_rely)==index_1)]
        combine_record=rbind(combine_record,c(index_1,index_2));combine_score_all=c(combine_score_all,sum(diag(adj_rely)))
      }

    }
  }
  combine_seq=cbind(combine_record,combine_score_all)
  colnames(combine_seq)=c("from","to","score")
  return(combine_seq)
}

