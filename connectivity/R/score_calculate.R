
#' score_calculate
#'
#' @description
#' calculate the target function based on score matrix and partition
#'
#' @param adj_rely score function
#' @param partition node partition list
#' @return score calculated according to the partition and
#' @export
#'
score_calculate=function(adj_rely,partition)
{

  community_scale=sapply(partition,length)
  temp_community_df=cbind(unlist(partition),rep(1:length(partition),community_scale))
  partition=temp_community_df[order(temp_community_df[,1]),][,2]
  #
  dummy_mat=nnet::class.ind(partition)
  score_calculate=sum(diag(t(dummy_mat)%*%adj_rely%*%dummy_mat))
  return(score_calculate)
}
