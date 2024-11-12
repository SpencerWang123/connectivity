dist_label=function(mat)
{
  dist_mat=matrix(rep(0,nrow(mat)^2),ncol=nrow(mat))
  N=ncol(mat)
  for (i in 1:nrow(mat))
  {
    temp_row=mat[i,]
    logic_table=(temp_row==t(mat))
    dist_mat[,i]=apply(logic_table,2,sum)
  }
  dist_mat=N-dist_mat
  rownames(dist_mat)=c(1:nrow(mat))
  colnames(dist_mat)=rownames(dist_mat)
  return(dist_mat)
}
