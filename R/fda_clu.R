fda_clu=function(connect_record_smooth,knot_gap,basis_order,miniest_cluster_scale)
{
  M=length(connect_record_smooth)
  x=1:M
  y=connect_record_smooth
  knots    = c(seq(1,M,knot_gap)) #Location of knots
  n_knots   = length(knots) #Number of knots
  n_order   = basis_order # order of basis functions: for cubic b-splines: order = 3 + 1
  n_basis   = length(knots) + n_order - 2;
  basis = fda::create.bspline.basis(rangeval = c(1,M), n_basis)
  fd=fda::Data2fd(argvals = x,y=y,basisobj = basis)
  D_fd=fda::deriv.fd(fd)
  eval_x=seq(1,M,0.01)
  eval_Dfd=as.vector(fda::eval.fd(D_fd,seq(1,M,0.01)))
  #求零点
  sign_Dfd_0=sign(sign(eval_Dfd)[1:(length(eval_Dfd)-1)]+sign(eval_Dfd)[2:length(eval_Dfd)])
  change_sign=which(sign_Dfd_0==0)
  cut_points_temp=c(floor(change_sign[which(sign_Dfd_0[change_sign-1]<0)]*0.01),M)
  #截断点之间差分小于最小社团规模剔除
  remain_cut=which(diff(c(1,cut_points_temp))>miniest_cluster_scale)
  cut_points_temp= cut_points_temp[remain_cut]
  return(cut_points_temp)
}
