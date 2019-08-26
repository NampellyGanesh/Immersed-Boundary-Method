  module variables
  implicit none

  integer                              ::     imax,jmax,k,i,j,no_of_cells,istop,restart,viscous
  integer                              ::     h=1,lbc,rbc,tbc,bbc,irest,iprint,ipurge,order,immersed
  real                                 ::     R=287.0,gama=1.4,kappa,mue
  real                                 ::     uinf,vinf,T_inf,pinf,rho_inf,CFL
  character(len=40)                    ::     file_name,filename_results
  
  !IBM variables
  integer                              ::     n_ibpts,closed_profile,imin_fc,imax_fc,jmin_fc,jmax_fc


  integer,dimension(:,:), allocatable  ::     cc,c_ib,bc,fc

  end module variables
  
