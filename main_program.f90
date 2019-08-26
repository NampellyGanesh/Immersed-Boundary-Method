program main_program 
   
use variables
use geometrical_properties 
use immersed_routine
use flux_calculation

implicit none

real        ::    Mach,angle
real,dimension(:,:,:), allocatable   ::     domain,i_norm,j_norm,cell_c
real,dimension(:,:),   allocatable   ::     i_length,j_length,volume,u,v,p,rho,T
real,dimension(:,:),   allocatable   ::     ib_norm,ibpts
! I face in vertical and j face in horizontal direction

 open(15,file = "Inputs.dat")

      read(15,*)
      read(15,*)
      read(15,*)file_name
      read(15,*)
      read(15,*)
      read(15,*)imax,jmax,n_ibpts
      read(15,*)
      read(15,*)
      read(15,*)pinf,Mach,T_inf,lbc,rbc,tbc,bbc,istop,restart
      read(15,*)
      read(15,*)
      read(15,*)CFL,angle,irest,iprint,ipurge,order,kappa,viscous,mue,closed_profile,immersed

 close(15)

 allocate(domain(imax,jmax,2))
 allocate(i_norm(imax,2:jmax,2))
 allocate(j_norm(2:imax,jmax,2))
 allocate(cell_c(imax,jmax,2))
 allocate(i_length(imax,2:jmax))
 allocate(j_length(2:imax,jmax))

 allocate(volume(2:imax,2:jmax))
 allocate(u(imax+1,jmax+1))
 allocate(v(imax+1,jmax+1))
 allocate(p(imax+1,jmax+1))
 allocate(T(imax+1,jmax+1))
 allocate(rho(imax+1,jmax+1))

 !IBM Variables
 allocate(ibpts(n_ibpts,2))
 allocate(ib_norm(n_ibpts,2))
 allocate(cc(imax+1,jmax+1))

 print*,"allocated"

 call i_normal(domain,i_norm,i_length,ibpts,ib_norm)
 call j_normal(domain,j_norm,j_length)
 call volume_cell(domain,volume)
 call cell_centre(domain,cell_c) 
! call ibm(cell_c,ib_norm,ibpts)
 print*,'Calling forcing'
! call forcing(cell_c,ib_norm,ibpts)
 call initialization(Mach,angle,restart,u,v,p,rho)
 call fluid_properties(cell_c,domain,u,v,p,T,rho,i_length,j_length,i_norm,j_norm,volume,ibpts,ib_norm)



  deallocate(domain)
  deallocate(i_norm)
  deallocate(j_norm)
  deallocate(cell_c,volume)
  deallocate(i_length)
  deallocate(j_length)
  deallocate(u,v,p,rho)
  deallocate(ibpts,ib_norm,cc)
print*,'deallocated'


end program main_program 
