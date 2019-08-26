!------------------------------------------------------
!**************Flux Reconstruction********************
!-----------------------------------------------------

  module flux_calculation 
        use variables
        use geometrical_properties
        use immersed_routine
 contains
     subroutine initialization(Mach,angle,restart,u,v,p,rho)
     implicit none
     
     integer   ::   restart
     real      ::   v_res,dummy1,dummy2,Mach,angle
     real,dimension(:,:)    ::  u,v,p,rho
     v_res=Mach*sqrt(gama*R*T_inf)   !mach number to velocity
     uinf=v_res*cos(angle)
     vinf=v_res*sin(angle)
     rho_inf=pinf/(R*T_inf)

     print*,'Flow Quantities:'
     print*,'Fress stream velocity is :',uinf
     print*,'Free stream Reynolds number is ',(rho_inf*uinf*1.0)/mue
     print*,'Density is :',rho_inf
     print*,'Mach Number is :',Mach

     if(restart.eq.0)then
          do j=2,jmax 
             do i=2,imax
                u(i,j)=uinf
                v(i,j)=vinf
                p(i,j)=pinf
                rho(i,j)=rho_inf
            !    u(i,j)=0.0
            !    v(i,j)=0.0
            !    p(i,j)=pinf
            !    rho(i,j)=rho_inf

             end do
          end do

     !restart the file
     else
          
          120 format('results/ycuttcell',i1,'.dat')
           write(filename_results,120)ipurge
           open(40,file=filename_results)
           print*,'restart' 
           read(40,*)
           read(40,*)
           do j=2,jmax
               do i=2,imax
                  read(40,*)dummy1,dummy2,u(i,j),v(i,j),p(i,j),rho(i,j)
                end do
           end do
           close(40)
      end if
      end subroutine initialization

      subroutine fluid_properties(cell_c,domain,u,v,p,T,rho,i_length,j_length,i_norm,j_norm,volume,ibpts,ib_norm)

        integer::  iteration
        real   ::  a_interface,M_i,M_ip,beta_p,beta_m,alpha_p,alpha_m,cvl_p,cvl_m,D_p,D_m
        real   ::  D_tilda_p,D_tilda_m,cls_p,cls_m,M_j,M_jp,b
        real   ::  u_int1,u_int2,u_int3,u_int4,v_int1,v_int2,v_int3,v_int4
        real   ::  u1,u2,u3,u4,a1,a2,a3,a4,lambda1,lambda2,lambda3,lambda4
        real   ::  scale1,scale2,scale3,scale4,resnorm,resnorm0
        real,dimension(4)   ::   Q_l,Q_r,inormal,jnormal,Q_c,residue
        real,dimension(:,:),allocatable     ::   delta_t,M,ho,a
        real,dimension(:,:,:),allocatable   ::  sum_flux,Fx,Fy
        real,dimension(:,:),allocatable     ::   u_n,v_n,p_n,rho_n,M_n
        real,dimension(:,:,:)    ::   domain
        real,dimension(:,2:,:)   ::   i_norm
        real,dimension(2:,:,:)   ::   j_norm
        real,dimension(:,:,:)    ::   cell_c
        real,dimension(:,2:)     ::   i_length
        real,dimension(2:,:)     ::   j_length
        real,dimension(2:,2:)    ::   volume
        real,dimension(:,:)      ::   u,v,p,rho,T

        !Higher order variables& Viscous flux
        real                       ::   k_tc,Pr=0.721
        real,dimension(5)          ::   p_vi,p_vip,p_vipp,p_vim,p_vl,p_vr
        real,dimension(:,:,:),allocatable      ::   grad_u,grad_ui,grad_uj
        real,dimension(:,:,:),allocatable      ::   grad_v,grad_vi,grad_vj
        real,dimension(:,:,:),allocatable      ::   grad_T,grad_Ti,grad_Tj
        real,dimension(:,:),allocatable        ::   Tau_xx,Tau_yy,Tau_xy
        real,dimension(:,:,:),allocatable      ::   Fv_x,Fv_y
        
        !Immersed Boundary
                                                    
        real,dimension(:,:)       ::  ib_norm,ibpts

        allocate(sum_flux(2:imax,2:jmax,4))
        allocate(delta_t(2:imax,2:jmax))
        allocate(M(imax+1,jmax+1))
        allocate(ho(imax+1,jmax+1))
        allocate(a(imax+1,jmax+1))
        allocate(u_n(imax+1,jmax+1))
        allocate(v_n(imax+1,jmax+1))
        allocate(p_n(imax+1,jmax+1))
        allocate(rho_n(imax+1,jmax+1))
        allocate(M_n(imax+1,jmax+1))
        allocate(Fx(imax,2:jmax,4))
        allocate(Fy(2:imax,jmax,4))

        !Visocus allocation
        allocate(grad_u(imax,jmax,2),grad_ui(imax,jmax,2),grad_uj(imax,jmax,2))
        allocate(grad_v(imax,jmax,2),grad_vi(imax,jmax,2),grad_vj(imax,jmax,2))
        allocate(grad_T(imax,jmax,2),grad_Ti(imax,jmax,2),grad_Tj(imax,jmax,2))
        allocate(Tau_xx(imax,jmax))
        allocate(Tau_yy(imax,jmax))
        allocate(Tau_xy(imax,jmax))
        allocate(Fv_x(imax,2:jmax,4))
        allocate(Fv_y(2:imax,jmax,4))

        no_of_cells=(imax-1)*(jmax-1) 
        residue=(/0.0,0.0,0.0,0.0/) 
 
        if(immersed .eq. 0)then
           do i = 1, imax+1
              do j = 1, jmax+1
                 cc(i,j) = 1

              end do
           end do
        else if(immersed .eq. 1)then
           call ibm(cell_c,ib_norm,ibpts)
        end if


        ! Enforce Boundary conditions on ghost cells
        !    1      =     subsonic inlet
        !    2      =     subsonic outlet
        !    3      =     slip wall
        !    4      =     No slip wall
        !    5      =     Supersonic inlet
        !    6      =     Supersonic outlet
  
        do iteration=1,istop                        
!---------------------------------------------------------------------        
!                 Flow initialization in ghost cells
!--------------------------------------------------------------------
            if(lbc.eq.1) then
                        u(1,2:jmax)=uinf
                        v(1,2:jmax)=vinf
                        p(1,2:jmax)=p(2,2:jmax)
                        rho(1,2:jmax)=rho(2,2:jmax)
                                                                                                                                                                 
             else if (lbc.eq.2) then
                        u(1,2:jmax)=u(2,2:jmax)
                        v(1,2:jmax)=v(2,2:jmax)
                        p(1,2:jmax)=pinf
                        rho(1,2:jmax)=rho_inf
                                                                                                                                                                 
              else if (lbc.eq.3) then
                        u(1,2:jmax)=u(2,2:jmax)-2*(u(2,2:jmax)*i_norm(1,2:jmax,1)&
                                       +v(2,2:jmax)*i_norm(1,2:jmax,2))*i_norm(1,2:jmax,1)
                        v(1,2:jmax)=v(2,2:jmax)-2*(u(2,2:jmax)*i_norm(1,2:jmax,1)&
                                       +v(2,2:jmax)*i_norm(1,2:jmax,2))*i_norm(1,2:jmax,2)
                        p(1,2:jmax)=p(2,2:jmax)
                        rho(1,2:jmax)=rho(2,2:jmax)
               
              else if (lbc.eq.4) then
                        u(1,2:jmax)=-u(2,2:jmax)  
                        v(1,2:jmax)=-v(2,2:jmax)
                        p(1,2:jmax)=p(2,2:jmax)
                        rho(1,2:jmax)=rho(2,2:jmax)
  
              else if(lbc.eq.5) then                                                                                                                                          
                        u(1,2:jmax)=uinf
                        v(1,2:jmax)=vinf
                        p(1,2:jmax)=pinf
                        rho(1,2:jmax)=rho_inf                                                                                                                                           
              else if (lbc.eq.6) then
                        u(1,2:jmax)=u(2,2:jmax)
                        v(1,2:jmax)=v(2,2:jmax)
                        p(1,2:jmax)=p(2,2:jmax)
                        rho(1,2:jmax)=rho(2,2:jmax)
  
              end if
  
              !Right boundary
                                                                                                                                                                 
              if(rbc.eq.1) then
                        u(imax+1,2:jmax)=uinf
                        v(imax+1,2:jmax)=vinf
                        p(imax+1,2:jmax)=p(imax,2:jmax)
                        rho(imax+1,2:jmax)=rho(imax,2:jmax)
                                                                                                
             else if (rbc.eq.2) then
                        u(imax+1,2:jmax)=u(imax,2:jmax)
                        v(imax+1,2:jmax)=v(imax,2:jmax)
                        p(imax+1,2:jmax)=pinf
                        rho(imax+1,2:jmax)=rho_inf
                                                                                                
             else if (rbc.eq.3) then
                        u(imax+1,2:jmax)=u(imax,2:jmax)-2*(u(imax,2:jmax)*i_norm(imax,2:jmax,1)&
                                         +v(imax,2:jmax)*i_norm(imax,2:jmax,2))*i_norm(imax,2:jmax,1)
                        v(imax+1,2:jmax)=v(imax,2:jmax)-2*(u(imax,2:jmax)*i_norm(imax,2:jmax,1)&
                                         +v(imax,2:jmax)*i_norm(imax,2:jmax,2))*i_norm(imax,2:jmax,2)
                        p(imax+1,2:jmax)=p(imax,2:jmax)
                        rho(imax+1,2:jmax)=rho(imax,2:jmax)
                                                                                                
             else if (rbc.eq.4) then
                       u(imax+1,2:jmax)=-u(imax,2:jmax)  
                       v(imax+1,2:jmax)=-v(imax,2:jmax)
                       p(imax+1,2:jmax)=p(imax,2:jmax)
                       rho(imax+1,2:jmax)=rho(imax,2:jmax)
  
             else if(rbc.eq.5) then                                                                                                                                          
                       u(imax+1,2:jmax)=uinf
                       v(imax+1,2:jmax)=vinf
                       p(imax+1,2:jmax)=pinf
                       rho(imax+1,2:jmax)=rho_inf
                                                                                                                                                                   
             else if (rbc.eq.6) then
                       u(imax+1,2:jmax)=u(imax,2:jmax)
                       v(imax+1,2:jmax)=v(imax,2:jmax)
                       p(imax+1,2:jmax)=p(imax,2:jmax)
                       rho(imax+1,2:jmax)=rho(imax,2:jmax)
                       
             end if 
                                                                                                                                                                 
             !top boundary
  
             if(tbc.eq.1) then                                                                                                                               
                       u(2:imax,jmax+1)=uinf
                       v(2:imax,jmax+1)=vinf
                       p(2:imax,jmax+1)=p(2:imax,jmax)
                       rho(2:imax,jmax+1)=rho(2:imax,jmax)
                                                                                                 
              else if (tbc.eq.2) then
                       u(2:imax,jmax+1)=u(2:imax,jmax)
                       v(2:imax,jmax+1)=v(2:imax,jmax)
                       p(2:imax,jmax+1)=pinf
                       rho(2:imax,jmax+1)=rho_inf
  
              else if (tbc.eq.3) then
                                                                                                 
                       u(2:imax,jmax+1)=u(2:imax,jmax)-2*(u(2:imax,jmax)*j_norm(2:imax,jmax,1)&
                                        +v(2:imax,jmax)*j_norm(2:imax,jmax,2))*j_norm(2:imax,jmax,1)
                       v(2:imax,jmax+1)=v(2:imax,jmax)-2*(u(2:imax,jmax)*j_norm(2:imax,jmax,1)&
                                        +v(2:imax,jmax)*j_norm(2:imax,jmax,2))*j_norm(2:imax,jmax,2)
                        p(2:imax,jmax+1)=p(2:imax,jmax)                                                                                                          
                        rho(2:imax,jmax+1)=rho(2:imax,jmax)                                                                                                          
  
              else if (tbc.eq.4) then
                                                                                                 
                       u(2:imax,jmax+1)=-u(2:imax,jmax)  
                       v(2:imax,jmax+1)=-v(2:imax,jmax)
                       p(2:imax,jmax+1)=p(2:imax,jmax)
                       rho(2:imax,jmax+1)=rho(2:imax,jmax)
  
  
              else if(tbc.eq.5) then                                                                                                                                          
                       u(2:imax,jmax+1)=uinf
                       v(2:imax,jmax+1)=vinf
                       p(2:imax,jmax+1)=pinf
                       rho(2:imax,jmax+1)=rho_inf
                                                                                                                                                                    
              else if (tbc.eq.6) then
                       u(2:imax,jmax+1)=u(2:imax,jmax)
                       v(2:imax,jmax+1)=v(2:imax,jmax)
                       p(2:imax,jmax+1)=p(2:imax,jmax)
                       rho(2:imax,jmax+1)=rho(2:imax,jmax)
  
              end if 
              
             !Bottom Boundary
                                                                                                                                                                 
             if(bbc.eq.1) then                                                                                                                               
                       u(2:imax,1)=uinf
                       v(2:imax,1)=vinf
                       p(2:imax,1)=p(2:imax,2)
                       rho(2:imax,1)=rho(2:imax,2)
                                                                                                      
             else if (bbc.eq.2) then
                       u(2:imax,1)=u(2:imax,2)
                       v(2:imax,1)=v(2:imax,2)
                       p(2:imax,1)=pinf
                       rho(2:imax,1)=rho_inf
                                                                                                      
             else if (bbc.eq.3) then
                       u(2:imax,1)=u(2:imax,2)-2*(u(2:imax,2)*j_norm(2:imax,1,1)&
                                     +v(2:imax,2)*j_norm(2:imax,1,2))*j_norm(2:imax,1,1)
                       v(2:imax,1)=v(2:imax,2)-2*(u(2:imax,2)*j_norm(2:imax,1,1)&
                                    +v(2:imax,2)*j_norm(2:imax,1,2))*j_norm(2:imax,1,2)
                       p(2:imax,1)=p(2:imax,2)                                                                                                          
                       rho(2:imax,1)=rho(2:imax,2)                                                                                                          
                                                                                                      
             else if (bbc.eq.4) then
                       u(2:imax,1)=-u(2:imax,2)  
                       v(2:imax,1)=-v(2:imax,2)
                       p(2:imax,1)=p(2:imax,2)
                       rho(2:imax,1)=rho(2:imax,2)
  
             else if(bbc.eq.5) then                                                                                                                                          
                       u(2:imax,1)=uinf
                       v(2:imax,1)=vinf
                       p(2:imax,1)=pinf
                       rho(2:imax,1)=rho_inf
                                                                                                                                                                   
             else if (bbc.eq.6) then
                       u(2:imax,1)=u(2:imax,2)
                       v(2:imax,1)=v(2:imax,2)
                       p(2:imax,1)=p(2:imax,2)
                       rho(2:imax,1)=rho(2:imax,2)
             end if           
            
  
             T(:,:)=p(:,:)/(rho(:,:)*R)
             ho(:,:)=gama*R*T(:,:)/(gama-1)+0.5*(u(:,:)**2+v(:,:)**2)
             a(:,:)=sqrt(gama*R*T(:,:))
  
             b=(3.0-kappa)/(1.0-kappa)

 !---------------------------------------------- ------------------------
 !                         Flux in x direction
 !------------------------------------------------------------------------
             do j=2,jmax
                do i=1,imax
            
                 if(cc(i,j) .eq. 1 .or. cc(i+1,j) .eq. 1)then
                    p_vi(:)=(/rho(i,j),u(i,j),v(i,j),p(i,j),ho(i,j)/)
                    p_vip(:)=(/rho(i+1,j),u(i+1,j),v(i+1,j),p(i+1,j),ho(i+1,j)/)
   
          !order 1- First order
          !       0- Second order

                  if(order .eq. 1)then
                     p_vl(:)=p_vi(:)
                     p_vr(:)=p_vip(:)
  
                  elseif(order .eq. 0)then
                     if(i .eq. 1)then
                        p_vim(:)=p_vi(:)
                     elseif(i .eq. imax)then
                        p_vipp(:)=p_vip(:)
                     else
                        p_vim(:)=(/rho(i-1,j),u(i-1,j),v(i-1,j),p(i-1,j),ho(i-1,j)/)
                        p_vipp(:)=(/rho(i+2,j),u(i+2,j),v(i+2,j),p(i+2,j),ho(i+2,j)/)
                      end if
  
                      p_vl(:)=p_vi(:)+(1-order)*(0.25*((1-kappa)*TVD(p_vi(:)-p_vim(:),&
                              b*(p_vip(:)-p_vi(:)))+(1+kappa)*&
                              TVD(b*(p_vi(:)-p_vim(:)),p_vip(:)-p_vi(:))))
                                                                                                                                                                                                                                        
                      p_vr(:)=p_vip(:)-(1-order)*(0.25*((1-kappa)*TVD(b*(p_vip(:)-p_vi(:)),&
                              p_vipp(:)-p_vip(:))+(1+kappa)*&
                              TVD(b*(p_vi(:)-p_vim(:)),p_vip(:)-p_vi(:))))
                  end if
  
  
                  Q_l(:)=(/1.0,p_vl(2),p_vl(3),p_vl(5)/) 
                  Q_r(:)=(/1.0,p_vr(2),p_vr(3),p_vr(5)/) 
            
                  inormal(:)=(/0.0,i_norm(i,j,1),i_norm(i,j,2),0.0/)
                  a_interface=(a(i,j)+a(i+1,j))*0.5
  
                  M_i=(p_vl(2)*i_norm(i,j,1)+p_vl(3)*i_norm(i,j,2))/a_interface
                  M_ip=(p_vr(2)*i_norm(i,j,1)+p_vr(3)*i_norm(i,j,2))/a_interface
   
  
                  beta_p=-max(0.0,1.0-int(abs(M_i))) 
                  beta_m=-max(0.0,1.0-int(abs(M_ip)))
                  alpha_p=0.5*(1.0+sign(1.0,M_i))
                  alpha_m=0.5*(1.0-sign(1.0,M_ip))
  
                  cvl_p=alpha_p*(1.0+beta_p)*M_i-(0.25*beta_p*(1.0+M_i)**2) 
                  cvl_m=alpha_m*(1.0+beta_m)*M_ip+(0.25*beta_m*(1.0-M_ip)**2) 
                  D_p=0.25*(2.0-M_i)*(1.0+M_i)**2
                  D_m=0.25*(2.0+M_ip)*(1.0-M_ip)**2
                  
                  D_tilda_p=alpha_p*(1.0+beta_p)-(beta_p*D_p)
                  D_tilda_m=alpha_m*(1.0+beta_m)-(beta_m*D_m)
                  
                  cls_p=max(0.0,cvl_p+cvl_m)
                  cls_m=min(0.0,cvl_p+cvl_m)
        
  
                  Fx(i,j,:)=(p_vl(1)*a_interface*cls_p*Q_l(:)+p_vr(1)*a_interface*cls_m*Q_r(:) &
                            +(D_tilda_p*p_vl(4)+D_tilda_m*p_vr(4))*inormal(:))*i_length(i,j)

                  
                  
                  !sanity check for NAN
                  if(isnan(Fx(i,j,4))) then
                    print*,Fx(i,j,:)
                    print*,'NAN is encountered'
                    stop
                 end if
               end if 
               end do
            end do           
!--------------------------------------------------------------------------------         
!                        flux in y direction
!----------------------------------------------------------------------------
        
            do j=1,jmax
               do i=2,imax
  
                   if(cc(i,j) .eq. 1 .or. cc(i,j+1) .eq. 1)then
                    p_vi(:)=(/rho(i,j),u(i,j),v(i,j),p(i,j),ho(i,j)/)                           
                    p_vip(:)=(/rho(i,j+1),u(i,j+1),v(i,j+1),p(i,j+1),ho(i,j+1)/)
  
  
                    if(order .eq. 1)then
                       p_vl(:)=p_vi(:)
                       p_vr(:)=p_vip(:)
  
                    elseif(order .eq. 0)then
                       if(j .eq. 1)then
                          p_vim(:)=p_vi(:)
                       elseif(j .eq. jmax)then
                          p_vipp(:)=p_vip(:)
                       else
                          p_vim(:)=(/rho(i,j-1),u(i,j-1),v(i,j-1),p(i,j-1),ho(i,j-1)/)
                          p_vipp(:)=(/rho(i,j+2),u(i,j+2),v(i,j+2),p(i,j+2),ho(i,j+2)/)
                       end if                                                                                            
                                                                                                                          
                       p_vl(:)=p_vi(:)+(1-order)*(0.25*((1-kappa)*TVD(p_vi(:)-p_vim(:),&
                               b*(p_vip(:)-p_vi(:)))+(1+kappa)*&
                               TVD(b*(p_vi(:)-p_vim(:)),p_vip(:)-p_vi(:))))
                                                                                                              
                       p_vr(:)=p_vip(:)-(1-order)*(0.25*((1-kappa)*TVD(b*(p_vip(:)-p_vi(:)),&
                               p_vipp(:)-p_vip(:))+(1+kappa)*&
                               TVD(b*(p_vi(:)-p_vim(:)),p_vip(:)-p_vi(:))))
                                                                
                     end if
                   
                   
                 Q_l(:)=(/1.0,p_vl(2),p_vl(3),p_vl(5)/) 
                 Q_r(:)=(/1.0,p_vr(2),p_vr(3),p_vr(5)/) 
                 jnormal(:)=(/0.0,j_norm(i,j,1),j_norm(i,j,2),0.0/)
                 
                 a_interface=(a(i,j)+a(i,j+1))*0.5                                           
                 M_j=(p_vl(2)*j_norm(i,j,1)+p_vl(3)*j_norm(i,j,2))/a_interface
                 M_jp=(p_vr(2)*j_norm(i,j,1)+p_vr(3)*j_norm(i,j,2))/a_interface
  
                 beta_p=-max(0.0,1.0-int(abs(M_j))) 
                 beta_m=-max(0.0,1.0-int(abs(M_jp)))
                 alpha_p=0.5*(1.0+sign(1.0,M_j))                                                                                                                                                                                         
                 alpha_m=0.5*(1.0-sign(1.0,M_jp))
  
                 cvl_p=alpha_p*(1.0+beta_p)*M_j-(0.25*beta_p*(1.0+M_j)**2) 
                 cvl_m=alpha_m*(1.0+beta_m)*M_jp+(0.25*beta_m*(1.0-M_jp)**2) 
                 D_p=0.25*(2.0-M_j)*(1.0+M_j)**2
                 D_m=0.25*(2.0+M_jp)*(1.0-M_jp)**2
  
                 D_tilda_p=alpha_p*(1.0+beta_p)-(beta_p*D_p)
                 D_tilda_m=alpha_m*(1.0+beta_m)-(beta_m*D_m)
                 cls_p=max(0.,cvl_p+cvl_m)
                 cls_m=min(0.,cvl_p+cvl_m)
   
  
                 Fy(i,j,:)=(p_vl(1)*a_interface*cls_p*Q_l(:)+p_vr(1)*a_interface*cls_m*Q_r(:)+& 
                           (D_tilda_p*p_vl(4)+D_tilda_m*p_vr(4))*jnormal(:))*j_length(i,j)
 


               end if 
               end do
            end do 
!----------------------------------------------------------------------------
!                      Viscous-Flux
!---------------------------------------------------------------------------
            if(viscous .eq. 1)then
              
               k_tc=(mue*1005.4)/Pr
               do j=2,jmax
                  do i=2,imax
                     u_int1=(u(i,j)+u(i+1,j))*0.5
                     u_int2=(u(i,j)+u(i-1,j))*0.5
                     u_int3=(u(i,j)+u(i,j+1))*0.5
                     u_int4=(u(i,j)+u(i,j-1))*0.5

                     v_int1=(v(i,j)+v(i+1,j))*0.5
                     v_int2=(v(i,j)+v(i-1,j))*0.5            ! Interface Velocity
                     v_int3=(v(i,j)+v(i,j+1))*0.5
                     v_int4=(v(i,j)+v(i,j-1))*0.5

                     T_int1=(T(i,j)+T(i+1,j))*0.5
                     T_int2=(T(i,j)+T(i-1,j))*0.5
                     T_int3=(T(i,j)+T(i,j+1))*0.5
                     T_int4=(T(i,j)+T(i,j-1))*0.5

               !Gradients at cell centres
                     grad_u(i,j,:)=(u_int1*i_length(i,j)*i_norm(i,j,:)-u_int2*i_length(i-1,j)*i_norm(i-1,j,:)+&
                                    u_int3*j_length(i,j)*j_norm(i,j,:)-u_int4*j_length(i,j-1)*j_norm(i,j-1,:))/volume(i,j)

                     grad_v(i,j,:)=(v_int1*i_length(i,j)*i_norm(i,j,:)-v_int2*i_length(i-1,j)*i_norm(i-1,j,:)+&
                                    v_int3*j_length(i,j)*j_norm(i,j,:)-v_int4*j_length(i,j-1)*j_norm(i,j-1,:))/volume(i,j)     

                     grad_T(i,j,:)=(T_int1*i_length(i,j)*i_norm(i,j,:)-T_int2*i_length(i-1,j)*i_norm(i-1,j,:)+&
                                    T_int3*j_length(i,j)*j_norm(i,j,:)-T_int4*j_length(i,j-1)*j_norm(i,j-1,:))/volume(i,j)

                  end do
               end do

               ! Gradients at the interface
               grad_ui(1,:,:)=grad_u(2,:,:)
               grad_vi(1,:,:)=grad_v(2,:,:)
               grad_Ti(1,:,:)=grad_T(2,:,:)

               grad_ui(imax,:,:)=grad_u(imax,:,:)
               grad_vi(imax,:,:)=grad_v(imax,:,:)
               grad_Ti(imax,:,:)=grad_T(imax,:,:)     !Left side face of Right ghost cell

               open(55,file='results/x-gradients.dat')
               write(55,*)'Variables=  x,   y, du_dx,   du_dy,   dv_dx,   dv_dy'
               write(55,*)'zone i=', imax-2,'    j='   ,  jmax-1

               do j=2,jmax
                  do i=2,imax-1
                     grad_ui(i,j,:)=(grad_u(i,j,:)+grad_u(i+1,j,:))*0.5
                     grad_vi(i,j,:)=(grad_v(i,j,:)+grad_v(i+1,j,:))*0.5
                     grad_Ti(i,j,:)=(grad_T(i,j,:)+grad_T(i+1,j,:))*0.5     ! x-Interface Gradients
                     write(55,*)domain(i,j,:),grad_ui(i,j,:),grad_vi(i,j,:)
                  end do
               end do
               close(55)
               
               do j=2,jmax
                  do i=1,imax
                     Tau_xx(i,j)=2*mue*(grad_ui(i,j,1)-(grad_ui(i,j,1)+grad_vi(i,j,2))/3.0)
                     Tau_xy(i,j)=mue*(grad_ui(i,j,2)+grad_vi(i,j,1))
                     Tau_yy(i,j)=2*mue*(grad_vi(i,j,2)-(grad_ui(i,j,1)+grad_vi(i,j,2))/3.0)
                                                                                               
                     Fv_x(i,j,:)=((/0.0,Tau_xx(i,j),Tau_xy(i,j),0.5*Tau_xx(i,j)*(u(i,j)+u(i+1,j))+&
                                  0.5*Tau_xy(i,j)*(v(i,j)+v(i+1,j))+k_tc*grad_Ti(i,j,1)/)*i_norm(i,j,1)+ &
                                  (/0.0,Tau_xy(i,j),Tau_yy(i,j),0.5*Tau_xy(i,j)*(u(i,j)+u(i+1,j))+&
                                  0.5*Tau_yy(i,j)*(v(i,j)+v(i+1,j))+k_tc*grad_Ti(i,j,2)/)*i_norm(i,j,2))*i_length(i,j)
                                                                                               
                  end do
               end do

               grad_uj(:,1,:)=grad_u(:,2,:)
               grad_vj(:,1,:)=grad_v(:,2,:)
               grad_Tj(:,1,:)=grad_T(:,2,:)
                                                                             
               grad_uj(:,jmax,:)=grad_u(:,jmax,:)
               grad_vj(:,jmax,:)=grad_v(:,jmax,:)   ! bottom face of top ghost cell
               grad_Tj(:,jmax,:)=grad_T(:,jmax,:)

               open(56,file='results/y-gradients.dat')
               write(56,*)'Variables=  x,   y, du_dx,   du_dy,   dv_dx,   dv_dy'
               write(56,*)'zone i=', imax-1,'    j='   ,  jmax-2

               do i=2,imax
                  do j=2,jmax-1
                     grad_uj(i,j,:)=(grad_u(i,j,:)+grad_u(i,j+1,:))*0.5
                     grad_vj(i,j,:)=(grad_v(i,j,:)+grad_v(i,j+1,:))*0.5     ! y-interface Gradients
                     grad_Tj(i,j,:)=(grad_T(i,j,:)+grad_T(i,j+1,:))*0.5 
                     write(56,*)domain(i,j,:),grad_uj(i,j,:),grad_vj(i,j,:)
                  end do
               end do
               close(56)
               
               do i=2,imax
                  do j=1,jmax
                     Tau_xy(i,j)=mue*(grad_uj(i,j,2)+grad_vj(i,j,1))
                     Tau_yy(i,j)=2*mue*(grad_vj(i,j,2)-(grad_uj(i,j,1)+grad_vj(i,j,2))/3.0)
                     Tau_xx(i,j)=2*mue*(grad_uj(i,j,1)-(grad_uj(i,j,1)+grad_vj(i,j,2))/3.0)

                     Fv_y(i,j,:)=((/0.0,Tau_xy(i,j),Tau_yy(i,j),0.5*Tau_xy(i,j)*(u(i,j)+u(i,j+1))+&
                                  0.5*Tau_yy(i,j)*(v(i,j)+v(i,j+1))+k_tc*grad_Tj(i,j,2)/)*j_norm(i,j,2)+ &
                                (/0.0,Tau_xx(i,j),Tau_xy(i,j),0.5*Tau_xx(i,j)*(u(i,j)+u(i,j+1))+&
                                  0.5*Tau_xy(i,j)*(v(i,j)+v(i,j+1))+k_tc*grad_Tj(i,j,1)/)*j_norm(i,j,1))*j_length(i,j)

                  end do
               end do

               Fx(:,2:,:)=Fx(:,2:,:)-Fv_x(:,2:,:)
               Fy(2:,:,:)=Fy(2:,:,:)-Fv_y(2:,:,:)
          end if

!----------------------------------------------------------------------------
  !                    Local time stepping
!----------------------------------------------------------------------------
   
           do j=2,jmax
              do i=2,imax
  
                 u_int1=(u(i,j)+u(i+1,j))*0.5                     
                 u_int2=(u(i,j)+u(i-1,j))*0.5
                 u_int3=(u(i,j)+u(i,j+1))*0.5
                 u_int4=(u(i,j)+u(i,j-1))*0.5
                                                                  
                 v_int1=(v(i,j)+v(i+1,j))*0.5    
                 v_int2=(v(i,j)+v(i-1,j))*0.5
                 v_int3=(v(i,j)+v(i,j+1))*0.5
                 v_int4=(v(i,j)+v(i,j-1))*0.5
                                                                  
                 u1=u_int1*i_norm(i,j,1)+v_int1*i_norm(i,j,2)
                 u2=u_int2*i_norm(i-1,j,1)+v_int2*i_norm(i-1,j,2)
                 u3=u_int3*j_norm(i,j,1)+v_int3*j_norm(i,j,2)
                 u4=u_int4*j_norm(i,j-1,1)+v_int4*j_norm(i,j-1,2)
  
                 a1=(a(i,j)+a(i+1,j))*0.5
                 a2=(a(i,j)+a(i-1,j))*0.5
                 a3=(a(i,j)+a(i,j+1))*0.5
                 a4=(a(i,j)+a(i,j-1))*0.5
  
                 lambda1=u1+a1
                 lambda2=u2+a2
                 lambda3=u3+a3
                 lambda4=u4+a4
  
                 
                 delta_t(i,j)=(0.5*(i_length(i,j)*lambda1+i_length(i-1,j)*lambda2+&
                                   j_length(i,j)*lambda3+j_length(i,j-1)*lambda4))/CFL
               end do
             end do
  !---------------------------------------------------------------------
  !                   Solution Update
  !----------------------------------------------------------------------
            do i=2,imax
               do j=2,jmax
  
                  ! primitive to conservative 
   
                  if(cc(i,j) .eq. 1)then
                 Q_c(:)=(/rho(i,j),rho(i,j)*u(i,j),rho(i,j)*v(i,j),p(i,j)/(gama-1)+0.5*rho(i,j)*(u(i,j)**2+v(i,j)**2)/)    
                 sum_flux(i,j,:)=Fx(i,j,:)-Fx(i-1,j,:)+Fy(i,j,:)-Fy(i,j-1,:)
                 Q_c(:)=Q_c(:)-(sum_flux(i,j,:))/delta_t(i,j)
                 residue(:)=residue(:)+abs(sum_flux(i,j,:))
 
                 !conservative to primitive
  
                 rho(i,j)=Q_c(1)
                 u(i,j)=Q_c(2)/Q_c(1)
                 v(i,j)=Q_c(3)/Q_c(1)
                 p(i,j)=(gama-1)*(Q_c(4)-0.5*(Q_c(2)**2+Q_c(3)**2)/Q_c(1))
                 T(i,j)=p(i,j)/(rho(i,j)*R)
                 M(i,j)=sqrt(u(i,j)**2+v(i,j)**2)/sqrt(gama*R*T(i,j))


                   else if(cc(i,j).eq. 0)then 
                        call forcing(i,j,cell_c,ib_norm,ibpts,u,v,p,T,rho) 
                  end if



               end do
            end do

!-----------------------------------------------------------------------
!                        Resnorm  Calculation
!----------------------------------------------------------------------
           open(47,file='results/Resnorm')
           open(51,file='results/resnorm.dat')
           resnorm=0.0
           scale1=rho_inf
           scale2=rho_inf*sqrt(uinf**2+vinf**2)
           scale3=scale2
           scale4=pinf/(gama-1)+0.5*rho_inf*sqrt(uinf**2+vinf**2)
           do j=2,jmax
              do i=2,imax
                 resnorm=resnorm  +  (sum_flux(i,j,1)/scale1)**2 +&
                                     (sum_flux(i,j,2)/scale2)**2 +&
                                     (sum_flux(i,j,3)/scale3)**2 +&
                                     (sum_flux(i,j,4)/scale4)**2 
              end do
           end do
           resnorm=sqrt(resnorm)
           if(iteration .eq. 1 .and. restart .eq.0)then
              resnorm0=resnorm
              write(47,209)resnorm0
           else if(iteration .eq. 1)then
              read(47,*)resnorm0
          end if
          209 format(f15.7)
          print*,iteration,resnorm/resnorm0
          write(51,215)iteration,resnorm/resnorm0
          215 format(i15,f15.7)
          close(47)
  !-------------------------------------------------------------------
  !                      writing data
  !-------------------------------------------------------------------
  
            residue(:)= residue(:)/no_of_cells
           
            open(2,file='results/residual.dat')
            100 format(1i10,4e30.15)
            write(2,100)iteration,residue
         !   print*,residue(:),iteration               
            120 format('results/ycuttcell',i1,'.dat')
            130 format('results/ycuttnode',i1,'.dat')
            if(mod(iteration,irest) .eq. 0.0)then
                  write(filename_results,120)h
                  open(16,file=filename_results)
                  write(16,*)'Variables=   x, y, u , v, p,  rho,   M'
                  write(16,*)'zone i=', imax-1,'    j='   ,  jmax-1
                  do j=2,jmax                                                                 
                                                                          
                     do i=2,imax
                                                                               
                        write(16,*)cell_c(i,j,1:2),u(i,j),v(i,j),p(i,j),rho(i,j),M(i,j)
                                                                               
                      end do
                   end do
                 
  !----------------------------------------------------
  !                       Node data
  !---------------------------------------------------
                  write(filename_results,130)h
                  open(17,file=filename_results)
                  write(17,*)'Variables=   x, y, u , v, p,  rho,   M'
                  write(17,*)'zone i=', imax,'    j='   ,  jmax
  
                  M(1,:)=sqrt(u(1,:)**2+v(1,:)**2)/sqrt(gama*R*T(1,:))
                  M(imax+1,:)=sqrt(u(imax+1,:)**2+v(imax+1,:)**2)/sqrt(gama*R*T(imax+1,:))
                  M(:,1)=sqrt(u(:,1)**2+v(:,1)**2)/sqrt(gama*R*T(:,1))
                  M(:,jmax+1)=sqrt(u(:,jmax+1)**2+v(:,jmax+1)**2)/sqrt(gama*R*T(:,jmax+1))
  
                  do j=1,jmax
                     do i=1,imax
                        if(i .eq. 1 .and. j .eq. 1)then
                           u_n(i,j)=(u(i+1,j+1)+u(i+1,j)+u(i,j+1))/3.0
                           v_n(i,j)=(v(i+1,j+1)+v(i+1,j)+v(i,j+1))/3.0
                           p_n(i,j)=(p(i+1,j+1)+p(i+1,j)+p(i,j+1))/3.0
                           rho_n(i,j)=(rho(i+1,j+1)+rho(i+1,j)+rho(i,j+1))/3.0
                           M_n(i,j)=(M(i+1,j+1)+M(i+1,j)+M(i,j+1))/3.0
                        elseif(i .eq. 1 .and. j .eq. jmax)then
                          u_n(i,j)=(u(i,j)+u(i+1,j)+u(i+1,j+1))/3.0
                          v_n(i,j)=(v(i,j)+v(i+1,j)+v(i+1,j+1))/3.0
                          p_n(i,j)=(p(i,j)+p(i+1,j)+p(i+1,j+1))/3.0
                          rho_n(i,j)=(rho(i,j)+rho(i+1,j)+rho(i+1,j+1))/3.0
                          M_n(i,j)=(M(i,j)+M(i+1,j)+M(i+1,j+1))/3.0
                        elseif(i .eq. imax .and. j .eq. 1)then
                          u_n(i,j)=(u(i,j)+u(i+1,j+1)+u(i,j+1))/3.0
                          v_n(i,j)=(v(i,j)+v(i+1,j+1)+v(i,j+1))/3.0
                          p_n(i,j)=(p(i,j)+p(i+1,j+1)+p(i,j+1))/3.0
                          rho_n(i,j)=(rho(i,j)+rho(i+1,j+1)+rho(i,j+1))/3.0
                          M_n(i,j)=(M(i,j)+M(i+1,j+1)+M(i,j+1))/3.0
                       elseif(i .eq. imax .and. j .eq. jmax)then
                         u_n(i,j)=(u(i,j)+u(i+1,j)+u(i,j+1))/3.0
                         v_n(i,j)=(v(i,j)+v(i+1,j)+v(i,j+1))/3.0
                         p_n(i,j)=(p(i,j)+p(i+1,j)+p(i,j+1))/3.0
                         rho_n(i,j)=(rho(i,j)+rho(i+1,j)+rho(i,j+1))/3.0
                         M_n(i,j)=(M(i,j)+M(i+1,j)+M(i,j+1))/3.0
                      else
                        u_n(i,j)=0.25*(u(i,j)+u(i+1,j)+u(i,j+1)+u(i+1,j+1))
                        v_n(i,j)=0.25*(v(i,j)+v(i+1,j)+v(i,j+1)+v(i+1,j+1))
                        p_n(i,j)=0.25*(p(i,j)+p(i+1,j)+p(i,j+1)+p(i+1,j+1))
                        rho_n(i,j)=0.25*(rho(i,j)+rho(i+1,j)+rho(i,j+1)+rho(i+1,j+1))
                        M_n(i,j)=0.25*(M(i,j)+M(i+1,j)+M(i,j+1)+M(i+1,j+1))
                    end if
                       
                    write(17,*)domain(i,j,1:2),u_n(i,j),v_n(i,j),p_n(i,j),rho_n(i,j),M_n(i,j)
  
                   end do
                end do
            
                 print*,   'continuity      x-momentum       y-momentum         Energy         Iteration'
                 close(16)
                 close(17)
                 h=h+1
                 if(h .gt. iprint)then
                    h=1
                 end if
             else 
                   continue
             end if 
             if (maxval(residue).lt.10e-24)then
                   print*, "solution is converged"
                   goto 110
             else 
                   continue                                        
             end if
                                      
        end do
        close(2)         
        close(51)
               
        110  write(filename_results,120)h
             open(16,file=filename_results)
             
             write(16,*)'Variables=   x, y, u , v, p,  rho,    M'
             write(16,*)'zone i=', imax-1,'    j='   ,  jmax-1
                                                               
         do j=2,jmax
                                                               
            do i=2,imax
                                                               
               write(16,*)cell_c(i,j,1:2),u(i,j),v(i,j),p(i,j),rho(i,j),M(i,j)
                                                               
            end do
         end do
         close(16)
  
         deallocate(sum_flux,Fx,Fy,delta_t,ho,a,M,u_n,v_n,p_n,rho_n,M_n)
         deallocate(grad_u,grad_ui,grad_uj,Fv_x,Fv_y,Tau_xx,Tau_xy,Tau_yy)
         deallocate(grad_v,grad_vi,grad_vj,grad_T,grad_Ti,grad_Tj)
    end subroutine fluid_properties

!--------------------------------------------------------------
!                MIN-MOD Limiter
!--------------------------------------------------------------

    function TVD(x,y)
    !    !  Takes (i plus-i),(i-i minus) as input
    !    !  And gives one output  ("MIN_MOD")
        implicit none
        real,dimension(5),intent(in) ::  x,y
        real, dimension(5)           ::  TVD
        integer                      ::  g
                                                    
        do g=1,5
          if(x(g)*y(g) .gt. 0.0)then
              if(abs(x(g)) .le. abs(y(g)))then
                 TVD(g)=x(g)
              else if(abs(x(g)) .gt. abs(y(g)))then
                 TVD(g)=y(g)
              end if
           else
             TVD(g)=0.0
           end if
       end do
     end function TVD
  
  
end module flux_calculation
