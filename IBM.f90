module immersed_routine
  use variables
  use geometrical_properties
  implicit none
  integer    ::   n_bc,n_fc
  integer    ::   ii,jj,bn = 1,fn = 1
contains 
  subroutine ibm(cell_c,ib_norm,ibpts)
     
     real    ::   xmax,xmin,ymax,ymin,phi
     real,dimension(:,:,:)     ::  cell_c
     real,dimension(:,:)       ::  ib_norm,ibpts
     real,dimension(:), allocatable  ::  dist

     allocate(dist(n_ibpts))
     xmax = maxval(ibpts(:,1))
     xmin = minval(ibpts(:,1))
     ymax = maxval(ibpts(:,2))
     ymin = minval(ibpts(:,2))

     xmax = 2*xmax
     xmin = 2*xmin
     ymax = 2*ymax
     ymin = 2*ymin
    
     !Algorithm to classify field cells

     open(21,file = 'results/Field_cells')    
     open(22,file = 'results/Interior_cells')    
     open(23,file = 'results/Band_cells')   

     if(closed_profile .eq. 1)then
          do i = imax,2,-1
             if(cell_c(i,2,1) .lt.  xmax)then
                imax_fc = i
                exit 
             end if
          end do

          do i = 2,imax
             if(cell_c(i,2,1) .gt. xmin)then
                imin_fc = i
                exit
             end if
          end do

          do j = jmax,2,-1
             if(cell_c(2,j,2) .lt.  ymax)then
                jmax_fc = j
                exit
             end if
         end do

         do j = 2,jmax
            if(cell_c(2,j,2) .gt. ymin)then
              jmin_fc = j
              exit
            end if
         end do

         !writing Field cells data
         do i = 2,imin_fc
            do j = 2,jmax
               cc(i,j) = 1
               n_fc = n_fc+1
               write(21,*)cell_c(i,j,:)
            end do
         end do

         do i = imax_fc+1,imax
            do j = 2,jmax
               cc(i,j) = 1
               n_fc = n_fc+1
               write(21,*)cell_c(i,j,:)
            end do
         end do

         do j = 2,jmin_fc
            do i = imin_fc+1,imax_fc-1
              cc(i,j) = 1
              n_fc = n_fc+1
              write(21,*)cell_c(i,j,:)
           end do
        end do

        do j = jmax_fc+1,jmax
           do i = imin_fc+1,imax_fc-1
              cc(i,j) = 1
              n_fc = n_fc+1
              write(21,*)cell_c(i,j,:)
           end do
        end do

        !Ghost cells classification
        cc(1,2:jmax) = 1
        cc(2:imax,1) = 1
        cc(imax,2:jmax) = 1
        cc(2:imax,jmax) = 1

   
        print*,'Classifying Field and Interior cells'
        allocate(c_ib(imin_fc:imax_fc,jmin_fc:jmax_fc))
            do i = imin_fc+1,imax_fc-1
               do j = jmin_fc+1,jmax_fc-1
                  do k = 1,n_ibpts
                     dist(k) = sqrt((ibpts(k,1)-cell_c(i,j,1))**2+(ibpts(k,2)-cell_c(i,j,2))**2)
                  end do
                  c_ib(i,j) = minloc(dist,dim = 1)
                  phi = (cell_c(i,j,1)-ibpts(c_ib(i,j),1))*ib_norm(c_ib(i,j),1)+&
                            (cell_c(i,j,2)-ibpts(c_ib(i,j),2))*ib_norm(c_ib(i,j),2)
                  if(phi .gt. 0.0)then
                       cc(i,j) = 1
                       n_fc = n_fc+1
                       write(21,*)cell_c(i,j,:)
                  elseif(phi.lt. 0.0)then
                       cc(i,j) = -1
                       write(22,*)cell_c(i,j,:)
                  end if
               end do
             end do
             print*,c_ib(140,34) 
             print*,'Band Cell Classification'

             ! Band cell Classification

             do i = imin_fc+1,imax_fc-1
                do j = jmin_fc+1,jmax_fc-1                
                   if(cc(i,j) .eq. 1)then
                       do ii = i-1,i+1
                          do jj = j-1,j+1
                             if(cc(ii,jj) .eq. -1)then
                                  cc(i,j) = 0
                                  n_bc = n_bc+1
                                  write(23,*)cell_c(i,j,:)
                                  goto 110
                              end if
                           end do
                        end do
             110     end if   
                end do
             end do

     else

       ! Field cells and interior cell classification for open Profile Geometry 

        allocate(c_ib(imax,jmax))

        print*,'Classifying Field and Interior cells for open geometry'

        do i = 1,imax
           do j = 1,jmax
             do k = 1,n_ibpts
                dist(k) = sqrt((ibpts(k,1)-cell_c(i,j,1))**2+(ibpts(k,2)-cell_c(i,j,2))**2)
             end do
             c_ib(i,j) = minloc(dist,dim = 1)
             phi = (cell_c(i,j,1)-ibpts(c_ib(i,j),1))*ib_norm(c_ib(i,j),1)+&              
                       (cell_c(i,j,2)-ibpts(c_ib(i,j),2))*ib_norm(c_ib(i,j),2)

             if(phi.gt. 0.0)then
                  cc(i,j) = 1
                  n_fc = n_fc+1
                  write(21,*)cell_c(i,j,:)
             elseif(phi .lt. 0.0)then
                   cc(i,j) = -1
                   write(22,*)cell_c(i,j,:)
             end if
           end do
        end do
       
        ! Band cells classification
        print*,'Band Cell Classification'
        do i = 1,imax
           do j = 1,jmax
              if(cc(i,j) .eq. 1)then
                 do ii = i-1,i+1
                    do jj = j-1,j+1
                       if(cc(ii,jj) .eq. -1)then
                          cc(i,j) = 0
                          n_bc = n_bc+1
                          write(23,*)cell_c(i,j,:)
                          goto 120
                       end if
                    end do
                 end do
          120  end if
           end do
        end do
      
      end if


      !Field cells index           
      allocate(fc(n_fc,2))
      
      do i = 2,imax
         do j = 2,jmax
            if(cc(i,j) .eq. 1)then
               fc(fn,1) = i
               fc(fn,2) = j
               fn = fn+1
            end if
         end do
      end do

      !  allocate band cells and assign their index of them
      allocate(bc(n_bc,2))
      do i = imin_fc+1,imax_fc-1
         do j = jmin_fc+1,jmax_fc-1
            if(cc(i,j) .eq. 0)then
                bc(bn,1) = i
                bc(bn,2) = j
                bn = bn+1
            end if
         end do
      end do
     !  Writing Cell Classification

      open(24,file = 'results/cells.dat')
      write(24,*)'Variables =   x,   y,   cell'
      write(24,*)'zone i = ', imax-1,'    j = '   ,  jmax-1
      do j = 2,jmax
         do i = 2,imax
           write(24,*)cell_c(i,j,:),cc(i,j)
         end do
     end do

     close(21)
     close(22)
     close(23)
     close(24)

 !    deallocate(fc,dist)         

  end subroutine ibm

  subroutine forcing(a,b,cell_c,ib_norm,ibpts,u,v,p,T,rho)
     
     integer                   ::  a,b
     real                      ::  sum2,u_ib,v_ib,g,d_hyp,d_par,d_per       
     real,dimension(5)         ::  sum1
     real,dimension(:,:)       ::  ibpts,ib_norm
     real,dimension(:,:,:)     ::  cell_c
     real,dimension(:,:)       ::  u,v,p,T,rho
     real,dimension(:,:,:),allocatable   ::  d_ip

     allocate(d_ip(500,500,5))
     ! Finding  Interpolation point
        sum1(:) = (/0.0,0.0,0.0,0.0,0.0/)
        sum2 = 0.0
      !  i = a
      !  j = b
        do ii = a-1,a+1
           do jj = b-1,b+1
              if(cc(ii,jj) .eq. 1)then
             !   print*,ii,jj,a,b,c_ib(a,b),ibpts(c_ib(a,b),:)
                d_par = -((ibpts(c_ib(a,b),1)-cell_c(ii,jj,1))*&
                            ib_norm(c_ib(a,b),1)+&
                          (ibpts(c_ib(a,b),2)-&
                          cell_c(ii,jj,2))*&
                          ib_norm(c_ib(a,b),2))
                d_hyp = sqrt((ibpts(c_ib(a,b),1)-cell_c(ii,jj,1))**2+(ibpts(c_ib(a,b),2)-cell_c(ii,jj,2))**2)
                d_per = sqrt(d_hyp**2-d_par**2)

                sum1(1) = sum1(1)+d_par/d_per
                sum2 = sum2+(1.0/d_per)    !Interpolation point

                sum1(2) = sum1(2)+u(ii,jj)/d_per   ! u at Interpolation

                sum1(3) = sum1(3)+v(ii,jj)/d_per   ! v at Interpolation
                
                sum1(4) = sum1(4)+p(ii,jj)/d_per   ! P at Interpolation

                sum1(5) = sum1(5)+rho(ii,jj)/d_per ! rho at Interpolation

              end if
           end do
        end do
        d_ip(i,j,:) = sum1(:)/sum2

     !Forcing the properties

     if(viscous .eq. 1)then
        u_ib = 0.0
        v_ib = 0.0 
           g = sqrt((ibpts(c_ib(i,j),1)-cell_c(i,j,1))**2+(ibpts(c_ib(i,j),2)-cell_c(i,j,2))**2)
           u(i,j) = u_ib+(g/d_ip(i,j,1))*(d_ip(i,j,2)-u_ib)
           v(i,j) = v_ib+(g/d_ip(i,j,1))*(d_ip(i,j,3)-v_ib)
           p(i,j) = d_ip(i,j,4)
           rho(i,j) = d_ip(i,j,5)
           T(i,j) = p(i,j)/(R*rho(i,j))
      else
           u(i,j) = d_ip(i,j,2)-d_ip(i,j,2)*ib_norm(c_ib(i,j),1)**2 
           v(i,j) = d_ip(i,j,3)-d_ip(i,j,3)*ib_norm(c_ib(i,j),2)**2 
           p(i,j) = d_ip(i,j,4)
           rho(i,j) = d_ip(i,j,5)
           T(i,j) = p(i,j)/(R*rho(i,j))    
      end if
     end subroutine forcing            
   end module immersed_routine
