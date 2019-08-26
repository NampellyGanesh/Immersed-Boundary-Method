!--------------------------------------------------------
!********Calculate Geometrical properties **************
!--------------------------------------------------------

  module geometrical_properties
  use variables

  implicit none
  ! i face in vertical and j  face in horizontal direction
  ! length of i faces
  integer             ::   n
  real                ::   magnitude
  real,dimension(3)   ::   z = [0,0,1]    
  real,dimension(2)   ::   diag1,diag2,j_res,i_res

  !--------------------------------------------------------------
  !                  Normals of domain
  !----------------------------------------------------------------
  contains 
    subroutine i_normal(domain,i_norm,i_length,ibpts,ib_norm)
           
          real,dimension(:,:,:)    ::  domain
          real,dimension(:,2:,:)   ::  i_norm
          real,dimension(:,2:)     ::  i_length
          real,dimension(:,:)      ::  ibpts,ib_norm
          real,dimension(:,:),allocatable      ::  res

          allocate(res(n_ibpts,2))

           open(10,file = file_name)
           n = imax*jmax
           do j = 1,jmax
              do i = 1,imax
                 read(10,*)domain(i,j,1:2)
              end do
           end do
           close(10)

           print*,'calculating geometric Quantities'     
           do j = 2,jmax
              do i = 1,imax
                 i_length(i,j) = sqrt((domain(i,j,1)-domain(i,j-1,1))**2+(domain(i,j,2)-domain(i,j-1,2))**2)  !length of I face

                 !normals of i face

                 i_res(1:2) = domain(i,j,1:2)-domain(i,j-1,1:2)
                 i_norm(i,j,1) = z(2)*0.0+z(3)*i_res(2)
                 i_norm(i,j,2) = -(z(1)*0.0+z(3)*i_res(1))
                 magnitude = sqrt(i_norm(i,j,1)**2+i_norm(i,j,2)**2)
                 i_norm(i,j,:) = i_norm(i,j,:)/magnitude
               end do 
           end do

           !writing into file 
           open(11,file = 'results/i_length.dat')
           write(11,*)'length, xnormal, ynormal, x,y'

           do j = 2,jmax
              do i = 1,imax
                 write(11,'(7e30.15)')i_length(i,j),i_norm(i,j,1)/50.0,i_norm(i,j,2)/50.0,domain(i,j,1),domain(i,j,2)
              end do
           end do
           close(11)
  !----------------------------------------------------------------------
  !                          Normals of Geometry
  !-----------------------------------------------------------------------

           open(15,file = 'ibsurfpts.dat')
           do i = 1,n_ibpts
              read(15,*)ibpts(i,1:2)
           end do
           close(15)

           !Finding Normals at intermediate points
           
          do i = 1,n_ibpts
             if(i .eq. 1)then
                res(i,:) = ibpts(i+1,:)-ibpts(n_ibpts,:)
             elseif(i .eq. n_ibpts)then
                res(i,:) = ibpts(1,:)-ibpts(i-1,:)
             else
                res(i,:) = ibpts(i+1,:)-ibpts(i-1,:)
             end if
          end do

          do i = 1,n_ibpts
             ib_norm(i,1) = z(2)*0.0+z(3)*res(i,2)
             ib_norm(i,2) = -(z(1)*0.0+z(3))*res(i,1)
             magnitude = sqrt(ib_norm(i,1)**2+ib_norm(i,2)**2)
             ib_norm(i,:) = ib_norm(i,:)/magnitude
          end do


          open(17,file = 'results/ibsurfpts')
          write(17,*)n_ibpts
          do i = 1,n_ibpts
             write(17,'(7e30.15)')ibpts(i,:),ib_norm(i,:)/10.0
          end do
          close(17)

          deallocate(res)

    !------------------------------------------------------------------------

     end subroutine i_normal

     subroutine j_normal(domain,j_norm,j_length)
         
           real,dimension(:,:,:)   ::   domain      
           real,dimension(2:,:,:)  ::   j_norm
           real,dimension(2:,:)    ::   j_length


           do j = 1,jmax
              do i = 2,imax
                 j_length(i,j) = sqrt((domain(i,j,1)-domain(i-1,j,1))**2+(domain(i,j,2)-domain(i-1,j,2))**2) !J Face length
 
                 j_res(1:2) = domain(i,j,1:2)-domain(i-1,j,1:2)
                 j_norm(i,j,1) = z(2)*0.0-z(3)*j_res(2)
                 j_norm(i,j,2) = -(z(1)*0.0-z(3)*j_res(1))
                 magnitude = sqrt(j_norm(i,j,1)**2+j_norm(i,j,2)**2)
                 j_norm(i,j,:) = j_norm(i,j,:)/magnitude
              end do
           end do

           !writing into file
           open(12,file = 'results/j_length.dat')
           write(12,*)'length,xnormal, ynormal,x,y'
 
           do j = 1,jmax
              do i = 2,imax
                 write(12,'(7e30.15)')j_length(i,j),j_norm(i,j,1)/50.0,j_norm(i,j,2)/50.0,domain(i,j,1),domain(i,j,2) 
              end do
           end do
           close(12)

     end subroutine j_normal

!---------------------------------------------------
!*******calculate volume of each cell*************
!--------------------------------------------------

subroutine volume_cell(domain,volume)

          real,dimension(:,:,:)    ::   domain
          real,dimension(2:,2:)    ::   volume
          real                     ::   cross_diag
          do j = 2,jmax
             do i = 2,imax
                diag1(1:2) = domain(i,j,1:2)-domain(i-1,j-1,1:2)
                diag2(1:2) = domain(i-1,j,1:2)-domain(i,j-1,1:2)
                cross_diag = diag1(1)*diag2(2)-diag1(2)*diag2(1)
                volume(i,j) = abs(cross_diag)/2.0
             end do
          end do
             
          open(13,file = 'results/Volume_cell.dat')
          write(13,*)'i,j,voume'
          do j = 2,jmax
             do i = 2,imax
                write(13,*)i,j,volume(i,j)
             end do
          end do

          close(13)

    end subroutine volume_cell
    
!---------------------------------------------------
!********Calculate cell_Centres********************
!---------------------------------------------------

subroutine cell_centre(domain,cell_c)
   
          real,dimension(:,:,:)      ::    domain
          real,dimension(:,:,:)      ::    cell_c
          do j = 2,jmax
             do i = 2,imax
                cell_c(i,j,:) = (domain(i,j,:)+domain(i,j-1,:)+domain(i-1,j-1,:)+domain(i-1,j,:))/4.0
             end do
          end do

          open(14,file = 'results/cell_centre.dat')
          write(14,*)'i,j,cell_centre'
          do j = 2,jmax
             do i = 2,imax
                 write(14,*),i,j,cell_c(i,j,1),cell_c(i,j,2)
             end do
          end do

         close(14)

   end subroutine cell_centre


end module geometrical_properties





