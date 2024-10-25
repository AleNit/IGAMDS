
      MODULE stimprot

      implicit none

      integer :: stnum
      integer, dimension(:), allocatable :: pID,dstim
      integer, dimension(:,:), allocatable :: sside 
      real, dimension(:), allocatable :: tiss,tise,Iappv,pot_prescr
      real, dimension(:,:), allocatable :: iua,iva 
      logical, dimension(:), allocatable :: potv

      CONTAINS

     !=====================================================
      
      SUBROUTINE read_stimprot

      implicit none
      integer :: i

      open(87,file='./input/stim_prot.in',status="old",action="read")
      read(87,*)
      read(87,*) stnum

      allocate(tiss(stnum),tise(stnum),Iappv(stnum),potv(stnum),pot_prescr(stnum))
      allocate(pID(stnum),dstim(stnum),sside(stnum,4))
      allocate(iua(stnum,2),iva(stnum,2))

      do i=1,stnum      
        read(87,*) 
        read(87,*)
        read(87,*) tiss(i)
        read(87,*) tise(i)
        read(87,*) pID(i)
        read(87,*) iua(i,1:2)
        read(87,*) iva(i,1:2)
        read(87,*) Iappv(i) 
        read(87,*) potv(i) 
        read(87,*) pot_prescr(i)
        read(87,*) sside(i,:)

        if ( (iua(i,1)==iua(i,2)) .and. (iva(i,1)==iva(i,2)) )  then  !point stimulus
          dstim(i)=1
        elseif ( (iua(i,1)==iua(i,2)) .or. (iva(i,1)==iva(i,2))  )  then  !edge stimulus
          dstim(i)=2
        else !area stimulus
          dstim(i)=3  
        endif  

      enddo

      close(87)
      

      ENDSUBROUTINE


      ENDMODULE

     !===================================================== 



     ! Allocates processing arrays and variables 
     ! defined in shared memory

      SUBROUTINE alloc(patch,m,nw,nc)

      USE refinement
      USE bcs
      USE std_mat_assembly

      implicit none
     !---------------------------------------------------------- 
      type(multipatch), dimension(n_patches) :: patch
      integer, intent(in) :: m, nw, nc
     !----------------------------------------------------------      
      integer :: i
     !---------------------------------------------------------- 


      do i=1, n_patches
        
!        allocate(patch(i)%mass_mat(ndof(i),ndof(i)))
!        allocate(patch(i)%stiff_mat(ndof(i),ndof(i)))
        
        allocate(patch(i)%pot_n(ndof(i)));      patch(i)%pot_n=0.0
        allocate(patch(i)%vel_n(ndof(i)));      patch(i)%vel_n=0.0      
        allocate(patch(i)%Iion_n(ndof(i)));     patch(i)%Iion_n=0.0
        allocate(patch(i)%wrec_n(ndof(i),nw+nc));     patch(i)%wrec_n=0.0
        allocate(patch(i)%Iapp_n(ndof(i)));     patch(i)%Iapp_n=0.0

        allocate(patch(i)%cp_pot_n(patch(i)%cp_u,patch(i)%cp_v))
        allocate(patch(i)%cp_wrec_n(patch(i)%cp_u,patch(i)%cp_v,nw+nc))

        patch(i)%cp_pot_n=0.0
        patch(i)%cp_wrec_n=0.0

      enddo
      
      allocate(stiff_mat_e(ndof_e,ndof_e), dof_loc(ndof_e))
      allocate(mass_mat_e(ndof_e,ndof_e))
      allocate(Ke(ndof_e,ndof_e))
      
!      allocate(stiff_mat_tot(ndof_tot,ndof_tot))
!      allocate(mass_mat_tot(ndof_tot,ndof_tot))
      
      allocate(pot_tot_0(ndof_tot))
      allocate(vel_tot_0(ndof_tot))
      allocate(Iion_tot_0(ndof_tot))
      allocate(Iapp_tot_0(ndof_tot))
      allocate(wrec_tot_0(ndof_tot,nw+nc))
      
      allocate(pot_tot_n(ndof_tot))
      allocate(vel_tot_n(ndof_tot))
      allocate(Iion_tot_n(ndof_tot))
      allocate(Iapp_tot_n(ndof_tot))
      allocate(wrec_tot_n(ndof_tot,nw+nc))

      allocate(pot_tot_np1(ndof_tot))
      allocate(vel_tot_np1(ndof_tot))
      allocate(Iion_tot_np1(ndof_tot))
      allocate(wrec_tot_np1(ndof_tot,nw+nc))

      allocate(Iion_tot_nm1(ndof_tot))


      ENDSUBROUTINE


