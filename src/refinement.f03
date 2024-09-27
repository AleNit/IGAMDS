      MODULE refinement
     
      USE utilities
    
      type multipatch
    
        integer :: p, q, qr, pr, du, dv, cp_u, cp_v, cp_u_r, cp_v_r
        integer :: refu, refv, insu, insv
        real, dimension(:), allocatable :: ur, vr
        real, dimension(:,:,:), allocatable :: cpr 
        real ,dimension(:), pointer :: u, v
        real, dimension(:,:,:), pointer :: cp
        
       !storage matrices
        real, dimension(:,:), allocatable :: emap
        real, dimension(:,:), allocatable :: NC_u, NC_v
        real, dimension(:,:), allocatable :: Gw
        
        real, dimension(:,:,:,:,:), allocatable :: R_s
        real, dimension(:,:,:,:,:,:), allocatable :: dR_s, ddR_s
        
        real, dimension(:,:,:,:,:), allocatable :: lg3_s
        real, dimension(:,:,:,:,:,:), allocatable :: g3_s, n_s, Gab_s
        real, dimension(:,:,:,:,:,:,:), allocatable :: G_s, Gc_s
       
        real, dimension(:,:), allocatable :: stiff_mat
        real, dimension(:,:), allocatable :: mass_mat
        real, dimension(:,:), allocatable :: dof
       
       !kinematic variables for dynamic analysis
        real, allocatable :: pot_n(:), vel_n(:),       &
                  Iion_n(:), i_s(:), wrec_n(:,:), Iapp_n(:), con_n(:,:)

        real, allocatable :: cp_pot_n(:,:), cp_wrec_n(:,:,:)
                   
      endtype
      
      CONTAINS

     !============================================================

     !It reads refined NURBS geometries from input file  
      
      SUBROUTINE read_geometry_refined(pa,k)

      implicit none
     !---------------------------------------------------------- 
      type(multipatch), intent(inout), target :: pa
      integer, intent(in) :: k
     !---------------------------------------------------------- 
      integer :: i, j
      character(1) :: patch_ID
      character(200) :: row, filename
     !---------------------------------------------------------- 

      write (patch_ID,'(I0)') k

      filename=trim('./input/igeo_p'//patch_ID//'.txt')
      open(99,file=filename,status="old",action="read")

      read(99,*) row
      read(99,*) pa%p
      read(99,*) pa%du
      allocate(pa%u(pa%du))
      read(99,*) pa%u

      read(99,*) row
      read(99,*) pa%q
      read(99,*) pa%dv
      allocate(pa%v(pa%dv))
      read(99,*) pa%v

      read(99,*) pa%cp_u, pa%cp_v
      allocate(pa%cp(pa%cp_u,pa%cp_v,4))
      do i=1, pa%cp_u
        do j=1, pa%cp_v
          read(99,*) pa%cp(i,j,:)
        enddo
      enddo

      close(99)

       
      ENDSUBROUTINE

     !============================================================

     !It performs some cheks on input geometry 
      
      SUBROUTINE check_geometry(pa,np)

      implicit none
     !----------------------------------------------------------
      integer, intent(in) :: np
      type(multipatch), intent(inout) :: pa
     !----------------------------------------------------------
      integer :: tmp, i, j, ii, jj
      real :: a
      real, dimension(3) :: tmp_array, diff_array
     !----------------------------------------------------------
 
      tmp=0

     !check input parameters matching
      if (pa%cp_u+pa%p+1/=pa%du) then
        write(*,*)
        write(*,*) 'WARNING: knot vector U, bf order p and cpts dont match'
        tmp=tmp+1
      endif

      if (pa%cp_v+pa%q+1/=pa%dv) then
        write(*,*)
        write(*,*) 'WARNING: knot vector V, bf order q and cpts dont match'
        tmp=tmp+1
      endif

     !check if there are coincidents control points
      do i=1, pa%cp_u
        do j=1, pa%cp_v

          tmp_array=pa%cp(i,j,1:3)

          do ii=1, pa%cp_u
            do jj=1, pa%cp_v

              if (i/=ii.and.j/=jj) then
                diff_array=tmp_array-pa%cp(ii,jj,1:3)
                a=dot_product(diff_array,diff_array)
                if (a==0.0) then

                  write(*,*) 'coincident control points '
                  write(*,*) pa%cp(i,j,:),' and ',pa%cp(ii,jj,:)
                  tmp=tmp+1

                endif
              endif

            enddo
          enddo

        enddo
      enddo

      if (tmp/=0) then
        write(*,*)
        write(*,*) 'PATCH', np
        write(*,*)
        stop
      endif


      ENDSUBROUTINE

     !============================================================
      

      ENDMODULE


                 
