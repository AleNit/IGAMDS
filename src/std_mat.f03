      
      MODULE std_mat_assembly

      USE refinement
      USE metrics
      USE bcs

      integer :: n_patches
      integer :: cp_u, cp_v, ndof_c, ndof_old
      integer :: enf
      logical :: u_edge

      real, dimension(2,2) :: DD

      integer, dimension(:), allocatable :: dof_loc
      real, dimension(:,:), allocatable :: stiff_mat_e, Ke, mass_mat_e
      real, dimension(:,:), allocatable :: stiff_mat_tot, mass_mat_tot

      real, allocatable :: pot_tot_0(:),     &
                           vel_tot_0(:),     &
                           Iion_tot_0(:),    &
                           Iapp_tot_0(:),    &
                           wrec_tot_0(:,:),  &
                           con_tot_0(:,:)
      
      real, allocatable :: pot_tot_n(:),     &
                           vel_tot_n(:),     &
                           Iion_tot_n(:),    &
                           Iapp_tot_n(:),    &
                           wrec_tot_n(:,:),  &
                           con_tot_n(:,:), i_s_tot(:)

      real, allocatable :: pot_tot_np1(:),      &
                           Iion_tot_np1(:),     &
                           wrec_tot_np1(:,:),   &
                           vel_tot_np1(:),      &
                           con_tot_np1(:,:)

      real, dimension(:), allocatable :: Iion_tot_nm1        


      CONTAINS

     !====================================================================
      
     !It assemble the linear mass matrix  
      
      SUBROUTINE mass_mat_ass(pa,np,VA,IA,JA,nd)

      implicit none
     !----------------------------------------------------------
      type(multipatch), intent(inout) :: pa
      integer, intent(in) :: np,nd
      real, intent(out) :: VA(nd)
      integer, dimension(nd), intent(out) :: IA,JA
     !----------------------------------------------------------
      integer :: i, j, k, cpi, cpj, kei, kej
      integer :: ii,jj,c,tmp(2),ind
      integer :: config
     !----------------------------------------------------------

      config=1
      c=0

      do j=1, pa%refv 
        if (pa%v(pa%qr+j+1)/=pa%v(pa%qr+j)) then
          do i=1, pa%refu
            if (pa%u(pa%pr+i+1)/=pa%u(pa%pr+i)) then
      
             mass_mat_e=0.0

              call e_mass_mat(i,j,pa,config)

             !insert element mass matrix into global [K] 
              dof_loc=0
              k=1
              do cpj=j, j+pa%qr
                do cpi=i, i+pa%pr
                  dof_loc(k)=pa%dof(cpi,cpj)
                  k=k+1
                enddo
              enddo

              do kei=1, ndof_e
                do kej=1, ndof_e
                        
!                 !full matrix assembly
!                  pa%mass_mat(dof_loc(kei),dof_loc(kej))=&
!                  mass_mat_e(kei,kej)+&
!                  pa%mass_mat(dof_loc(kei),dof_loc(kej))

                 !matrix assembly in COO format 
                  ii = dof_loc(kei)
                  jj = dof_loc(kej)
                  CALL indint(ind,c,nd,IA,JA,ii,jj)
                  IA(ind) = ii
                  JA(ind) = jj
                  VA(ind) = VA(ind) + mass_mat_e(kei,kej)

                enddo
              enddo

            endif
          enddo
        endif
      enddo


      ENDSUBROUTINE 

     !====================================================================

     !It assemble the linear stiffness matrix 
      
      SUBROUTINE stiff_mat_ass(pa,np,VA,IA,JA,nd)

      implicit none
     !----------------------------------------------------------
      type(multipatch), intent(inout) :: pa
      integer, intent(in) :: np, nd
      real, intent(out) :: VA(nd)
      integer, intent(out) :: IA(nd),JA(nd)      
     !----------------------------------------------------------
      integer :: i, j, k, cpi, cpj, kei, kej
      integer :: ii,jj,c,tmp(2),ind,co
      real, dimension(pa%cp_u*pa%cp_v,3) :: cplin
     !----------------------------------------------------------

      co=1
      c=0
      
      do j=1, pa%refv 
        if (pa%v(pa%qr+j+1)/=pa%v(pa%qr+j)) then
          do i=1, pa%refu
            if (pa%u(pa%pr+i+1)/=pa%u(pa%pr+i)) then

              stiff_mat_e=0.0

             !Curv coord implementation
              CALL e_stiff_mat(i,j,pa,co)
             
             !insert element stiffness matrix into global [K] 
              dof_loc=0;  k=1
              do cpj=j, j+pa%qr
                do cpi=i, i+pa%pr
                  dof_loc(k)=pa%dof(cpi,cpj)
                  k=k+1
                enddo
              enddo

              do kei=1, ndof_e
                do kej=1, ndof_e

!                 !full matrix assembly
!                  pa%stiff_mat(dof_loc(kei),dof_loc(kej))=&
!                  stiff_mat_e(kei,kej)+&
!                  pa%stiff_mat(dof_loc(kei),dof_loc(kej))

                 !matrix assembly in COO format 
                  ii = dof_loc(kei)
                  jj = dof_loc(kej)
                  CALL indint(ind,c,nd,IA,JA,ii,jj)
                  IA(ind) = ii
                  JA(ind) = jj
                  VA(ind) = VA(ind) + stiff_mat_e(kei,kej)

               enddo
              enddo

            endif
          enddo
        endif
      enddo


      ENDSUBROUTINE 

     !====================================================================  
      
     !It evaluate the element mass matrix 
      
      SUBROUTINE e_mass_mat(i,j,pa,co)

      implicit none
     !----------------------------------------------------------
      integer, intent(in):: i, j, co
      type(multipatch), intent(inout) :: pa
     !----------------------------------------------------------
      integer :: ku, kv, s, t, a,b,c,k
      integer :: ii, jj
      real :: dAe, gwl, un, vn, map
      real, dimension(ne) :: R
      real, dimension(ndof_e,ndof_e) :: R_mat, mea
      real, dimension(4,pa%pr+1) :: BSb_u 
      real, dimension(4,pa%qr+1) :: BSb_v      
     !----------------------------------------------------------

      map=pa%emap(i,j)

      do kv=1, n_gpv
        do ku=1, n_gpu 

         !recover useful metrics variables
          dAe=pa%lg3_s(co,i,j,ku,kv)
          
          R=pa%R_s(i,j,ku,kv,:)
          gwl=pa%Gw(ku,kv)

          do t=1, ndof_e
            do s=1, ndof_e
              R_mat(t,s)=R(s)*R(t)
            enddo
          enddo

          mea=R_mat*dAe
          mass_mat_e=mass_mat_e+mea*gwl*map

        enddo
      enddo

      ENDSUBROUTINE 

     !====================================================================

     !It evaluate the element stiffness matrix for a shell topology
      
      SUBROUTINE e_stiff_mat(i,j,pa,co)

      implicit none
     !----------------------------------------------------------
      integer, intent(in):: i, j, co
      type(multipatch), intent(in) :: pa
     !----------------------------------------------------------
      integer :: k, ku, kv, r, s
      integer :: al,be,ga,de
      real :: dAe,gwl,tmp,sc1,sc2
      real, dimension(3,2) :: G,Gc,eca
      real, dimension(2,2) :: DDcu
      real, dimension(ne,2) :: dR
     !----------------------------------------------------------


      eca=reshape([1.0,0.0,0.0,0.0,1.0,0.0],shape(eca))

      map=pa%emap(i,j)

      do kv=1, n_gpv
        do ku=1, n_gpu 
         
         !recover useful metrics variables
          dAe=pa%lg3_s(co,i,j,ku,kv)
          G=pa%G_s(co,i,j,ku,kv,:,:)
          Gc=pa%Gc_s(co,i,j,ku,kv,:,:)

          dR=pa%dR_s(i,j,ku,kv,:,:)
          gwl=pa%Gw(ku,kv)

         !bring the conductivity tensor in local curvilinear coordinates
         !DD is provided in the parametric space, i.e., in a 2D Cartesian frame
          do ga=1, 2
            do de=1, 2
              tmp=0.0 
              do al=1, 2 
                do be=1, 2

                  sc1=dot_product(Gc(:,ga),eca(:,al))
                  sc2=dot_product(Gc(:,de),eca(:,be))

                  tmp=tmp+DD(al,be)*sc1*sc2
                enddo
              enddo
              DDcu(ga,de)=tmp
            enddo
          enddo

          do r=1, ndof_e
            do s=1, ndof_e
              Ke(r,s)=dR(s,1)*dR(r,1)*DDcu(1,1)+   &
                      dR(s,1)*dR(r,2)*DDcu(1,2)+   &
                      dR(s,2)*dR(r,1)*DDcu(2,1)+   &
                      dR(s,2)*dR(r,2)*DDcu(2,2)
            enddo
          enddo

          stiff_mat_e=stiff_mat_e + Ke*dAe*map*gwl

        enddo
      enddo

      ENDSUBROUTINE 

     !====================================================================


      ENDMODULE



