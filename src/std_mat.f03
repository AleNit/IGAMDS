      
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

      integer :: connl
      integer, dimension(:,:), allocatable :: CONN


      CONTAINS

     !====================================================================
      
     !It assemble the linear mass matrix  
      
      SUBROUTINE mass_mat_ass(pa,np,VA,IA,JA,nd,ndofm)

      implicit none
     !----------------------------------------------------------
      type(multipatch), intent(inout) :: pa
      integer, intent(in) :: np,nd,ndofm
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
                  ii = dof_loc(kei) + ndofm
                  jj = dof_loc(kej) + ndofm
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
      
      SUBROUTINE stiff_mat_ass(pa,ip,VA,IA,JA,nd,ndofm)

      implicit none
     !----------------------------------------------------------
      type(multipatch), intent(inout) :: pa
      integer, intent(in) :: ip,nd,ndofm
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

             !element stiffness matrix
              CALL e_stiff_mat(i,j,pa,ip,co)
             
             !insert element stiffness matrix into global matrix
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
                  ii = dof_loc(kei) + ndofm
                  jj = dof_loc(kej) + ndofm
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
      
      SUBROUTINE e_stiff_mat(i,j,pa,ip,co)

      implicit none
     !----------------------------------------------------------
      integer, intent(in):: i, j, co, ip
      type(multipatch), intent(in) :: pa
     !----------------------------------------------------------
      integer :: k, ku, kv, r, s
      integer :: al,be,ga,de
      real :: dAe,gwl,tmp,sc1,sc2
      real, dimension(3,2) :: G,Gc,eca
      real, dimension(2,2) :: DDcu,Gabcm
      real, dimension(3) :: e1,e2
      real, dimension(ne,2) :: dR
     !----------------------------------------------------------


     !conductivity coefficients w.r.t. the parametric space
      e1=[1.0,0.0,0.0]
      e2=[0.0,1.0,0.0]
      eca(:,1)=e1
      eca(:,2)=e2 


      map=pa%emap(i,j)

      do kv=1, n_gpv
        do ku=1, n_gpu 
         
         !recover useful metrics variables
          dAe=pa%lg3_s(co,i,j,ku,kv)
          G=pa%G_s(co,i,j,ku,kv,:,:)
          Gc=pa%Gc_s(co,i,j,ku,kv,:,:)

         !for fibers oriented along the mesh line directions
!          eca(:,1)=G(:,1)/norm2(G(:,1))
!          eca(:,2)=G(:,2)/norm2(G(:,2)) 

         !for Epicardium simulaiton 
!          CALL epi_fibers(co,i,j,pa,ip,ku,kv,e1,e2)
!          eca(:,1)=e1
!          eca(:,2)=e2


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

          stiff_mat_e = stiff_mat_e + Ke*dAe*map*gwl

        enddo
      enddo

      ENDSUBROUTINE 

     !====================================================================

      SUBROUTINE epi_fibers(co,i,j,pa,ip,ku,kv,a1,a2)

      implicit none
     !--------------------------------------------------------------------
      integer, intent(in) :: ip,i,j,ku,kv,co
      type(multipatch), intent(in) :: pa
      real, dimension(3), intent(out) :: a1,a2
     !--------------------------------------------------------------------
      real :: th,th2
      real, parameter :: alp=acos(-1.0)/3.0
      real, dimension(3) :: n,ta
      real, dimension(3,2) :: G,eca
     !--------------------------------------------------------------------
      
      G=pa%G_s(co,i,j,ku,kv,:,:)

      if ( ip==1 .or. ip==2 .or. ip==5 ) then 
          a1 = G(:,1)/norm2(G(:,1))
      elseif ( ip==3 .or. ip==4 ) then 
          a1 = -G(:,1)/norm2(G(:,1))
      endif
      a2 = g(:,2)/norm2(g(:,2))
      
     !find normal array
      CALL cross(a1,a2,n)
      n = n/sqrt(dot_product(n,n))
      
     !rotate the first base vecor in the tangent plane until its z component becomes almost null 
     !use the Rodrigues' rotation formula 
      if ( a1(3) < 0.0 ) then
        th=-acos(-1.0)/500.0
      else
        th=acos(-1.0)/500.0
      endif
      do while ( abs(a1(3)) > 1.0e-2 )
        CALL cross(n,a1,ta)  
        a1 = a1*cos(th) + ta*sin(th) + n*(dot_product(n,a1))*(1.0-cos(th))
      enddo
      
     !rotate again the base vector; rotation angle depend on the local surface orientation 
      th2=alp*abs(cos(n(3)))
      CALL cross(n,a1,ta)
      a1 = a1*cos(th2) + ta*sin(th2) + n*(dot_product(n,a1))*(1.0-cos(th2))
      
     !overwrite bitangent vector by a cross product
      CALL cross(a1,n,ta)
      a2 = ta/sqrt(dot_product(a2,a2))


      ENDSUBROUTINE

     !==================================================================== 

      ENDMODULE

