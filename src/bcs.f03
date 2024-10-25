      
      MODULE bcs

      USE refinement
      USE utilities
      USE metrics

      integer :: n_con_tot
      integer, dimension(:), allocatable :: n_con, con_de

      CONTAINS


     !====================================================================

     !Integrate ionic current over the computational domain by a 
     !State Variable Interpolation (SVI) method

      SUBROUTINE ICSVI(co,pa,mmod,nmodc,modc,nc,nw)

      implicit none
     !----------------------------------------------------------
      integer, intent(in) :: co, mmod, nmodc, nc, nw
      real, dimension(50), intent(in) :: modc
      type(multipatch), intent(inout) :: pa
     !----------------------------------------------------------
      integer :: i, j, k, j1, j2, i1, i2, b, c, idir, ku, kv, n 
      integer :: ii, jj
      real :: dA, map, gwl, ww(nc+nw), vv
      real :: un, vn, Iint
      real, dimension(pa%pr+1,pa%qr+1) :: f_el
      real, dimension(:,:), allocatable :: f_n
     !----------------------------------------------------------


      allocate(f_n(pa%cp_u,pa%cp_v)); f_n=0.0

      do j=1, pa%refv    
        if (pa%v(pa%qr+j+1)/=pa%v(pa%qr+j)) then
          do i=1, pa%refu 
            if (pa%u(pa%pr+i+1)/=pa%u(pa%pr+i)) then

              map=pa%emap(i,j)

              do kv=1, n_gpv 
                do ku=1, n_gpu

                  dA=pa%lg3_s(co,i,j,ku,kv)
                  gwl=pa%Gw(ku,kv)

                 !Gauss points parametric coordinates
                  ii=pa%pr+i
                  jj=pa%qr+j
                  un=(pa%u(ii+1)+pa%u(ii)+gp_u(ku)*(pa%u(ii+1)-pa%u(ii)))/2.0
                  vn=(pa%v(jj+1)+pa%v(jj)+gp_v(kv)*(pa%v(jj+1)-pa%v(jj)))/2.0
                  
                 !state variables interpolation 
                  do n=1, nc+nw
                    CALL get_point_eval(vv,ww(n),un,vn,pa%cp,       &
                            pa%cp_pot_n,pa%cp_wrec_n(:,:,n),pa)
                  enddo

                  CALL Iioneval(vv,ww,Iint,mmod,nmodc,modc,nc,nw)

                  k=0
                  do c=1, pa%qr+1
                    do b=1, pa%pr+1
                      k=k+1
                      f_el(b,c)=Iint*map*gwl*pa%R_s(i,j,ku,kv,k)*dA 
                    enddo
                  enddo

                  f_n(i:i+pa%pr,j:j+pa%qr)=f_n(i:i+pa%pr,j:j+pa%qr)+f_el

                enddo
              enddo

            endif
          enddo
        endif
      enddo

     !store all nodal force components in a 1D array 
      k=1
      do j=1, pa%cp_v
        do i=1, pa%cp_u
          pa%Iion_n(k)=f_n(i,j)+pa%Iion_n(k)
          k=k+1
        enddo
      enddo

      deallocate(f_n)

      ENDSUBROUTINE 
      
     !====================================================================

     !Integrate ionic current over the computational domain by a 
     !State Variable Interpolation (SVI) method
     !A patchwise integration is adopted

      SUBROUTINE ICSVI_patchwise(co,pa,mmod,nmodc,modc,nc,nw)

      implicit none
     !----------------------------------------------------------
      integer, intent(in) :: co, mmod, nmodc, nc, nw
      real, dimension(50), intent(in) :: modc
      type(multipatch), intent(inout) :: pa
     !----------------------------------------------------------
      integer :: i, j, k, j1, j2, i1, i2, b, c, idir, ku, kv, n 
      integer :: ki, kj
      real :: dA, map, gwl, ww(nc+nw), vv
      real :: un, vn, Iint
      real, dimension(pa%pr+1,pa%qr+1) :: f_el
      real, dimension(:,:), allocatable :: f_n
      real, dimension(ne) :: BFV
     !----------------------------------------------------------


      allocate(f_n(pa%cp_u,pa%cp_v)); f_n=0.0


      do j=1, pa%refv    
        if (pa%v(pa%qr+j+1)/=pa%v(pa%qr+j)) then
          do i=1, pa%refu 
            if (pa%u(pa%pr+i+1)/=pa%u(pa%pr+i)) then

              map=pa%emap(i,j)

              do kv=1, ngpv_pi(j)
                do ku=1, ngpu_pi(i)

                  ki=sum(ngpu_pi(1:i-1))+ku
                  kj=sum(ngpv_pi(1:j-1))+kv
                  
                  dA=pa%dA_pi(ki,kj)
                  gwl=gwu_pi(ki)*gwv_pi(kj)
                  un=un_pi(ki)
                  vn=vn_pi(kj)
                 
                 !state variables interpolation 
                  do n=1, nc+nw
                    CALL get_point_eval(vv,ww(n),un,vn,pa%cp,       &
                            pa%cp_pot_n,pa%cp_wrec_n(:,:,n),pa)
                  enddo

                  CALL Iioneval(vv,ww,Iint,mmod,nmodc,modc,nc,nw)

                 !retrieve basis function 
                  BFV=pa%R_s_pi(ki,kj,:)

                  k=0
                  do c=1, pa%qr+1
                    do b=1, pa%pr+1
                      k=k+1
                      f_el(b,c)=Iint*gwl*BFV(k)*dA
                   enddo
                  enddo

                  f_n(i:i+pa%pr,j:j+pa%qr)=f_n(i:i+pa%pr,j:j+pa%qr)+f_el

                enddo
              enddo

            endif
          enddo
        endif
      enddo


     !store all nodal force components in a 1D array 
      k=1
      do j=1, pa%cp_v
        do i=1, pa%cp_u
          pa%Iion_n(k)=f_n(i,j)+pa%Iion_n(k)
          k=k+1
        enddo
      enddo
      

      deallocate(f_n)

      ENDSUBROUTINE 
     
     !=========================================================

      SUBROUTINE areastim(F,ub,vb,co,pa)

     ! returns the consistent nodal current distribution on a surface
     ! portion; ub,vb are load extension in parametric coord

      implicit none
     !----------------------------------------------------------
      integer, intent(in) :: co
      real, intent(in) :: F
      real, dimension(2), intent(in) :: ub, vb
      type(multipatch), intent(inout) :: pa
     !----------------------------------------------------------
      integer :: i, j, k, j1, j2, i1, i2, b, c, idir, ku, kv 
      real :: dA, map, gwl
      real, dimension(pa%pr+1,pa%qr+1) :: f_el
      real, dimension(:,:), allocatable :: f_n
     !----------------------------------------------------------

      allocate(f_n(pa%cp_u,pa%cp_v)); f_n=0.0

      i1=findspan(ub(1),pa%u,pa%cp_u)-pa%pr
      i2=findspan(ub(2),pa%u,pa%cp_u)-pa%pr
      j1=findspan(vb(1),pa%v,pa%cp_v)-pa%qr
      j2=findspan(vb(2),pa%v,pa%cp_v)-pa%qr

      do j=j1, j2      
        if (pa%v(pa%qr+j+1)/=pa%v(pa%qr+j)) then
          do i=i1, i2     
            if (pa%u(pa%pr+i+1)/=pa%u(pa%pr+i)) then

              map=pa%emap(i,j)

              do kv=1, n_gpv 
                do ku=1, n_gpu

                  dA=pa%lg3_s(co,i,j,ku,kv)
                  gwl=pa%Gw(ku,kv)

                  k=0
                  do c=1, pa%qr+1
                    do b=1, pa%pr+1
                      k=k+1
                      f_el(b,c)=F*map*gwl*pa%R_s(i,j,ku,kv,k)*dA 
                    enddo
                  enddo

                  f_n(i:i+pa%pr,j:j+pa%qr)=&
                        f_el(:,:)+f_n(i:i+pa%pr,j:j+pa%qr)

                enddo
              enddo

            endif
          enddo
        endif
      enddo

     !store all nodal force components in a 1D array 
      k=1
      do j=1, pa%cp_v
        do i=1, pa%cp_u
          pa%Iapp_n(k)=pa%Iapp_n(k)+f_n(i,j)
          k=k+1
        enddo
      enddo

      deallocate(f_n)

      ENDSUBROUTINE 
    
     !====================================================================

      SUBROUTINE areastim2(F,co,pa)

     ! returns the consistent nodal current distribution on a surface
     ! portion; ub,vb are load extension in parametric coord

      implicit none
     !----------------------------------------------------------
      integer, intent(in) :: co
      real, intent(in) :: F
      type(multipatch), intent(inout) :: pa
     !----------------------------------------------------------
      integer :: i, j, k, j1, j2, i1, i2, b, c, idir, ku, kv 
      real :: dA, map, gwl
      real :: ui, uf, vi, vf
      real, dimension(pa%pr+1,pa%qr+1) :: f_el
      real, dimension(:,:), allocatable :: f_n
     !----------------------------------------------------------

      allocate(f_n(pa%cp_u,pa%cp_v)); f_n=0.0

     !give load parametric coordinates
      ui=0.0
      vi=0.0
      vf=0.0999

      i1=findspan(ui,pa%u,pa%cp_u)-pa%pr
!      i2=findspan(uf,pa%u,pa%cp_u)-pa%pr
      j1=findspan(vi,pa%v,pa%cp_v)-pa%qr
      j2=findspan(vf,pa%v,pa%cp_v)-pa%qr

      do j=j1, j2      
        if (pa%v(pa%qr+j+1)/=pa%v(pa%qr+j)) then


          i2=j2-j+1


          do i=i1, i2     
            if (pa%u(pa%pr+i+1)/=pa%u(pa%pr+i)) then

              map=pa%emap(i,j)

              do kv=1, n_gpv 
                do ku=1, n_gpu

                  dA=pa%lg3_s(co,i,j,ku,kv)
                  gwl=pa%Gw(ku,kv)

                  k=0
                  do c=1, pa%qr+1
                    do b=1, pa%pr+1
                      k=k+1
                      f_el(b,c)=F*map*gwl*pa%R_s(i,j,ku,kv,k)*dA 
                    enddo
                  enddo

                  f_n(i:i+pa%pr,j:j+pa%qr)=&
                        f_el(:,:)+f_n(i:i+pa%pr,j:j+pa%qr)

                enddo
              enddo

            endif
          enddo
        endif
      enddo

     !store all nodal force components in a 1D array 
      k=1
      do j=1, pa%cp_v
        do i=1, pa%cp_u
          pa%Iapp_n(k)=pa%Iapp_n(k)+f_n(i,j)
          k=k+1
        enddo
      enddo

      deallocate(f_n)

      ENDSUBROUTINE 

     !====================================================================

      SUBROUTINE edgestim(F,ub,vb,co,pa)

     ! returns the consistent nodal current as a line load
     
      implicit none
     !----------------------------------------------------------
      integer, intent(in) :: co
      real, intent(in) :: F
      real, dimension(2), intent(in) :: ub, vb
      type(multipatch), intent(inout) :: pa
     !----------------------------------------------------------
      integer :: i1, i2, j1, j2, ku, kv, b, c, i, j, k
      integer :: ugauss, vgauss, xys, uv, uvf
      real :: un, vn, gwu, gwv, mapu, mapv
      real :: map, dsdu, gw, l_g, l_g_con
      real, dimension(:,:,:), allocatable :: ds
      real, dimension(:,:), allocatable :: f_n, f_el_l
     !----------------------------------------------------------

      allocate(f_n(pa%cp_u,pa%cp_v)); f_n=0.0
      allocate(f_el_l(pa%pr+1,pa%qr+1)); f_el_l=0.0
      allocate(ds(pa%pr+1,pa%qr+1,3)); ds=0.0

      i1=findspan(ub(1),pa%u,pa%cp_u)
      
      if (ub(1)==ub(2)) then
        un=ub(1)
        i2=i1
        ugauss=1
        gwu=1
        mapu=1
        uv=2
      else
        i2=findspan(ub(2),pa%u,pa%cp_u)
        if (ub(2)/=pa%u(size(pa%u))) then;  i2=i2-1;  endif
        ugauss=n_gpu
      endif

      j1=findspan(vb(1),pa%v,pa%cp_v)
      
      if (vb(1)==vb(2)) then
        vn=vb(1)
        j2=j1
        vgauss=1
        gwv=1
        mapv=1
        uv=1
      else
        j2=findspan(vb(2),pa%v,pa%cp_v)
        if (vb(2)/=pa%v(size(pa%v))) then;  j2=j2-1;  endif
        vgauss=n_gpv
      endif


      do j=j1, j2
        do i=i1, i2

          do kv=1, vgauss
            do ku=1, ugauss

              if (ub(1)/=ub(2)) then
                un=(pa%u(i+1)+pa%u(i)+gp_u(ku)*(pa%u(i+1)-pa%u(i)))/2.0
                mapu=(pa%u(i+1)-pa%u(i))/2.0
                gwu=wg_u(ku)
              endif
              if (vb(1)/=vb(2)) then
                vn=(pa%v(j+1)+pa%v(j)+gp_v(kv)*(pa%v(j+1)-pa%v(j)))/2.0
                mapv=(pa%v(j+1)-pa%v(j))/2.0
                gwv=wg_v(kv)
              endif

              map=mapu*mapv;  gw=gwu*gwv
              
             !----------- compute integration lenght ds
              CALL spare_metrics(i,j,un,vn,pa%cp,pa)

              dsdu=sqrt(dot_product(G(:,uv),G(:,uv)))
              ds=0.0
              k=0
              do c=0, pa%qr
                do b=0, pa%pr
                  k=k+1
                  ds(b+1,c+1,1)=R(k)*dsdu
                enddo
              enddo
             !-----------------------------------------
                  
              f_el_l(:,:)=F*ds(:,:,1)*gw*map
              
              f_n(i-pa%pr:i,j-pa%qr:j)=&
                     f_el_l(:,:)+f_n(i-pa%pr:i,j-pa%qr:j)

            enddo
          enddo
        enddo
      enddo

      k=1
      do j=1, pa%cp_v
        do i=1, pa%cp_u
          pa%Iapp_n(k)=pa%Iapp_n(k)+f_n(i,j)
          k=k+1
        enddo
      enddo

      deallocate(f_n,ds,f_el_l)

      ENDSUBROUTINE
      
     !====================================================================
      
      SUBROUTINE pointstim(F,ub,vb,co,pa)

      implicit none
     !----------------------------------------------------------
      integer, intent(in) :: co
      real, intent(in) :: F, ub, vb
      type(multipatch), intent(inout) :: pa
     !----------------------------------------------------------
      integer :: i, j, c, b, k
      real :: SumNw
      real, dimension(:,:), allocatable :: Bsb_u, Bsb_v
      real, dimension(:,:), allocatable :: f_n
     !----------------------------------------------------------
      
      i=findspan(ub,pa%u,pa%cp_u)
      j=findspan(vb,pa%v,pa%cp_v)

      allocate(f_n(pa%cp_u,pa%cp_v)); f_n=0.0

      allocate(Bsb_u(4,pa%pr+1), Bsb_v(4,pa%qr+1))
      call bspline_basis(i,pa%pr,ub,pa%u,size(pa%u),Bsb_u)
      call bspline_basis(j,pa%qr,vb,pa%v,size(pa%v),Bsb_v)

      SumNw=0.0
      do c=0, pa%qr
        do b=0, pa%pr
          SumNw=Bsb_u(1,b+1)*Bsb_v(1,c+1)*pa%cp(i-pa%pr+b,j-pa%qr+c,4)+SumNw
        enddo
      enddo

      do c=0, pa%qr
        do b=0, pa%pr
          f_n(i-pa%pr+b,j-pa%qr+c)=F*Bsb_u(1,b+1)*Bsb_v(1,c+1)*&
                 pa%cp(i-pa%pr+b,j-pa%qr+c,4)/SumNw
        enddo
      enddo

     !store all nodal force components in a 1D array 
      k=1
      do j=1, pa%cp_v
        do i=1, pa%cp_u
          pa%Iapp_n(k)=pa%Iapp_n(k)+f_n(i,j)
          k=k+1
        enddo
      enddo

      deallocate(f_n,Bsb_u,Bsb_v)

      ENDSUBROUTINE 
     
     !====================================================================

      SUBROUTINE read_patchwise_int

      implicit none
     !------------------------------------------------------------
      integer :: i,neu,nev,ngpu,ngpv
     !------------------------------------------------------------

      open(87,file='./input/qp_p1.in',status="old",action="read")    

      read(87,*) neu
      allocate(ngpu_pi(neu))
      do i=1,neu
        read(87,*) ngpu_pi(i)
      enddo
      read(87,*)
      read(87,*) ngpu
      allocate(un_pi(ngpu),gwu_pi(ngpu))
      do i=1,ngpu
        read(87,*) un_pi(i),gwu_pi(i)
      enddo
      read(87,*)


      read(87,*) nev
      allocate(ngpv_pi(nev))
      do i=1,nev
        read(87,*) ngpv_pi(i)
      enddo
      read(87,*)
      read(87,*) ngpv
      allocate(vn_pi(ngpv),gwv_pi(ngpv))
      do i=1,ngpv
        read(87,*) vn_pi(i),gwv_pi(i)
      enddo


      close(87)
 
      ENDSUBROUTINE

     !============================================================


      ENDMODULE

