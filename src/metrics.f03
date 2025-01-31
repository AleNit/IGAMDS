      
      MODULE metrics

      USE utilities
      USE refinement

      integer :: ne, ndof_e, ndof_tot
      integer, dimension(:), allocatable :: ndof

      real :: lg3, lg1, lg_con2, invdetgab, map
      real, dimension(3) :: g3, n, Bv, Gab, Gab_con 
      real, dimension(2,2) :: eg, eg_c
      real, dimension(2,3) :: e_t
      real, dimension(3,2)  :: G, e, G_con 
      real, dimension(3,3) :: H, Tm, H_t, T_g_e, Tm_2
      real, dimension(3,4) :: H2

      real, dimension(:), allocatable :: R
      real, dimension(:,:), allocatable :: dR, ddR

      CONTAINS

     !============================================================
     
     !It performs a metrics mapping on the refined geometry w.r.t.
     !the integration points and stores the metrics parameters in 
     !large rank arrays

      SUBROUTINE NURBS_mapping(pa)

      implicit none
     !----------------------------------------------------------
      type(multipatch), intent(inout) :: pa 
     !----------------------------------------------------------
      integer :: i, j, ku, kv, ru, rv, k, ki, kj, ii, jj
      real :: map, un, vn
      real, dimension(4,pa%pr+1) :: BSb_u 
      real, dimension(4,pa%qr+1) :: BSb_v
     !----------------------------------------------------------

      ru=pa%refu;  rv=pa%refv

     !allocate global mapping matrices
      allocate(pa%emap(ru,rv))
      allocate(pa%NC_u(ru,n_gpu))
      allocate(pa%NC_v(rv,n_gpv))
      allocate(pa%Gw(n_gpu,n_gpv))
     
     !allocate global storage NURBS basis functions matrices
      allocate(pa%R_s(ru,rv,n_gpu,n_gpv,ne))
      allocate(pa%dR_s(ru,rv,n_gpu,n_gpv,ne,2))
      allocate(pa%ddR_s(ru,rv,n_gpu,n_gpv,ne,3))

     !allocate storage matricx of basis functions needed for patchwise integration
      kj=sum(pa%ngpv_pi)
      ki=sum(pa%ngpu_pi)
      allocate(pa%R_s_pi(ki,kj,ne))
      allocate(pa%dR_s_pi(ki,kj,ne,2))
      allocate(pa%dA_pi(ki,kj))
      
      do j=pa%qr+1, size(pa%v)-pa%qr-1
        do i=pa%pr+1, size(pa%u)-pa%pr-1
          if (pa%u(i+1)/=pa%u(i)) then
          if (pa%v(j+1)/=pa%v(j)) then

           !element map (Jacobian)
            map=(pa%u(i+1)-pa%u(i))*(pa%v(j+1)-pa%v(j))/4.0
            pa%emap(i-pa%pr,j-pa%qr)=map


           !for standard Gauss-Legendre quadrature 
            do kv=1, n_gpv
              do ku=1, n_gpu

               !NURBS coordinates un,vn from Gauss element coordinates
                un=(pa%u(i+1)+pa%u(i)+gp_u(ku)*(pa%u(i+1)-pa%u(i)))/2.0
                vn=(pa%v(j+1)+pa%v(j)+gp_v(kv)*(pa%v(j+1)-pa%v(j)))/2.0

                call bspline_basis(i,pa%pr,un,pa%u,size(pa%u),BSb_u)
                call bspline_basis(j,pa%qr,vn,pa%v,size(pa%v),BSb_v)

                call NURBS_basis(i,j,BSb_u,BSb_v,pa)

               !fill NURBS derivative storage matrices
                pa%NC_u(i-pa%pr,ku)=un
                pa%NC_v(j-pa%qr,kv)=vn
                pa%R_s(i-pa%pr,j-pa%qr,ku,kv,:)=R
                pa%dR_s(i-pa%pr,j-pa%qr,ku,kv,:,:)=dR
                pa%ddR_s(i-pa%pr,j-pa%qr,ku,kv,:,:)=ddR

              enddo
            enddo


           !for reduced (patchwise) quadrature
            ii=i-pa%pr
            jj=j-pa%qr
            do kv=1, pa%ngpv_pi(jj) 
              do ku=1, pa%ngpu_pi(ii)

                ki=sum(pa%ngpu_pi(1:ii-1))+ku
                kj=sum(pa%ngpv_pi(1:jj-1))+kv

                un=pa%un_pi(ki)
                vn=pa%vn_pi(kj)
            
                call bspline_basis(i,pa%pr,un,pa%u,size(pa%u),BSb_u)
                call bspline_basis(j,pa%qr,vn,pa%v,size(pa%v),BSb_v)

                call NURBS_basis(i,j,BSb_u,BSb_v,pa)

               !fill NURBS derivative storage matrices
                pa%R_s_pi(ki,kj,:)=R
                pa%dR_s_pi(ki,kj,:,:)=dR

              enddo
            enddo


          endif
          endif
        enddo
      enddo

     !map Gauss weights
      do kv=1, n_gpv
        do ku=1, n_gpu
          pa%Gw(ku,kv)=wg_u(ku)*wg_v(kv)
        enddo
      enddo

      allocate(pa%lg3_s(2,ru,rv,n_gpu,n_gpv))
      allocate(pa%g3_s(2,ru,rv,n_gpu,n_gpv,3))
      allocate(pa%n_s(2,ru,rv,n_gpu,n_gpv,3))
      allocate(pa%Gab_s(2,ru,rv,n_gpu,n_gpv,3))
      allocate(pa%G_s(2,ru,rv,n_gpu,n_gpv,3,2))
      allocate(pa%Gc_s(2,ru,rv,n_gpu,n_gpv,3,2))

      ENDSUBROUTINE 

     !============================================================

     !It computes metrics parameters in configuration
     !passed as input argument
      
      SUBROUTINE metrics_objs(pa,co,cpp)

      implicit none
     !----------------------------------------------------------
      type(multipatch), intent(inout) :: pa
      integer, intent(in) :: co
      real, dimension(pa%cp_u,pa%cp_v,4), intent(in) :: cpp
     !----------------------------------------------------------
      integer :: i, j, k, a, b, c, ku, kv, ki,kj
     !----------------------------------------------------------
      
      do j=1, pa%refv
        do i=1, pa%refu
          if (pa%u(pa%pr+i+1)/=pa%u(pa%pr+i)) then
          if (pa%v(pa%qr+j+1)/=pa%v(pa%qr+j)) then
            

            do kv=1, n_gpv
              do ku=1, n_gpu
       
               !covariant base vectors matrix G and hessian matrix H
                G=0.0     ! G(1)=dxdu;    G(2)=dxdv

                k=0
                do c=0, pa%qr
                  do b=0, pa%pr
                    k=k+1
                    do a=1, 3
                      G(a,:)=pa%dR_s(i,j,ku,kv,k,:)*cpp(i+b,j+c,a)+G(a,:)
                    enddo
                  enddo
                enddo

               !basis vector g3
                CALL cross (G(:,1),G(:,2),g3)

               !length of g3 (= area dA)
                lg3=sqrt(dot_product(g3,g3))
               
               !normal vector n
                n=g3/lg3

               !covariant metric gab in Voigt notation
               !Gab(1)=Gab_11;Gab(2)=Gab_22;Gab(3)=Gab_12 
                Gab(1)=dot_product(G(:,1),G(:,1))
                Gab(2)=dot_product(G(:,2),G(:,2))
                Gab(3)=dot_product(G(:,1),G(:,2))
                
               !contravariant metric gab_con and base vectors g_con
                invdetgab=1.0/(Gab(1)*Gab(2)-Gab(3)*Gab(3))
                Gab_con(1)=invdetgab*Gab(2)
                Gab_con(3)=-invdetgab*Gab(3)
                Gab_con(2)=invdetgab*Gab(1)
                
                G_con(:,1)=G(:,1)*Gab_con(1)+G(:,2)*Gab_con(3)
                G_con(:,2)=G(:,1)*Gab_con(3)+G(:,2)*Gab_con(2)

               !fill metrics storage matrices
                pa%G_s(co,i,j,ku,kv,:,:)=G
                pa%Gc_s(co,i,j,ku,kv,:,:)=G_con
                pa%g3_s(co,i,j,ku,kv,:)=g3
                pa%Gab_s(co,i,j,ku,kv,:)=Gab
                pa%lg3_s(co,i,j,ku,kv)=lg3
                pa%n_s(co,i,j,ku,kv,:)=n
                
              enddo
            enddo
          

           !for reduced (patchwise) quadrature
            do kv=1, pa%ngpv_pi(j) 
              do ku=1, pa%ngpu_pi(i)

                ki=sum(pa%ngpu_pi(1:i-1))+ku
                kj=sum(pa%ngpv_pi(1:j-1))+kv

                G=0.0
                k=0
                do c=0, pa%qr
                  do b=0, pa%pr
                    k=k+1
                    do a=1, 3
                      G(a,:)=pa%dR_s_pi(ki,kj,k,:)*cpp(i+b,j+c,a)+G(a,:)
                    enddo
                  enddo
                enddo

               !basis vector g3
                CALL cross (G(:,1),G(:,2),g3)

               !length of g3 (= area dA)
                lg3=sqrt(dot_product(g3,g3))
           
               !fill NURBS derivative storage matrices
                pa%dA_pi(ki,kj)=lg3

              enddo
            enddo


          endif
          endif
        enddo
      enddo

      ENDSUBROUTINE 

     !============================================================

     !It computes metrics parameters on a specific location, out
     !of the integration points net
      
      SUBROUTINE spare_metrics(i,j,uu,vv,cpp,pa)

      implicit none
     !----------------------------------------------------------
      type(multipatch), intent(in) :: pa
      integer, intent(in) :: i, j
      real, intent(in) :: uu, vv
      real, dimension(pa%cp_u,pa%cp_v,4), intent(in)  :: cpp
     !----------------------------------------------------------
      integer :: c, b, k, a
      real, dimension(4,pa%pr+1) :: BSb_u
      real, dimension(4,pa%qr+1) :: BSb_v
     !----------------------------------------------------------


      CALL bspline_basis(i,pa%pr,uu,pa%u,size(pa%u),BSb_u)
      CALL bspline_basis(j,pa%qr,vv,pa%v,size(pa%v),BSb_v)

      CALL NURBS_basis(i,j,Bsb_u,Bsb_v,pa)

      G=0.0;  H=0.0;  H2=0.0    
      
      k=0
      do c=0, pa%qr
        do b=0, pa%pr
          k=k+1
          do a=1, 3
            G(a,:)=dR(k,:)*cpp(i-pa%pr+b,j-pa%qr+c,a)+G(a,:)
            H(a,:)=ddR(k,:)*cpp(i-pa%pr+b,j-pa%qr+c,a)+H(a,:)
          enddo
        enddo
      enddo

      CALL cross (G(:,1),G(:,2),g3)
      lg3=sqrt(g3(1)**2+g3(2)**2+g3(3)**2)
      
      n=g3/lg3
      
      H_t=transpose(H)
      Bv=matmul(H_t,n)
      
      Gab(1)=dot_product(G(:,1),G(:,1))
      Gab(2)=dot_product(G(:,2),G(:,2))
      Gab(3)=dot_product(G(:,1),G(:,2))

      invdetgab=1.0/(Gab(1)*Gab(2)-Gab(3)*Gab(3))
      Gab_con(1)=invdetgab*Gab(2)
      Gab_con(3)=-invdetgab*Gab(3)
      Gab_con(2)=invdetgab*Gab(1)

      G_con(:,1)=G(:,1)*Gab_con(1)+G(:,2)*Gab_con(3)
      G_con(:,2)=G(:,1)*Gab_con(3)+G(:,2)*Gab_con(2)

      lg1= sqrt(G(1,1)**2+G(2,1)**2+G(3,1)**2)
      e(:,1)=G(:,1)/lg1
      lg_con2= sqrt(G_con(1,2)**2+G_con(2,2)**2+G_con(3,2)**2)
      e(:,2)=G_con(:,2)/lg_con2

      e_t=transpose(e)
      eg=matmul(e_t,g_con)

     !trasformation matrix from G_con to E (for strain in voight notation) 
      Tm(1,1)=eg(1,1)*eg(1,1);  Tm(2,1)=eg(2,1)*eg(2,1)
      Tm(1,2)=eg(1,2)*eg(1,2);  Tm(2,2)=eg(2,2)*eg(2,2)
      Tm(3,1)=2*eg(1,1)*eg(2,1);  Tm(3,2)=2*eg(1,2)*eg(2,2)
      Tm(1,3)=2*eg(1,1)*eg(1,2);  Tm(2,3)=2*eg(2,1)*eg(2,2)
      Tm(3,3)=2*(eg(1,1)*eg(2,2)+eg(1,2)*eg(2,1))
     
     !trasformation matrix from E to G (for PK2 stress)
      Tm_2(1,1)=eg(1,1)*eg(1,1);    Tm_2(2,1)=eg(1,2)*eg(1,2)
      Tm_2(1,2)=eg(2,1)*eg(2,1);    Tm_2(2,2)=eg(2,2)*eg(2,2)
      Tm_2(3,1)=eg(1,1)*eg(1,2);    Tm_2(3,2)=eg(2,1)*eg(2,2)
      Tm_2(1,3)=2*eg(1,1)*eg(2,1);  Tm_2(2,3)=2*eg(1,2)*eg(2,2)
      Tm_2(3,3)=eg(1,1)*eg(2,2)+eg(1,2)*eg(2,1)

      eg_c=matmul(e_t,G)

     !trasformation matrix from G to e (for Chaucy stress)
      T_g_e(1,1)=  eg_c(1,1)*eg_c(1,1);  T_g_e(2,1)=  eg_c(2,1)*eg_c(2,1)
      T_g_e(1,2)=  eg_c(1,2)*eg_c(1,2);  T_g_e(2,2)=  eg_c(2,2)*eg_c(2,2)
      T_g_e(3,1)=  eg_c(1,1)*eg_c(2,1);  T_g_e(3,2)=  eg_c(1,2)*eg_c(2,2)
      T_g_e(1,3)=2*eg_c(1,1)*eg_c(1,2);  T_g_e(2,3)=2*eg_c(2,1)*eg_c(2,2)
      T_g_e(3,3)=eg_c(1,1)*eg_c(2,2)+eg_c(1,2)*eg_c(2,1)


      ENDSUBROUTINE 
     
     !============================================================ 

     !It evaluates BSpline basis functions, I, II and III derivatives
     !for a given location c, knot vector arr and polynomial order p
      
      SUBROUTINE bspline_basis (i,p,c,arr,dim_arr,Bsb)
      
      implicit none
     !----------------------------------------------------------
      integer, intent(in) :: i, p, dim_arr
      real, intent(in) :: c
      real, dimension(dim_arr), intent(in) :: arr
     !----------------------------------------------------------
      integer :: jj, r, k, pk, rk, j1, j2, s1, s2
      real :: saved, temp, d
      real, dimension(4,p+1) :: Bsb
      real, dimension(2,p+1) :: a
      real, dimension(p+1) :: left, right
      real, dimension(p+1,p+1) :: nd
     !----------------------------------------------------------

      Bsb=0.0;  nd=0.0;  a=0.0;  left=0.0;  right=0.0
    
      nd(1,1)=1.0
      do jj=1, p
        left(jj+1)=c-arr(i+1-jj)
        right(jj+1)=arr(i+jj)-c
        saved=0.0
        do r=0, jj-1
          nd(jj+1,r+1)=right(r+2)+left(jj-r+1)
          temp=nd(r+1,jj)/nd(jj+1,r+1)
          nd(r+1,jj+1)=saved+right(r+2)*temp
          saved=left(jj-r+1)*temp
        enddo
        nd(jj+1,jj+1)=saved
      enddo
      
     !load basis functions
      do jj=0, p;  Bsb(1,jj+1)=nd(jj+1,p+1);  enddo
    
     !compute derivatives
      do r=0, p             !loop over function index
        s1=0
        s2=1            !alternate rows in array a
        a(1,1)=1.0      !loop to compute 1st derivative
        do k=1, 3
          d=0.0
          rk=r-k
          pk=p-k
          if (r>=k) then
             a(s2+1,1)=a(s1+1,1)/nd(pk+2,rk+1)
             d=a(s2+1,1)*nd(rk+1,pk+1)
          endif
    
          if (rk>=-1) then;  j1=1;  else;  j1=-rk;  endif
          if ((r-1)<=pk) then;  j2=k-1;  else;  j2=p-r;  endif
    
          do jj=j1, j2
            a(s2+1,jj+1)=(a(s1+1,jj+1)-a(s1+1,jj))/nd(pk+2,rk+jj+1)
            d=d+a(s2+1,jj+1)*nd(rk+jj+1,pk+1)
          enddo
    
          if (r<=pk) then
            a(s2+1,k+1)=-a(s1+1,k)/nd(pk+2,r+1)
            d=d+a(s2+1,k+1)*nd(r+1,pk+1)
          endif
    
          Bsb(k+1,r+1)=d;  jj=s1;  s1=s2;  s2=jj     !switch rows
    
        enddo
      enddo
     
     !multiply through by the correct factors
      r=p
      do k=1, 3
        do jj=0, p
          Bsb(k+1,jj+1)=Bsb(k+1,jj+1)*r
        enddo
        r=r*(p-k)
      enddo
      
      ENDSUBROUTINE

     !============================================================

     !It returns NURBS basis function and I, II, III derivatives w.r.t
     !NURBS coordinates un and vn
           
      SUBROUTINE NURBS_basis(i,j,Bsb_u,Bsb_v,pa)

      implicit none
     !----------------------------------------------------------
      type(multipatch), intent(in) :: pa
      integer, intent(in) :: i, j
      real, dimension(4,pa%pr+1), intent(in) :: Bsb_u
      real, dimension(4,pa%qr+1), intent(in) :: Bsb_v
     !----------------------------------------------------------
      integer :: k, c, b
      real :: sumr
      real, dimension(2) :: dsumr
      real, dimension(3) :: ddsumr
      real, dimension(4) :: dddsumr
     !----------------------------------------------------------
     
      sumr=0.0;  dsumr=0.0;  ddsumr=0.0;  dddsumr=0.0

      k=0
      do c=0, pa%qr
        do b=0, pa%pr
          k=k+1
          R(k)=Bsb_u(1,b+1)*Bsb_v(1,c+1)*pa%cp(i-pa%pr+b,j-pa%qr+c,4)
          sumr=sumr+R(k)
                           
         !first derivatives
          dR(k,1)=Bsb_u(2,b+1)*Bsb_v(1,c+1)*pa%cp(i-pa%pr+b,j-pa%qr+c,4)
          dsumr(1)=dsumr(1)+dR(k,1)
          dR(k,2)=Bsb_u(1,b+1)*Bsb_v(2,c+1)*pa%cp(i-pa%pr+b,j-pa%qr+c,4)
          dsumr(2)=dsumr(2)+dR(k,2)
                            
         !second derivatives: 1-du^2, 2-dv^2, 3-dudv 
          ddR(k,1)=Bsb_u(3,b+1)*Bsb_v(1,c+1)*pa%cp(i-pa%pr+b,j-pa%qr+c,4)
          ddsumr(1)=ddsumr(1)+ddR(k,1)
          ddR(k,2)=Bsb_u(1,b+1)*Bsb_v(3,c+1)*pa%cp(i-pa%pr+b,j-pa%qr+c,4)
          ddsumr(2)=ddsumr(2)+ddR(k,2)
          ddR(k,3)=Bsb_u(2,b+1)*Bsb_v(2,c+1)*pa%cp(i-pa%pr+b,j-pa%qr+c,4)
          ddsumr(3)=ddsumr(3)+ddR(k,3)

        enddo
      enddo

     !divide through by sum
      do k=1, ne

        ddR(k,1)=ddR(k,1)/sumr-2*dR(k,1)*dsumr(1)/sumr**2-& 
                 R(k)*ddsumr(1)/sumr**2+2*R(k)*dsumr(1)**2/sumr**3
        ddR(k,2)=ddR(k,2)/sumr- 2*dR(k,2)*dsumr(2)/sumr**2-&
                 R(k)*ddsumr(2)/sumr**2+2*R(k)*dsumr(2)**2/sumr**3
        ddR(k,3)=ddR(k,3)/sumr-dR(k,1)*dsumr(2)/sumr**2-dR(k,2)*dsumr(1)/sumr**2-&
                 R(k)*ddsumr(3)/sumr**2+2*R(k)*dsumr(1)*dsumr(2)/sumr**3
        
        dR(k,1)=dR(k,1)/sumr-R(k)*dsumr(1)/sumr**2
        dR(k,2)=dR(k,2)/sumr-R(k)*dsumr(2)/sumr**2

        R(k)=R(k)/sumr

      enddo

      ENDSUBROUTINE 

     !============================================================

     !It evaluates the NURBS surface area in the passed configuration 
      
      SUBROUTINE garea(pa,np,co,area)

      implicit none
     !----------------------------------------------------------
      integer, intent(in) :: co, np
      type (multipatch), intent(in) :: pa
     !----------------------------------------------------------
      integer :: i, j, ku, kv
      real :: area
      real :: map, dA, gwl
     !----------------------------------------------------------

      if (np==1) then
        area=0.0
      endif


      do j=1, pa%refv
        do i=1, pa%refu

          map=pa%emap(i,j)

          do kv=1, n_gpv
            do ku=1, n_gpu

               dA=pa%lg3_s(co,i,j,ku,kv)
               gwl=pa%Gw(ku,kv)

               area=area+dA*map*gwl

            enddo
          enddo

        enddo
      enddo


      ENDSUBROUTINE

     !============================================================

     !It interpolates the transmembrane potential ant the recovery variable                                    
     !at the parametric location [uc,vc]

      SUBROUTINE get_point_eval(vv,ww,uc,vc,cpp,cpvv,cpww,pa)
             
      USE utilities, ONLY: findspan 
      
      implicit none
     !----------------------------------------------------------                                               
      type(multipatch), intent(in) :: pa
      real, intent(in) :: uc ,vc 
      real, dimension(pa%cp_u,pa%cp_v,4), intent(in) :: cpp
      real, dimension(pa%cp_u,pa%cp_v), intent(in) :: cpvv, cpww
      real, intent(out) :: vv, ww
     !----------------------------------------------------------                                               
      integer :: c, b, i, j, k  
      real :: sumnw
      real, dimension(4,pa%pr+1) :: Nu    
      real, dimension(4,pa%qr+1) :: Nv 
      real, dimension(ne) :: R
     !----------------------------------------------------------                                               
           
      i=findspan(uc,pa%u,pa%cp_u)  
      j=findspan(vc,pa%v,pa%cp_v) 
      
      call bspline_basis(i,pa%pr,uc,pa%u,size(pa%u),Nu)  
      call bspline_basis(j,pa%qr,vc,pa%v,size(pa%v),Nv) 
          
      SumNw=0.0;  k=0     
      do c=0, pa%qr
        do b=0, pa%pr     
          k=k+1          
          sumNw=Nu(1,b+1)*Nv(1,c+1)*cpp(i-pa%pr+b,j-pa%qr+c,4)+SumNw
        enddo
      enddo   
                
      k=0       
      do c=0, pa%qr
        do b=0, pa%pr
          k=k+1    
          R(k)=Nu(1,b+1)*Nv(1,c+1)*cpp(i-pa%pr+b,j-pa%qr+c,4)/SumNw
        enddo      
      enddo      
                  
      vv=0.0     
      ww=0.0  
                  
      k=0                   
      do c=0, pa%qr       
        do b=0, pa%pr
          k=k+1          
          vv=vv+R(k)*cpvv(i-pa%pr+b,j-pa%qr+c)  
          ww=ww+R(k)*cpww(i-pa%pr+b,j-pa%qr+c) 
        enddo       
      enddo           
                      
      ENDSUBROUTINE 

     !===================================================================


      ENDMODULE





