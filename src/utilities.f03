      
      MODULE utilities
     
      integer :: n_gpu, n_gpv
      real :: t_lsolve_i, t_lsolve_f
      real, dimension(:), allocatable :: gp_u, wg_u, gp_v, wg_v
      
      CONTAINS
     
     !============================================================
     
     !It returns knot span where u lies; rounding range must be 
     !considered for the last knot

      FUNCTION findspan (u,vec,n)
      
      implicit none
     !----------------------------------------------------------
      integer, intent(in) :: n
      real, intent(in) :: u
      real, dimension(n+1), intent(in) :: vec
     !---------------------------------------------------------- 
      integer :: i
      real :: eps
      integer :: findspan
     !----------------------------------------------------------
      
      eps=1e-10
      
      if (abs(u-vec(n+1))<eps) then  !spacial case; last knot
        findspan=n
        return
      endif
       
      do i=1, size(vec)-1
        if (u<vec(i+1)) then
          findspan=i
          return
        endif
      enddo
      
      return

      ENDFUNCTION
      
     !============================================================

     !It computes the cross product between two arrays[3] 
      
      SUBROUTINE cross(a,b,c)
   
      implicit none
     !---------------------------------------------------------- 
      real, dimension(3), intent(in)  :: a, b
      real, dimension(3), intent(out) :: c
     !----------------------------------------------------------
      
      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=a(3)*b(1)-a(1)*b(3)
      c(3)=a(1)*b(2)-a(2)*b(1)
       
      ENDSUBROUTINE
      
     !============================================================
      
     !It computes binomial coefficient of (cb_1,cb_2) 
      
      FUNCTION coeff_bin(cb_1,cb_2)
     
      implicit none
     !----------------------------------------------------------
      integer, intent(in) :: cb_1, cb_2
     !----------------------------------------------------------
      real :: coeff_bin
      integer :: cb_1_f, cb_2_f, diff_f, i
     !----------------------------------------------------------   
      
      cb_1_f=1;  cb_2_f=1;  diff_f=1
      
      do i=1, cb_1;  cb_1_f=cb_1_f*i;  enddo
      
      do i=1, cb_2;  cb_2_f=cb_2_f*i;  enddo
      
      do i=1, (cb_1-cb_2);  diff_f=diff_f*i;  enddo
      
      coeff_bin=cb_1_f/(diff_f*cb_2_f)
       
      ENDFUNCTION
     
     !============================================================

     !Given the number of Gauss points, it compute 2 arrays of coordinates
     !and weights for Gauss-Legendre numerical quadrature
           
      SUBROUTINE gauss_legendre_pi_wi(n_gp,xabsc,weig)
 
      implicit none
     !---------------------------------------------------------- 
      integer, intent(in) :: n_gp
     !----------------------------------------------------------
      integer :: i, j, m
      real ::  p1, p2, p3, pp, z, z1, eps, pi
      real, dimension(n_gp) :: xabsc, weig
     !----------------------------------------------------------
   
      eps=1.e-20 !precision of Newton-Raphson iterative method
      pi=acos(-1.0)
      m=int((n_gp+1)/2) !roots of Legendre polynomial are symmetric  
      
      do i = 1, m
        z=cos(pi*(i-0.25)/(n_gp+0.5))  !initial approximation 
        z1=0.0
        
        do while (abs(z-z1).gt.eps) !main loop of refinement    
          p1=1.0
          p2=0.0
          
          do j=1, n_gp
        !loop up the recurrence relation to get the 
        !Legendre polynomial evaluated at z
            p3=p2
            p2=p1
            p1=((2.d0*j-1.d0)*z*p2-(j-1.d0)*p3)/float(j)
          enddo
          
          pp=n_gp*(z*p1-p2)/(z*z-1.0) !Legendre polynomial derivative
          z1=z
          z=z1-p1/pp
        enddo
        
        xabsc(i)=-z
        xabsc(n_gp+1-i)=z                    !symmetry 
        weig(i)=2.0/((1.0-z*z)*pp*pp)
        weig(n_gp+1-i)=weig(i)               !symmetry
      enddo
        
      if (mod(n_gp,2).ne.0) then;  xabsc((n_gp+1)/2)=0.0;  endif
      
      ENDSUBROUTINE
     
     !============================================================
     
     !compute the norm to infinite of a matrix 
      
      SUBROUTINE norm_inf(A,n)

      implicit none
     !------------------------------------------------------------
      integer, intent(in) :: n
      real, dimension(n,n), intent(in) :: A
     !------------------------------------------------------------
      integer :: i, j
      real :: k, ni
      real, dimension(n) :: store
     !------------------------------------------------------------
      
      store=0.0

      do i=1, n
        k=0.0
        do j=1, n
          k=k+A(i,j)
        enddo
        store(i)=k
      enddo

      ni=maxval(store)

      ENDSUBROUTINE

     !============================================================ 

     !It computes an array multiplication: s=0.5*{b}^T[A]{b} 
      
      SUBROUTINE tri_prod(n,b,A,s)

      implicit none
     !------------------------------------------------------------
      integer, intent(in) :: n
      real, dimension(n), intent(in) :: b
      real, dimension(n,n), intent(in) :: A
      real, intent(out) :: s
     !------------------------------------------------------------
      integer :: i
      real, dimension(n) :: int_array
     !------------------------------------------------------------

      do i=1, n
         int_array(i)=dot_product(b,A(:,i))
      enddo
      s=0.5*dot_product(int_array,b)

      
      ENDSUBROUTINE

     !============================================================

      ENDMODULE


      
     !============================================================ 

     !It reads simulation parameters in input
           
      SUBROUTINE read_input_par(mp,is,np,Cm,mmod,       &
                         Diso,chi,dt,rt,pint,trest,redint)
      
      implicit none
     !----------------------------------------------------------
      logical, intent(out) :: mp, is, redint
      integer, intent(out) :: np, pint, mmod
      real, intent(out) :: Cm, Diso, chi, dt, rt
      real, intent(out) :: trest
     !----------------------------------------------------------
      integer :: i
     !----------------------------------------------------------

      
      open(87,file='./input/modelpar.in',status="old",action="read")

      read(87,*)
      read(87,*) np
      read(87,*) is
      read(87,*) pint
      read(87,*) trest

      do i=1, 3; read(87,*); enddo
      read(87,*) mmod  
      read(87,*) Cm
      read(87,*) Diso
      read(87,*) chi

      do i=1, 3; read(87,*); enddo
      read(87,*) dt 
      read(87,*) rt
      read(87,*) redint
   
      close(87)


      if (np==1) then
        mp=.false.
      else
        mp=.true.
        write(*,*)
        write(*,*) '... the multipatch option is temporary not available'
        write(*,*)
        stop
      endif


      ENDSUBROUTINE

     !============================================================

     !It reads parameters for the membrane model; this file is 
     !located in the code folder
           
      SUBROUTINE read_mem_par(mmod,nmodc,modc,nw,nc,vr,vp,scalev)

      implicit none
     !----------------------------------------------------------
      integer, intent(in) :: mmod
      integer, intent(out) :: nmodc, nw, nc
      logical, intent(out) :: scalev
      real, intent(out) :: vr,vp
      real, dimension(50), intent(inout) :: modc
     !----------------------------------------------------------
      integer :: i,lt
      character(512) :: exepath
     !----------------------------------------------------------

      CALL get_command_argument(0, exepath)
      lt=len_trim(exepath)

      modc=0.0

      selectcase(mmod)

        case(1) !Aliev-Panfilov model
          open(87,file=exepath(1:lt-5)//'membrane/Aliev_Panfilov.in',status="old",action="read")

        case(2) !Rogers-McCulloch model
          open(87,file=exepath(1:lt-5)//'membrane/Rogers-McCulloch.in',status="old",action="read")
       
        case(3) !Beeler-Reuter model
          open(87,file=exepath(1:lt-5)//'membrane/Beeler-Reuter.in',status="old",action="read")         

      endselect

     !read from file 
      read(87,*)
      read(87,*) nw,nc,nmodc
      do i=1, nmodc
        read(87,*) modc(i)
      enddo
      read(87,*) scalev
      read(87,*) vr
      read(87,*) vp

      close(87)


      ENDSUBROUTINE

     !============================================================


