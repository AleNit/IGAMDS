
     !solve the recovery variable equation by a 4th order
     !explicit Runge-Kutta method
      
      SUBROUTINE solgw(mmod,m,w,wo,v,dt,t,nwc)

      implicit none
     !---------------------------------------------------------- 
      integer, intent(in) :: mmod, nwc
      real, intent(in) :: m(50), dt, v, t
      real, intent(out) :: w(nwc), wo(nwc)
     !---------------------------------------------------------- 
      real, dimension(nwc) :: k1, k2, k3, k4, tmp1, tmp2, tmp3
     !---------------------------------------------------------- 

 
     !RK4 integration steps 
      CALL wfun(wo,v,t,mmod,m,k1,nwc)

      tmp1=wo+k1*dt/2.0
      CALL wfun(tmp1,v,t+dt/2.0,mmod,m,k2,nwc)

      tmp2=wo+k2*dt/2.0
      CALL wfun(tmp2,v,t+dt/2.0,mmod,m,k3,nwc)

      tmp3=wo+k3*dt
      CALL wfun(tmp3,v,t+dt,mmod,m,k4,nwc) 

      w=wo+dt*(k1/6.0+k2/3.0+k3/3.0+k4/6.0)


      ENDSUBROUTINE

     !====================================================================

      SUBROUTINE wfun(wo,v,t,mmod,m,k,nw)

      implicit none
     !----------------------------------------------------------
      integer, intent(in) :: mmod, nw
      real, intent(in) :: m(50), t, wo(nw), v
      real, intent(out) :: k(nw)
     !---------------------------------------------------------- 
      real :: alpha_x1, alpha_m, alpha_h, alpha_j, alpha_d, alpha_f
      real :: beta_x1, beta_m, beta_h, beta_j, beta_d, beta_f
      real :: i_s, V_Ca
     !---------------------------------------------------------- 
          
      selectcase(mmod)

        case(1)   !Aliev-Panfilov model
        k(1)=(m(6)+m(4)*wo(1)/(m(5)+v))*(-wo(1)-m(1)*v*(v-m(3)-1.0))

        case(2)   !Mitchell-Schaeffer model
        if ((v-m(5))<=1.0e-12) then
          k(1)=(1.0-wo(1))/m(3)
        else
          k(1)=-wo(1)/m(4)
        endif

        case(3)   !Rogers-McCulloch model
        k(1)=m(5)*(v-m(3)*wo(1))

        case(4)   !Beeler-Reuter model
        alpha_x1=0.0005*exp(0.083*(v+50.0))/(exp(0.057*(v+50.0))+1.0)
        beta_x1=0.0013*exp(-0.06*(v+20.0))/(exp(-0.04*(v+20.0))+1.0)

        alpha_m=-(v+47.0)/(exp(-0.1*(v+47.0))-1.0)
        beta_m=40.0*exp(-0.056*(v+72.0))

        alpha_h=0.126*exp(-0.25*(v+77.0))
        beta_h=1.7/(exp(-0.082*(v+22.5))+1.0)

        alpha_j=0.055*exp(-0.25*(V+78.0))/(exp(-0.2*(V+78.0))+1.0)
        beta_j=0.3/(exp(-0.1*(V+32.0))+1.0)

        alpha_d=0.095*exp(-0.01*(v-5.0))/(exp(-0.072*(v-5.0))+1.0)
        beta_d=0.07*exp(-0.017*(v+44.0))/(exp(0.05*(v+44.0))+1.0)

        alpha_f=0.012*exp(-0.008*(v+28.0))/(exp(0.15*(v+28.0))+1.0)
        beta_f=0.0065*exp(-0.02*(v+30.0))/(exp(-0.2*(V+30.0))+1.0)

       !concentration 
        V_Ca=-82.3-13.0287*log(wo(1)) 
        i_s=m(9)*wo(6)*wo(7)*(v-V_Ca)
        
        k(1)=-m(1)*i_s+m(3)*(m(2)-wo(1))
        
       !rec vars 
        k(2)=alpha_x1*(1.0-wo(2))-beta_x1*wo(2)
        k(3)=alpha_m*(1.0-wo(3))-beta_m*wo(3)
        k(4)=alpha_h*(1.0-wo(4))-beta_h*wo(4)
        k(5)=alpha_j*(1.0-wo(5))-beta_j*wo(5)
        k(6)=alpha_d*(1.0-wo(6))-beta_d*wo(6)
        k(7)=alpha_f*(1.0-wo(7))-beta_f*wo(7)

      endselect


      ENDSUBROUTINE

     !====================================================================

      SUBROUTINE Iioneval(v,w,Iint,mmod,nmodc,modc,nw,nc)

      implicit none
     !---------------------------------------------------------- 
      integer, intent(in) :: mmod, nmodc, nw, nc
      real, intent(in) :: v, w(nw+nc), modc(50)
      real, intent(out) :: Iint
     !---------------------------------------------------------- 
      real :: V_Ca, i_k1, i_Na, i_x1, i_s
     !---------------------------------------------------------- 

      i_s=0.0

      selectcase(mmod)

        case(1)   !Aliev-Panfilov model
          Iint=modc(1)*v*(v-modc(2))*(1.0-v)-v*w(1)

        case(2)   !Mitchell-Shaeffer model
          Iint=v/modc(2)+w(1)/modc(1)*v**2*(v-1.0)

        case(3)   !Rogers-McCulloch model  
          Iint=modc(1)*v*(v-modc(4))*(1.0-v)-modc(2)*v*w(1)

        case(4)   !Beeler-Reuter model
          
          V_Ca=-82.3-13.0287*log(w(1))

          i_K1=modc(4)*(4.0*(exp(0.04*(v+85.0))-1.0)/    &
               (exp(0.08*(v+53.0))+exp(0.04*(v+53.0)))+   &
               0.2*(v+23.0)/(1.0-exp(-0.04*(v+23.0))))
          i_x1=modc(5)*w(2)*(exp(0.04*(v+77.0))-1.0)/exp(0.04*(v+35.0))
          i_Na=(modc(7)*w(3)**3*w(4)*w(5)+modc(8))*(v-modc(6))
          i_s=modc(9)*w(6)*w(7)*(v-V_Ca)

          Iint=-(i_k1+i_x1+i_Na+i_s)


      endselect


      ENDSUBROUTINE

     !==================================================================== 


