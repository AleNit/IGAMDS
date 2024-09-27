
!     ********************************************
!     ********************************************
!     **                                        **
!     **                                        **
!     **       IGA FINITE ELEMENTS SOLVER       **
!     **         FOR MONODOMAIN PROBLEM         **
!     **          OVER NURBS SURFACES           **
!     **                                        **
!     **                                        **
!     ********************************************
!     ********************************************
!
!
!--------------------------------------------------------------------
! Solve the dimensional monodomain quation (with ionic current model):
!
!      ∂v  
!  C_m -- - div( D * grad(v) ) + chi * i_ion = chi * i_a
!      ∂t  
!  
!   x        : [cm] space   
!   t        : [ms] time   
!   v        : [mV] action potential
!   C_m      : [mF/cm^3] tissue capacity
!   D        : [mA/(mV*cm)] conductivity tensor
!   chi      : [1/cm] surface to volume ratio
!   i_ion    : [mA/cm^2] ionic current density
!   i_a      : [mA/cm^2] applied current density
!  
!--------------------------------------------------------------------
!
!     Delevolep by Alessandro Nitti, Politechnic University of Bari
!     November 2019
!
!     Libreries needed:
!     HDF5, any release      
!     Lapack 3.6.0,  BLAS 3.6.0
!
!--------------------------------------------------------------------


      PROGRAM igamds
     
      USE utilities
      USE refinement
      USE metrics
      USE boundary_conditions
      USE std_mat_assembly
      USE postprocessing
      USE stimprot

      implicit none
     !----------------------------------------------------------
      integer :: i, j, k, l, c
      integer :: np, nss
      integer :: config
      integer :: n_steps, ns
      integer :: mmod, nmodc, kplt, nww, ncc
      integer :: iitr, eitr, nnz, itt, pint 

      logical :: dir_e, verb, scalev
      logical :: mp, isr

      real :: t_start, t_finish           
      real :: dt, time, tfin 
      real :: Cm, Diso, chi, vrest, vpeak, vamp
      real :: potch, iach, iich 
      real :: tol_abs, tol_rel, resid, trest
      real :: t1,t2
      
      real, dimension(50) :: modc

      real, dimension(:), allocatable :: RHS,LHS_CSR,RHSmat_CSR,mass_CSR,stiff_CSR
      integer, dimension(:), allocatable :: IA, JA

      character(50) :: newdir
      character(16) :: dir_g
      character(1) :: patch_ID

      type (multipatch), dimension(:), allocatable :: patch
     !---------------------------------------------------------- 

      
      CALL cpu_time(t_start)

      write(*,*)

      write(*,*) '*****************************************'
      write(*,*) '*                                       *'
      write(*,*) '*    ISOGEOMETRIC SOLVER FOR THE        *'
      write(*,*) '*       MONODOMAIN MODEL OVER           *'
      write(*,*) '*           NURBS SURFACES              *'
      write(*,*) '*                                       *'
      write(*,*) '*****************************************'



     !########################  READING INPUT DATA
      write(*,*)
      write(*,*) '... reading input data'  

     !read input parameters 
      CALL read_input_par(mp,isr,np,Cm,mmod,Diso,chi,dt,tfin,     &
                          pint,trest)

      
     !read membrane model data       
      CALL read_mem_par(mmod,nmodc,modc,nww,ncc,vrest,vpeak,scalev)

      vamp=vpeak-vrest

         
      if (.not.mp) then;  n_patches=1;  endif




     !########################  OUTPUT DIRECTORIES
      write(*,*)
      write(*,*) '... creating output directories'  

      dir_g='out'
      inquire(file=dir_g, exist=dir_e)
      if (.not.dir_e) then
        CALL system('mkdir '//dir_g)
      endif

      do i=1, n_patches
        write (patch_ID,'(I0)') i
        newdir=trim(dir_g)//'/p'//patch_ID//'_kin'
        inquire(file=trim(newdir), exist=dir_e)
        if (.not.dir_e) then
          CALL system('mkdir '//newdir)
        endif
      enddo




     !########################  READING GEOMETRY
      write(*,*)
      write(*,*) '... reading refined geometry'  

      allocate(patch(n_patches), ndof(n_patches))

      do i=1, n_patches
        CALL read_geometry_refined(patch(i),i)
        patch(i)%pr=patch(i)%p
        patch(i)%qr=patch(i)%q
        patch(i)%refu=patch(i)%cp_u-patch(i)%pr
        patch(i)%refv=patch(i)%cp_v-patch(i)%qr
      enddo

      ne=(patch(1)%pr+1)*(patch(1)%qr+1)    !basis function per element
      ndof_e=ne                             !degrees of freedom per element
     
     !compute the number of degrees of freedom and estimate the 
     !number of nonzero elements
      ndof_tot=0
      do i=1, n_patches           
        ndof(i)=patch(i)%cp_u*patch(i)%cp_v
        ndof_tot=ndof(i)+ndof_tot
      enddo



   
     !########################  NURBS MAPPING ON INTEGRATION POINTS
      n_gpu=patch(1)%pr+1; n_gpv=patch(1)%qr+1

      allocate(gp_u(n_gpu), wg_u(n_gpu), gp_v(n_gpv), wg_v(n_gpv))
      CALL gauss_legendre_pi_wi (n_gpu,gp_u,wg_u)
      CALL gauss_legendre_pi_wi (n_gpv,gp_v,wg_v)

     !allocate local NURBS bf matrices
      allocate(R(ne), dR(ne,2), ddR(ne,3))
      
      do i=1, n_patches
        CALL NURBS_mapping(patch(i))
      enddo



     
     !########################  REFERENCE CONFIGURATION METRICS MAPPING 
      write(*,*)
      write(*,*) '... computing reference metrics'  

      config=1    !1-reference config; 2-actual config
      do i=1, n_patches
        CALL metrics_objs(patch(i),config,patch(i)%cp) 
      enddo




     
     !########################  PRE-PROCESSING DATA STRUCTURE
      write(*,*)
      write(*,*) '... building processing data structure'  

     !the conductivity tensor is given in the patchwise parametric
     !space, i.e., in a 2D Cartesian frame
      DD=reshape([1.0,0.0,0.0,1.0],shape(DD))
      DD=Diso*DD

      
     !assign 1 dof to every cp and build dof matrix
      do l=1, n_patches
        allocate(patch(l)%dof(patch(l)%cp_u,patch(l)%cp_v))
        k=1
        do j=1, patch(l)%cp_v
          do i=1, patch(l)%cp_u
            patch(l)%dof(i,j)=k
            k=k+1
          enddo
        enddo
      enddo


     !allocate some processing arrays       
      n_steps=ceiling(tfin/dt)
      tfin=n_steps*dt

      CALL alloc(patch,mmod,nww,ncc)
      
      allocate(RHS(ndof_tot))





     !########################  MASS MATRIX
      write(*,*)
      write(*,*) '... assembling mass matrix'  


     !mass matrix assembly
      mass_mat_tot=0.0
      do i=1, n_patches
        CALL mass_mat_ass(patch(i),i)
      enddo

      k=0
      do i=1, n_patches
        mass_mat_tot(k+1:k+ndof(i),k+1:k+ndof(i))=patch(i)%mass_mat
        k=k+ndof(i)
      enddo

     !store the matrix in CSR format! 
      nnz=count(mass_mat_tot/=0.0)
  
      allocate(mass_CSR(nnz))
      allocate(IA(ndof_tot+1), JA(nnz))

      CALL CRSarr(ndof_tot,nnz,mass_mat_tot,mass_CSR,IA,JA)

      do i=1, n_patches
        deallocate(patch(i)%mass_mat)
      enddo
      deallocate(mass_mat_tot)





     !################################# RESTART      
      if (.not.isr) then

        write(*,*)
        write(*,*) '... reading initial solution'  
  
       !initial solution
        pot_tot_0=0.0
        vel_tot_0=0.0

        if (mmod==4) then
       !for Beeler-Reuter model 
        pot_tot_0=-83.3
        wrec_tot_0(:,1)=1.87e-7
        wrec_tot_0(:,2)=0.1644
        wrec_tot_0(:,3)=0.01
        wrec_tot_0(:,4)=0.9814
        wrec_tot_0(:,5)=0.9673
        wrec_tot_0(:,6)=0.0033
        wrec_tot_0(:,7)=0.9884
        endif

        Iion_tot_0=0.0
        Iapp_tot_0=0.0
        stiff_mat_tot=0.0
  
  
        k=0
        do c=1, n_patches
              patch(c)%pot_n=pot_tot_0(k+1:k+ndof(c))
              patch(c)%vel_n=vel_tot_0(k+1:k+ndof(c))
              patch(c)%Iion_n=Iion_tot_0(k+1:k+ndof(c))
              patch(c)%wrec_n=wrec_tot_0(k+1:k+ndof(c),:)
              patch(c)%Iion_n=Iion_tot_0(k+1:k+ndof(c))
              patch(c)%Iapp_n=Iapp_tot_0(k+1:k+ndof(c))
              k=k+ndof(c)
        enddo
  
       
        do i=1, n_patches
          
          write (patch_ID,'(I0)') i
          newdir=trim(dir_g)//'/p'//patch_ID//'_kin'
      
         !plot refined, undeformed geometry 
          CALL plot_geo(patch(i)%pr,patch(i)%qr,&
               patch(i)%u,patch(i)%v,patch(i)%cp_u,patch(i)%cp_v,&
               patch(i)%cp,newdir,i,patch(i))
        
        enddo
  
      
  
       !print output results
        do i=1, n_patches
          write (patch_ID,'(I0)') i
          newdir=trim(dir_g)//'/p'//patch_ID//'_kin'

         !apply affine transformations to get dimensional output
          if (scalev) then  
            patch(i)%pot_n=patch(i)%pot_n*vamp+vrest
            patch(i)%vel_n=patch(i)%vel_n*vamp
            patch(i)%wrec_n=patch(i)%wrec_n*vamp
            patch(i)%Iapp_n=patch(i)%Iapp_n*vamp
            patch(i)%Iion_n=patch(i)%Iion_n*vamp          
          endif
 
          CALL plot_time_n(ndof(i),patch(i)%pot_n,patch(i)%vel_n,   &
                         patch(i)%Iion_n,patch(i)%wrec_n,           &
                         patch(i)%Iapp_n,newdir,i,0,nww+ncc)  
        enddo
  
  
        CALL print_gq(0.0,0,maxval(pot_tot_0),maxval(vel_tot_0),     &
                      maxval(wrec_tot_0,1),maxval(Iion_tot_0),       &
                      maxval(Iapp_tot_0),dir_g,n_steps,0)
  
        pot_tot_n=pot_tot_0
        vel_tot_n=vel_tot_0
        wrec_tot_n=wrec_tot_0
        Iion_tot_n=Iion_tot_0
        Iion_tot_nm1=0.0
        Iapp_tot_n=Iapp_tot_0

      endif






     !########################  STIFFNESS MATRIX
      write(*,*)
      write(*,*) '... assemling stiffness matrix'  

     !asseble stiffness matrix
      config=1
      do i=1, n_patches
        patch(i)%stiff_mat=0.0
        CALL stiff_mat_ass(patch(i),i,config)
      enddo
  
      k=0
      do i=1, n_patches
        stiff_mat_tot(k+1:k+ndof(i),k+1:k+ndof(i))=patch(i)%stiff_mat
        k=k+ndof(i)
      enddo

      allocate(stiff_CSR(nnz))
      allocate(LHS_CSR(nnz))
      allocate(RHSmat_CSR(nnz))

      CALL CRSarr(ndof_tot,nnz,stiff_mat_tot,stiff_CSR,IA,JA)

      do i=1, n_patches
        deallocate(patch(i)%stiff_mat)
      enddo
      deallocate(stiff_mat_tot)





     !############################ LHS and RHS 
      write(*,*)
      write(*,*) '... generating LHS and static part of the RHS in CSR format'  

      LHS_CSR=Cm/dt*mass_CSR+0.5*stiff_CSR
      RHSmat_CSR=Cm/dt*mass_CSR-0.5*stiff_CSR

      deallocate(stiff_CSR,mass_CSR) 



     !assign tolerances and iteration limits for the mgmres solver
      iitr=10   !inner iteration
      eitr=5    !outer iteration
      tol_abs=1.0e-14
      tol_rel=1.0e-8


      


     !################################# RESTART
      if (isr) then

        write(*,*)
        write(*,*) '... reading restart files'

        CALL readrest(ndof_tot,dir_g,time,nss,     &
                  pot_tot_n,wrec_tot_n,Iion_tot_nm1)
      else
        time=0.0
        nss=1
      endif





     !################################# STIMULATION PROTOCOL
      write(*,*)
      write(*,*) '... reading stimulation protocol'

      CALL read_stimprot




      write(*,*)
      write(*,*) '>>> start time loop' 
      write(*,*) '    one step every',pint,'is displayed and written to file,' 
      write(*,*) '    the problem has',ndof_tot,'degrees of freedom'
      write(*,*)
   
     !print results at time 0 
      potch=maxval(pot_tot_0)


     


     !###################################### START TIME LOOP       
      kplt=1
      do ns=nss, n_steps

        CALL cpu_time(t1)

      
        Iapp_tot_n=0.0
        do i=1, n_patches;   patch(i)%Iapp_n=0.0;    enddo


        do i=1, stnum   !apply stimulation          
          if ( time>=tiss(i) .and. time<=tise(i) ) then
            
            if (potv(i)) then  !force potential value
              
              if (sside(i,1)==1) then                  
                k=1
                do j=1, patch(pID(i))%cp_v
                  pot_tot_n(k)=pot_prescr(i)
                  k=k+patch(pID(i))%cp_u
                enddo
              endif
              if (sside(i,3)==1) then                  
                k=patch(pID(i))%cp_u+1
                do j=1, patch(pID(i))%cp_v
                  pot_tot_n(k)=pot_prescr(i)
                  k=k+patch(pID(i))%cp_u
                enddo
              endif
              if (sside(i,2)==1) then                  
                k=(patch(pID(i))%cp_v-1)*patch(pID(i))%cp_u
                do j=1, patch(pID(i))%cp_u
                  pot_tot_n(k)=pot_prescr(i)
                  k=k+1
                enddo
              endif
              if (sside(i,4)==1) then                  
                k=1
                do j=1, patch(pID(i))%cp_u
                  pot_tot_n(k)=pot_prescr(i)
                  k=k+1
                enddo
              endif


            else  ! apply current

              selectcase (dstim(i))
                case (1)
                  CALL pointstim(Iappv(i),iua(i,1),iva(i,1),1,patch(pID(i)))
                case (2)
                  CALL edgestim(Iappv(i),iua(i,:),iva(i,:),1,patch(pID(i)))
                case (3)
                  CALL areastim(Iappv(i),iua(i,:),iva(i,:),1,patch(pID(i)))
              endselect    

              k=0
              do j=1, n_patches
                Iapp_tot_n(k+1:k+ndof(j))=patch(j)%Iapp_n
                k=k+ndof(j)
              enddo

            endif

          endif
        enddo

     
       !get cp matrices configuration at time n
        k=1
        do i=1, n_patches
          do j=1, patch(i)%cp_v
            do l=1, patch(i)%cp_u
              patch(i)%cp_pot_n(l,j)=pot_tot_n(k)
              patch(i)%cp_wrec_n(l,j,:)=wrec_tot_n(k,:)
              k=k+1
            enddo
          enddo
        enddo


       !integrate ionic current
        config=1 
        do i=1, n_patches
          patch(i)%Iion_n=0.0
          CALL ICSVI(config,patch(i),mmod,nmodc,modc,ncc,nww)
        enddo
  
        k=0
        do i=1, n_patches
          Iion_tot_n(k+1:k+ndof(i))=patch(i)%Iion_n
          k=k+ndof(i)
        enddo


       !assemble right hand side
        CALL amux(ndof_tot,pot_tot_n,RHS,RHSmat_CSR,JA,IA)
        RHS=RHS-0.5*chi*Iion_tot_nm1+1.5*chi*Iion_tot_n+chi*Iapp_tot_n

      
       !solve linear system (iterative solver)
        verb=.false.
        pot_tot_np1=pot_tot_n  !initial guess
        CALL pmgmres_ilu_cr(ndof_tot,nnz,IA,JA,LHS_CSR,pot_tot_np1,  &
               RHS,iitr,eitr,tol_abs,tol_rel,resid,itt,verb)
  
  
       !retrieve potential time derivative
        vel_tot_np1=(pot_tot_np1-pot_tot_n)/dt
  
  
       !solve the gating variable equation at each node
        do i=1, ndof_tot
          CALL solgw(mmod,modc,wrec_tot_np1(i,:),wrec_tot_n(i,:),   &
                          pot_tot_n(i),dt,time,nww+ncc) 
        enddo
  
  
        k=0
        do i=1, n_patches
          patch(i)%pot_n=pot_tot_np1(k+1:k+ndof(i))
          patch(i)%vel_n=vel_tot_np1(k+1:k+ndof(i))
          patch(i)%Iion_n=Iion_tot_n(k+1:k+ndof(i))
          patch(i)%wrec_n=wrec_tot_np1(k+1:k+ndof(i),:)
          k=k+ndof(i)
        enddo
  
        
       !print results on output files
        if (mod(ns,pint)==0) then
          do i=1, n_patches
            write (patch_ID,'(I0)') i
            newdir=trim(dir_g)//'/p'//patch_ID//'_kin'
    
           !apply affine transformations to get dimensional output
            if (scalev) then  
              patch(i)%pot_n=patch(i)%pot_n*vamp+vrest
              patch(i)%vel_n=patch(i)%vel_n*vamp
              patch(i)%wrec_n=patch(i)%wrec_n*vamp
              patch(i)%Iapp_n=patch(i)%Iapp_n*vamp
              patch(i)%Iion_n=patch(i)%Iion_n*vamp          
            endif
    
            CALL plot_time_n(ndof(i),patch(i)%pot_n,patch(i)%vel_n,   &
                           patch(i)%Iion_n,patch(i)%wrec_n,           &
                           patch(i)%Iapp_n,newdir,i,kplt,nww+ncc)  
          enddo
          kplt=kplt+1
        endif
  
  
       !write restart files
        if ( (ns>nss).and.(mod(time,trest).lt.dt) ) then
          write(*,*)
          write(*,*) ' ...writing restart'
          CALL writerest(ns,ndof_tot,time,pot_tot_np1,wrec_tot_np1,   &
                    Iion_tot_n,dir_g)
        endif
  
  
        CALL print_gq(time,ns,maxval(pot_tot_np1),maxval(vel_tot_np1),   &
                      maxval(wrec_tot_np1,1),maxval(Iion_tot_n),         &
                      maxval(Iapp_tot_n),dir_g,n_steps,nss)
  
       
       !screen output at time step n
        CALL cpu_time(t2)

        if (scalev) then
          potch=maxval(pot_tot_np1)*vamp+vrest
          iich=maxval(Iion_tot_n)
          iach=maxval(Iapp_tot_n)
        else
          potch=maxval(pot_tot_np1)
          iich=maxval(Iion_tot_n)
          iach=maxval(Iapp_tot_n)
        endif


       ! print step header
        if (mod(ns,pint)==0) then       

          write(*,*)
          write(*,"(1x,'step ',i8,', time = ',f9.5,' [ms]')") ns, time
          write(*,"(1x,'step exec. time = ',f9.5,' [s]')") t2-t1

          write(*,123) itt,resid
 123      format(1x,'mgmres total iterations: ',i5,', residual norm =',es0.4)        
  
          write(*,124) potch, iich, iach 
 124      format(1x,'maxval{pot,Iion,Iapp} = ',es0.4,' [mV], ',    &
                 es0.4,' [-], ',es0.4,' [-]')        
          write(*,*) '_________________________________________________'
     
        endif
       
       !update time variables for next iteration
        pot_tot_n=pot_tot_np1
        wrec_tot_n=wrec_tot_np1
        Iion_tot_nm1=Iion_tot_n
  
        time=time+dt
       


      enddo  ! TIME LOOP END


      CALL cpu_time(t_finish)
      
      write(*,*)
      write(*,"(' ... solution finalized, execution time: ',f16.8,' [min]')") (t_finish-t_start)/60.0
      write(*,*)
 

      ENDPROGRAM


