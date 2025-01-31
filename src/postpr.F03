      
      MODULE postprocessing

      USE refinement
      USE metrics

      integer :: stpr
      logical :: pstress

      CONTAINS

     !====================================================================

     !It prints reference or actual geometry on output file 
      
      SUBROUTINE plot_geo(p,q,u,v,cp_u,cp_v,cpp,dir,np,pa)

      implicit none
     !----------------------------------------------------------
      integer, intent(in) :: p, q, cp_u, cp_v, np
      real, dimension(cp_u+p+1), intent(in) :: u
      real, dimension(cp_v+q+1), intent(in) :: v
      real, dimension(cp_u,cp_v,4), intent(in) :: cpp
      character(*), intent(in) :: dir
      type(multipatch), intent(in) :: pa
     !----------------------------------------------------------
      integer :: i, j, k
      character(100) :: filename
      character(1) :: p_ID
     !----------------------------------------------------------
    
      write (p_ID,'(I0)') np

      filename=trim(dir)//'/refgeo.txt'
      
      open(97, file=filename, action="write")

      write(97,*) p;  write(97,*) q
      write(97,*) size(u);  write(97,*) size(v)
      write(97,*) cp_u;  write(97,*) cp_v
      do i=1, size(u);  write(97,*) u(i);  enddo
      do i=1, size(v);  write(97,*) v(i);  enddo
      
      do i=1, 4;  do j=1, cp_v;  do k=1, cp_u
          write(97,*) cpp(k,j,i)
      enddo;  enddo;  enddo

      close(97)

    
      ENDSUBROUTINE 
     

     !===================================================================

      ENDMODULE
     
     
     
     !===================================================================

     !It prints the potential, velocity to output files
      
      SUBROUTINE plot_time_n(ndof,pot,vel,Iion,wrec,Iapp,dir,np,ns,nw)

      USE hdf5

      implicit none
     !----------------------------------------------------------
      integer, intent(in) :: ndof, np, ns, nw
      real, dimension(ndof), intent(in) :: pot,vel,Iion,Iapp
      real, dimension(ndof,nw), intent(in) :: wrec
      character(*), intent(in) :: dir
     !----------------------------------------------------------
      integer :: i 
      integer :: hdferr, nd, ndw
      integer(HSIZE_T) :: dimu(1), dimw(2)
      integer(HID_T) :: filespace, memspace
      integer(HID_T) :: dpot, dvel, diion, diapp, dwrec
      character(100) :: filename
      character(6) :: time
      character(1) :: p_ID
     !----------------------------------------------------------

      write (p_ID,'(I0)') np
      write (time,'(I6.6)') ns

       
      nd=1
      dimu(1)=ndof
      ndw=2
      dimw=[ndof,nw]

      filename=trim(dir)//'/t'//trim(time)//'.h5'
      
      CALL h5open_f(hdferr)

      CALL h5fcreate_f(filename,H5F_ACC_TRUNC_F,filespace,hdferr)
      CALL h5screate_simple_f(nd,dimu,memspace, hdferr)  
     
     !--------------- 
      CALL h5dcreate_f(filespace,'pot',H5T_NATIVE_DOUBLE,  &
                       memspace,dpot,hdferr)
      CALL h5dwrite_f(dpot,H5T_NATIVE_DOUBLE,pot,   &
                        dimu,hdferr)
      CALL h5dclose_f(dpot, hdferr)

     !---------------    
      CALL h5dcreate_f(filespace,'iion',H5T_NATIVE_DOUBLE,  &
                       memspace,diion,hdferr)
      CALL h5dwrite_f(diion,H5T_NATIVE_DOUBLE,iion,   &
                        dimu,hdferr)
      CALL h5dclose_f(diion, hdferr) 

      !---------------   
      CALL h5dcreate_f(filespace,'iapp',H5T_NATIVE_DOUBLE,  &
                       memspace,diapp,hdferr)
      CALL h5dwrite_f(diapp,H5T_NATIVE_DOUBLE,iapp,   &
                        dimu,hdferr)
      CALL h5dclose_f(diapp, hdferr) 

      CALL h5sclose_f(memspace, hdferr)       

     !---------------      
      CALL h5screate_simple_f(ndw,dimw,memspace, hdferr)  
      CALL h5dcreate_f(filespace,'wrec',H5T_NATIVE_DOUBLE,  &
                       memspace,dwrec,hdferr)
      CALL h5dwrite_f(dwrec,H5T_NATIVE_DOUBLE,wrec,   &
                        dimw,hdferr)
      CALL h5dclose_f(dwrec, hdferr) 

     !--------------- 
      CALL h5sclose_f(memspace, hdferr)        
      CALL h5fclose_f(filespace, hdferr)        


      ENDSUBROUTINE 

     !===================================================================

     !It prints global quantities on putput file 
            
      SUBROUTINE print_gq(t,ns,pot,vel,wrec,iion,iapp,dir,n_steps,nss,itt,resid)

      implicit none
     !----------------------------------------------------------
      integer, intent(in) :: ns, n_steps, nss, itt
      real, intent(in) :: pot, vel, t, wrec, iion, iapp, resid
      character(*), intent(in) :: dir
     !----------------------------------------------------------
      character(100) :: filename
     !----------------------------------------------------------

      if (ns==nss) then
        filename=trim(dir)//'/globq.out'
        open(193, file=filename, action="write")
      endif

      write(193,128) t, ns, pot, wrec, iion, iapp, itt, resid
 128  format(1x,e16.9,1x,i5,1x,e16.9,1x,e16.9,1x,e16.9,1x,e16.9,1x,i3,1x,e16.9)     

      if (ns==n_steps)  close(193)

      ENDSUBROUTINE 
     
     !===================================================================

     !writes restart files
            
      SUBROUTINE writerest(ns,ndof,t,pot,wrec,Iionnm1,dir)

      USE hdf5

      implicit none
     !----------------------------------------------------------
      integer, intent(in) :: ns, ndof
      real, intent(in) :: t, pot(ndof), wrec(ndof), Iionnm1(ndof)
      character(*), intent(in) :: dir
     !----------------------------------------------------------
      integer :: i
      integer :: hdferr, nd
      integer(HSIZE_T) :: dimu(1)
      integer(HID_T) :: filespace, memspace
      integer(HID_T) :: dpot, diion, dwrec
      character(100) :: filename
     !----------------------------------------------------------

      filename=trim(dir)//'/restgq.out'
      open(93, file=filename, action="write", status='unknown')

      write(93,*) t
      write(93,*) ns
      write(93,*) ndof
      
      close(93)
      
       
      nd=1
      dimu(1)=ndof

      filename=trim(dir)//'/restvar.h5'
      
      CALL h5open_f(hdferr)

      CALL h5fcreate_f(filename,H5F_ACC_TRUNC_F,filespace,hdferr)
     !---------------
      CALL h5screate_simple_f(nd,dimu,memspace, hdferr)  
     !--------------- 
      
      CALL h5dcreate_f(filespace,'pot',H5T_NATIVE_DOUBLE,  &
                       memspace,dpot,hdferr)
      CALL h5dwrite_f(dpot,H5T_NATIVE_DOUBLE,pot,   &
                        dimu,hdferr)
      CALL h5dclose_f(dpot, hdferr)

      CALL h5dcreate_f(filespace,'iionnm1',H5T_NATIVE_DOUBLE,  &
                       memspace,diion,hdferr)
      CALL h5dwrite_f(diion,H5T_NATIVE_DOUBLE,Iionnm1,   &
                        dimu,hdferr)
      CALL h5dclose_f(diion, hdferr) 

      CALL h5dcreate_f(filespace,'wrec',H5T_NATIVE_DOUBLE,  &
                       memspace,dwrec,hdferr)
      CALL h5dwrite_f(dwrec,H5T_NATIVE_DOUBLE,wrec,   &
                        dimu,hdferr)
      CALL h5dclose_f(dwrec, hdferr) 
     
     !--------------- 
      CALL h5sclose_f(memspace, hdferr)        
     !--------------- 
      CALL h5fclose_f(filespace, hdferr)        


      ENDSUBROUTINE 

     !===================================================================

     !writes restart files
            
      SUBROUTINE readrest(na,dir,t,ns,pot,wrec,Iionnm1)

      USE hdf5

      implicit none
     !----------------------------------------------------------
      integer, intent(in) :: na
      character(*), intent(in) :: dir
      integer, intent(out) :: ns
      real, intent(out) :: pot(na), wrec(na), Iionnm1(na)
      real, intent(out) :: t
     !----------------------------------------------------------
      integer :: i, ndof, k
      integer :: hdferr, nd
      logical :: dire1, dire2
      integer(HSIZE_T) :: dims(1)
      integer(HID_T) :: file_ID, memspace
      integer(HID_T) :: dset_id
      character(100) :: filename, filen
     !----------------------------------------------------------
      
      filename=trim(dir)//'/restvar.h5'
       

     !check restart files
      inquire(file=trim(filename), exist=dire1)

      filen=trim(dir)//'/restgq.txt'
      inquire(file=trim(filename), exist=dire2)
      
      if ((.not.dire1).or.(.not.dire2)) then
        write(*,*)
        write(*,*) ' ...a restart file is missing; '// &
                'SIMULATION ABORTED'
        stop
      endif

     
      open(93, file=filen, action="read", status='old')
      read(93,*) t
      read(93,*) ns
      read(93,*) ndof
      close(93)

    
     !check if the current NURBS refinement corresponds 
      if (ndof/=na) then
        write(*,*)
        write(*,*) ' ...current refinement does not corresponds to'// &
          'restarted values; SIMULATION ABORTED'
        stop
      endif
      
      
      nd=1
      dims=na

      CALL h5open_f(hdferr)

      CALL h5fopen_f(filename,H5F_ACC_RDONLY_F,file_id,hdferr)
      
      CALL h5dopen_f(file_id,'pot',dset_id,hdferr)
      CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE,pot,dims,hdferr)
      CALL h5dclose_f(dset_id,hdferr)

      CALL h5dopen_f(file_id,'wrec',dset_id,hdferr)
      CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE,wrec,dims,hdferr)
      CALL h5dclose_f(dset_id,hdferr)

      CALL h5dopen_f(file_id,'iionnm1',dset_id,hdferr)
      CALL h5dread_f(dset_id,H5T_NATIVE_DOUBLE,iionnm1,dims,hdferr)
      CALL h5dclose_f(dset_id,hdferr)

      CALL h5fclose_f(file_id,hdferr)


      ENDSUBROUTINE 


