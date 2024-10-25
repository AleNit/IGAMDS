
     !fill Compressed Sparse Row (CSR) format arrays of the
     !input matrix; the input matrix is given as dense matrix

      SUBROUTINE CRSarr(n,nnz,A,Aa,IA,JA)

      implicit none
     !----------------------------------------------------- 
      integer, intent(in) :: n, nnz
      real, intent(in) :: A(n,n)
      integer, intent(out) :: IA(n+1), JA(nnz)
      real, intent(out) :: Aa(nnz)
     !----------------------------------------------------- 
      integer :: i, j, k, kk, nnzj
     !----------------------------------------------------- 

      k=1
      kk=1
      IA(1)=1
      nnzj=0
      do i=1, n
        do j=1, n
          if (A(i,j)/=0.0) then
            Aa(k)=A(i,j)
            JA(k)=j
            nnzj=nnzj+1
            k=k+1
          endif
        enddo
        IA(kk+1)=IA(kk)+nnzj
        kk=kk+1
        nnzj=0
      enddo


      ENDSUBROUTINE

     !=====================================================

     !generate CSR arrays from COO arrays; exploit bubble sort

      SUBROUTINE COO2CSR(nd,ndof,nnz,VAcoo,IAcoo,JAcoo,VA,IA,JA)

      implicit none
     !-----------------------------------------------------
      integer, intent(in) :: ndof,nnz,nd
      integer, intent(inout) :: IAcoo(nd),JAcoo(nd)
      real, intent(inout) :: VAcoo(nd)
      integer, intent(out) :: IA(ndof+1),JA(nnz)
      real, intent(out) :: VA(nnz)
     !-----------------------------------------------------
      integer :: i,c,cc,k,IAtmp(nnz)
      real :: t1,t2,t3
     !-----------------------------------------------------

      CALL cpu_time(t1)

     !sort COO arrays by row inclusing zeros 
      CALL quicksort(IAcoo,1,nd,VAcoo,JAcoo)

      VA    = VAcoo(nd-nnz+1:nd)
      JA    = JAcoo(nd-nnz+1:nd)
      IAtmp = IAcoo(nd-nnz+1:nd)

      CALL cpu_time(t2)

     !for each row sort the COO arrays by column; IAtmp will be destroyed
      c=1
      cc=0
      k=2
      IA(1)=1
      do i=2, nnz        
        cc=cc+1
        if (IAtmp(i)/=IAtmp(i-1)) then
          CALL Bubble_Sort(VA(c:i-1),JA(c:i-1),IAtmp(c:i-1),i-c)
          IA(k)=IA(k-1)+cc
          c=i
          k=k+1
          cc=0
        endif
      enddo
      IA(ndof+1)=nnz+1

      CALL cpu_time(t3)

!      write(*,*) 'outer sort time',t2-t1
!      write(*,*) 'inner sort time',t3-t2
      

      ENDSUBROUTINE
      
     !===================================================== 

      SUBROUTINE amux ( n, x, y, a, ja, ia )
      
      !*****************************************************************************80
      !
      !! AMUX multiplies a CSR matrix A times a vector.
      !
      !  Discussion:
      !
      !    This routine multiplies a matrix by a vector using the dot product form.
      !    Matrix A is stored in compressed sparse row storage.
      !
      !  Modified:
      !
      !    07 January 2004
      !
      !  Author:
      !
      !    Youcef Saad
      !
      !  Parameters:
      !
      !    Input, integer ( kind = 4 ) N, the row dimension of the matrix.
      !
      !    Input, real X(*), and array of length equal to the column dimension 
      !    of A.
      !
      !    Input, real A(*), integer ( kind = 4 ) JA(*), IA(NROW+1), the matrix in CSR
      !    Compressed Sparse Row format.
      !
      !    Output, real Y(N), the product A * X.
      !
        implicit none
      
        integer :: n
      
        real ::  a(*)
        integer :: i
        integer :: ia(*)
        integer :: ja(*)
        integer :: k
        real ::  t
        real ::  x(*)
        real ::  y(n)
      
        do i = 1, n
      
          t = 0.0D+00
          do k = ia(i), ia(i+1)-1
            t = t + a(k) * x(ja(k))
          end do
      
          y(i) = t
      
        end do
      
      ENDSUBROUTINE

     !===================================================== 

      SUBROUTINE indint(ind,c,nd,IA,JA,ii,jj)

      implicit none
     !-----------------------------------------------------
      integer, intent(in) :: nd,ii,jj,IA(nd),JA(nd)
      integer, intent(inout) :: c
      integer, intent(out) :: ind
     !-----------------------------------------------------
      integer :: i
     !-----------------------------------------------------
      
      do i=c,1,-1
        if ( IA(i)==ii .and. JA(i)==jj ) then
          ind=i
          return
        endif
      enddo

      c=c+1
      ind=c

      ENDSUBROUTINE

     !=====================================================

      SUBROUTINE Bubble_Sort(a,indi,indj,k)

      implicit none
     !----------------------------------------------------- 
      integer, intent(in) :: k
      REAL, INTENT(in out), DIMENSION(k) :: a
      integer, INTENT(in out), DIMENSION(k) :: indi,indj
     !----------------------------------------------------- 
      REAL :: temp  
      INTEGER :: i, j, tempi, tempj
      LOGICAL :: swapped
     !----------------------------------------------------- 
 

      DO j = SIZE(indi)-1, 1, -1
        swapped = .FALSE.
        DO i = 1, j
          IF (indi(i) > indi(i+1)) THEN
    
            tempi = indi(i)
            indi(i) = indi(i+1)
            indi(i+1) = tempi
            
            tempj = indj(i)
            indj(i) = indj(i+1)
            indj(i+1) = tempj
    
            temp = a(i)
            a(i) = a(i+1)
            a(i+1) = temp
    
            swapped = .TRUE.
    
          END IF
        END DO
        IF (.NOT. swapped) EXIT
      END DO

   
      ENDSUBROUTINE

!     !===================================================== 
!      
!     ! Sorts an array arr into ascending numerical order using the Quicksort algorithm. arr is
!     ! replaced on output by its sorted rearrangement.
!     ! Parameters: NN is the size of subarrays sorted by straight insertion and NSTACK is the
!     ! required auxiliary storage.
!
!      SUBROUTINE quicksort(n,arr)
!
!      IMPLICIT NONE
!     !-----------------------------------------------------
!      INTEGER, INTENT(IN) :: n
!      REAL, DIMENSION(n), INTENT(INOUT) :: arr
!      INTEGER, PARAMETER :: NN=15, NSTACK=50
!      REAL :: a
!      INTEGER :: k,i,j,jstack,l,r
!      INTEGER, DIMENSION(NSTACK) :: istack
!     !-----------------------------------------------------
!
!      jstack=0
!      l=1
!      r=n
!      
!      do
!
!        if (r-l < NN) then !Insertion sort when subarray small enough.
!          do j=l+1,r
!            a=arr(j)
!            do i=j-1,l,-1
!              if (arr(i) <= a) exit
!              arr(i+1)=arr(i)
!            enddo
!            arr(i+1)=a
!          enddo
!          if (jstack == 0) RETURN
!          r=istack(jstack) !Pop stack and begin a new round of partitionl=
!          istack(jstack-1) .
!          jstack=jstack-2
!
!        else !Choose median of left, center, and right elements
!             !as partitioning element a. Also rearrange so
!             !that a(l) ≤ a(l+1) ≤ a(r).
!          k=(l+r)/2
!          call swap(arr(k),arr(l+1))
!          call swap(arr(l),arr(r),arr(l)>arr(r))
!          call swap(arr(l+1),arr(r),arr(l+1)>arr(r))
!          call swap(arr(l),arr(l+1),arr(l)>arr(l+1))
!          i=l+1 !Initialize pointers for partitioning.
!          j=r
!          a=arr(l+1) !Partitioning element.
!      
!          do !Here is the meat.
!            do !Scan up to find element >= a.
!              i=i+1
!              if (arr(i) >= a) exit
!            enddo
!            do !Scan down to find element <= a.
!              j=j-1
!              if (arr(j) <= a) exit
!            enddo
!            if (j < i) exit !Pointers crossed. Exit with partitioning complete.
!            call swap(arr(i),arr(j)) !Exchange elements.
!          enddo
!          arr(l+1)=arr(j) !Insert partitioning element.
!          arr(j)=a
!          jstack=jstack+2
!      
!      !Push pointers to larger subarray on stack; process smaller subarray immediately.
!          if (jstack > NSTACK) then
!            write(*,*) '... QuickSort algorithm failed'
!            stop
!          endif  
!          if (r-i+1 >= j-l) then
!            istack(jstack)=r
!            istack(jstack-1)=i
!            r=j-1
!          else
!            istack(jstack)=j-1
!            istack(jstack-1)=l
!            l=i
!          endif
!        endif
!      
!      enddo
!
!
!      ENDSUBROUTINE
!
!     !=====================================================
!
!      SUBROUTINE swap_i(a,b)
!
!      implicit none
!      integer, intent(inout) :: a,b
!      integer :: dum
!
!      dum=a
!      a=b
!      b=dum
!
!      ENDSUBROUTINE
!     
!     !=====================================================

      ! QuickSort algorithm on integer array a, with two slave arrays,
      ! one real and one integer
      ! Author: t-nissie
      ! License: GPLv3
      ! Gist: https://gist.github.com/t-nissie/479f0f16966925fa29ea

       RECURSIVE SUBROUTINE quicksort(a,first,last,slaver,slavei)
       
       implicit none
!       real*8  a(*), x, t
       integer :: a(*), x, ti, slavei(*)
       real :: slaver(*), tr
       integer :: first, last
       integer :: i, j
     
       x = a( (first+last) / 2 )
       i = first
       j = last
       do
          do while (a(i) < x)
             i=i+1
          end do
          do while (x < a(j))
             j=j-1
          end do
          if (i >= j) exit          
          ti = a(i);       a(i) = a(j);            a(j) = ti   !value swap
          ti = slavei(i);  slavei(i) = slavei(j);  slavei(j) = ti
          tr = slaver(i);  slaver(i) = slaver(j);  slaver(j) = tr
          i=i+1
          j=j-1
       end do
       if (first < i-1) call quicksort(a,first,i-1,slaver,slavei)
       if (j+1 < last)  call quicksort(a,j+1,last,slaver,slavei)
     
       ENDSUBROUTINE

     !=====================================================

