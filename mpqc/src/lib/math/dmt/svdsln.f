	subroutine svdsln( A, lda, n, b, x, lsvec, rsvec,
     1			   gamma, rsdnrm, rcndtl, debug)
c
c------------------------------------------------------------------------------
c	svdsln
c
c	compute the solution to Ax=b by using the SVD of A
c
c	Parameters
c	  A      --> matrix A
c	  lda    --> leading dimension of A
c	  n      --> number of columns of A
c	  b      --> right hand side
c	  x      <-- on return contains the solution to Ax=b
c	  lsvec  <-- left singular vectors of A
c	  rsvec  <-- right singular vectors of A
c	  gamma  <-- value of the component in the deflated solution
c	  rsdnrm <-- norm of the residual (b-Ax)
c	  rcndtl--> tolerance on the reciprocal condition number
c	  debug  --> flag which controls debug printing statements
c
c
c       Notes
c         Maximum problem size is determined by WKMAX.  
c         To increase the size change the value of WKMAX
c
c 	Uses
c	  DSVDC      LINPACK routine, computes an SVD decompostion
c
c------------------------------------------------------------------------------
c
	integer lda, n, debug
c
        integer WKMAX
        parameter (WKMAX = 20)
	real*8  A(lda,n), Vwork(WKMAX,WKMAX)
	real*8  b(n), x(n), lsvec(n), rsvec(n)
	real*8	rcndtl, rsdnrm, gamma
c
	integer jobsvd, info, i, p
        real*8  svdwk1(WKMAX), svdwk2(WKMAX), sing(WKMAX)
	real*8  rconda, alpha
        real*8  ddot
c
c       Check sizes
c
        if (n .gt. WKMAX) then
            write(6,*) 'The matrix is too big. n = ', n
            write(6,*) 'The maximum size is      = ', WKMAX
            return
        endif
c
	jobsvd = 11
        ldv = WKMAX
	call dsvdc( A, lda, n, n, sing, svdwk1, A, lda, 
     1              Vwork, ldv, svdwk2, jobsvd, info )
c
c     Compute the condition number of the matrix
c
	p = 0
	do 50 i=1,n
	   rconda = sing(i)/sing(1)
	   if( rconda .lt. rcndtl ) then
		p = i-1
	        goto 100
	   endif
   50	continue
	p = n
  100	continue
c
	do 150 i=1,n
	   x(i)    = 0.0	
	   rsvec(i) = Vwork(i,n)
	   lsvec(i) = A(i,n)
  150	continue
c
C	write(6,*) 'Rank of matrix = ', p
	do 200 i=1,p
	   alpha = 1.0/sing(i) * ddot(n, A(1,i), 1, b, 1)
	   call daxpy( n, alpha, Vwork(1,i), 1, x, 1)
 200	continue
c
   	gamma = 1.0/sing(n) * ddot(n, A(1,n), 1, b, 1)
Cc
C	call amap(n, x, svdwk1)
C	do 300 i=1,n
C           svdwk2(i) = b(i) - svdwk1(i)
C  300	continue
Cc
C	rsdnrm = dnrm2(n,svdwk2,1)
Cc
	if( debug .gt. 4) then
	write(6,*) 'SVDSLN: Deflated solution using SVD'
	write(6,*) '   i	x		Rsd         sing. values'
	do 9000 i=1,n
	   write(6,9001) i, x(i), svdwk2(i), sing(i)
 9001	   format( 1x, I3, 3(1Pe13.6, 3x))
 9000	continue
	write(6,*) '|| Residual || = ', rsdnrm
	endif
c
	return
	end
