      subroutine ljfcn(nangles, angle, energy)
      integer nangles
      real*8  angle(nangles)
      real*8  energy
c------------------------------------------------------------------------------
c
c     ljfcn
c
c     purpose
c
c       Compute the energy associated with a 2D linear chain polymer
c       using a Lennard-Jones potential.
c       
c     parameters
c
c       nangles    ---> the number of angles
c       angle      ---> the angles describing the configuration 
c                       of the molecule
c
c     On exit
c
c       energy    <---  contains the energy associated with the molecule
c
c
c     Subroutine used
c
c       dist            computes the distance between two atoms
c       angtoc          converts from internal coordinates specified by
c                       angles to a cartesian coordinate system
c
c------------------------------------------------------------------------------

      integer MAXATOMS
      parameter (MAXATOMS = 100)

      integer NDIM
      parameter (NDIM     = 3)

      integer i, j, natoms
      integer nd
      real*8  coord(NDIM,MAXATOMS)
      real*8  dist
      real*8  r, rbar, rsub, r6, r12
      real*8  a, a6, a12
      real*8  y, f, width

      real*8  zero, one, two
      data    zero /0.0d0/, one /1.0d0/, two /2.0d0/

      common /fcncom/ nfeval
      integer nfeval
      save fcncom 
c
c     Set local variables
c
      natoms = nangles + 2
      a      = one
      nd     = NDIM
      rbar   = .06d0
      width  = .1d0
      energy = zero

c
c     convert angles to cartesian coordinates
c
      call angtoc(nangles, angle, coord, nd)

c-------------------------------------------------------------------------
c
c     Lennard-Jones Potential with a cutoff
c
c-------------------------------------------------------------------------

      do 200 i=2, natoms
         do 100 j=1, i-1

            r = dist(nd, coord(1,i), coord(1,j) )
            y = (r - rbar) / width

            if (y .lt. zero) then
               f = zero
            else if (y .gt. one) then
               f = one
            else
               f = y**3 * (10.0d0 + y*( 6.0d0*y - 15.0d0))
            endif
            
            if (r .lt. rbar) then
               rsub = rbar
            else if (rbar .le. r .and. r .le. (rbar+width)) then
               rsub = rbar + f*width
            else if (r .gt. (rbar+width)) then
               rsub = r
            endif
            
            r6   = rsub**6
            r12  = r6**2
c
c           Treat special case of a = 1
c
            if (a .eq. one) then
               energy = energy + (one/r12 - two/r6)
            else
               a6   = a**6
               a12  = a6**2
               energy = energy + (a12/r12 - two*a6/r6)
            endif
 100     continue
 200  continue

      nfeval = nfeval + 1
      return
      end

      real*8 function dist(n, x1, x2)
      integer n
      real*8  x1(n), x2(n)
c------------------------------------------------------------------------------
c
c     Compute the distance between two points in n-dimensional space
c
c------------------------------------------------------------------------------

      integer i

      dist = 0.0d0
      do 100 i=1,n
         dist = dist + (x1(i) - x2(i))**2
 100  continue

      dist = sqrt(dist)

      return
      end

      subroutine angtoc(nangles, angle, coord, ndim)
      integer nangles, ndim
      real*8  angle(nangles)
      real*8  coord(ndim,nangles+2)
c----------------------------------------------------------------------------
c
c
c     convert angles to cartesian coordinates
c 
c     Angles are assumed to be in radians
c
c----------------------------------------------------------------------------

      integer ix, iy, iz
      integer natoms
      real*8  v1(3), v2(3), angv1
      integer i, j, k
      real*8  zero, one, two
      data    zero /0.0d0/, one /1.0d0/, two /2.0d0/


      natoms = nangles + 2

      ix = 1
      iy = 2
      iz = 3

      coord(ix,1) =  zero
      coord(iy,1) =  zero
      coord(iz,1) =  zero

      coord(ix,2) =  zero
      coord(iy,2) = -one
      coord(iz,2) =  zero

      v1(ix) =  zero
      v1(iy) = -one
      v1(iz) =  zero

      v2(iz) =  zero

      do 100 i = 3, natoms

         v2(ix) = cos(angle(i-2))
         v2(iy) = sin(angle(i-2))

         angv1  = -atan2(-v1(iy), -v1(ix))
         v1(ix) =  cos(angv1)*v2(ix) + sin(angv1)*v2(iy)
         v1(iy) = -sin(angv1)*v2(ix) + cos(angv1)*v2(iy)

         coord(ix,i) =  coord(ix,i-1) + v1(ix)
         coord(iy,i) =  coord(iy,i-1) + v1(iy)
         coord(iz,i) =  coord(iz,i-1) + v1(iz)

 100     continue

      return
      end
      subroutine ljjac (n, x, fx, g)
      integer n
      real*8  x(n), g(n)
      real*8  fx
c------------------------------------------------------------------------------
c
c     ljjac
c     routine to evaluate function (f) and gradient (g) of 
c     the Lennard-Jones function at the point x using finite-differences
c
c------------------------------------------------------------------------------
c
      integer i, imeth
      real*8  eps, mcheps, hi, xtmp, typx
      real*8  fplus, fminus

      data    mcheps /1.1d-16/
      data    typx   /1.0d0/
c
c     Compute function at x
c
      call ljfcn(n,x,fx)
c
c     Compute gradient at x
c     For now hard-code forward-differences
c
      imeth  = 1

      if (imeth .eq. 1) then
c
c     Using first order forward finite-differences
c
         eps = sqrt(mcheps)
         do 10 i = 1,n
            hi    = eps*max(abs(x(i)),typx)
            xtmp  = x(i)
            x(i)  = xtmp + hi
            call ljfcn(n,x,fplus)
            g(i)  = (fplus - fx) / hi
            x(i)  = xtmp
 10      continue
      
      else
c
c     Using first order central finite-differences
c

         eps = mcheps**(1.0d0/3.0d0)
         do 20 i   = 1,n
            hi     = eps*max(abs(x(i)),typx)
            xtmp   = x(i)
            x(i)   = xtmp + hi
            call ljfcn(n,x,fplus)
            x(i)   = xtmp - hi
            call ljfcn(n,x,fminus)
            g(i)   = (fplus - fminus) / (2.0d0*hi)
            x(i)   = xtmp
 20      continue
      
      endif

      return
      end
      subroutine ljjac2 (n, x, fx, g)
      integer n
      real*8  x(n), g(n)
      real*8  fx
c------------------------------------------------------------------------------
c
c     ljjac2
c     routine to evaluate function (f) and gradient (g) of 
c     the Lennard-Jones function at the point x 
c     This version uses the analytic gradients supplied from adifor
c
c------------------------------------------------------------------------------
c
      integer i, j, imeth
      real*8  eps, mcheps, hi, xtmp, typx
      real*8  fplus, fminus

      real*8  zero, one, two
      data    zero /0.0d0/, one /1.0d0/, two /2.0d0/

      logical first
      data    first /.true./
      save    first

      integer ldg$x, ldg$energy
      real*8  energy, g$x(100,100)
      real*8  xsave(100), g$energy(100)
c
c     Compute function at x
c
c      call ljfcn(n,x,fx)
c
c     Compute gradient at x
c     Interface to g$energy$6 
c
      ldg$x = 100
      ldg$energy = 100

      do i=1,n
         xsave(i) = x(i)
      enddo

      do i=1,n
         do j=1,n
            g$x(i,j) = zero
            if (i .eq. j) g$x(i,j) = one
         enddo
      enddo
      call g$ljfcn$6(n, n, xsave, g$x, ldg$x,
     $               energy, g$energy, ldg$energy)

      fx = energy
      do i=1,n
         g(i) = g$energy(i)
      enddo

      return
      end
