************************************
* Friction and mobility            *
* for a hard-sphere configuration  *
* Module PBC                       *
*                                  *
* K. Hinsen                        *
* Last revision: December 7, 1993  *
************************************


*********************************************************************
* Determine coefficients for the calculation of the Hasimoto tensor *
*********************************************************************

* differentiation

      subroutine difx(c,n,twice)
      implicit none
      integer n
      double precision c(0:n,0:n,0:n,4)
      logical twice

      integer i,j,k,l

      do 100 l = 1,4
         do 110 j = 0,n
            do 120 k = 0,n
               do 130 i = 0,n-1
                  c(i,j,k,l) = 2*(i+1)*c(i+1,j,k,l)
                  if (twice) c(i,j,k,l) = (2*i+1)*c(i,j,k,l)
130            continue
               c(n,j,k,l) = 0.d0
120         continue
110      continue
100   continue

      end

      subroutine dify(c,n,twice)
      implicit none
      integer n
      double precision c(0:n,0:n,0:n,4)
      logical twice

      integer i,j,k,l

      do 100 l = 1,4
         do 110 i = 0,n
            do 120 k = 0,n
               do 130 j = 0,n-1
                  c(i,j,k,l) = 2*(j+1)*c(i,j+1,k,l)
                  if (twice) c(i,j,k,l) = (2*j+1)*c(i,j,k,l)
130            continue
               c(i,n,k,l) = 0.d0
120         continue
110      continue
100   continue

      end

      subroutine difz(c,n,twice)
      implicit none
      integer n
      double precision c(0:n,0:n,0:n,4)
      logical twice

      integer i,j,k,l

      do 100 l = 1,4
         do 110 i = 0,n
            do 120 j = 0,n
               do 130 k = 0,n-1
                  c(i,j,k,l) = 2*(k+1)*c(i,j,k+1,l)
                  if (twice) c(i,j,k,l) = (2*k+1)*c(i,j,k,l)
130            continue
               c(i,j,n,l) = 0.d0
120         continue
110      continue
100   continue

      end


* copy coefficient arrays

      subroutine cocopy(c1,c2,n)
      implicit none
      integer n
      double precision c1(0:n,0:n,0:n,4),c2(0:n,0:n,0:n,4)

      integer i,j,k,l

      do 100 l = 1,4
         do 110 i = 0,n
            do 120 j = 0,n
               do 130 k = 0,n
                  c1(i,j,k,l) = c2(i,j,k,l)
130            continue
120         continue
110      continue
100   continue

      end


* compress coefficients

      subroutine compr(cc,c,n0,n1,n2,n3)
      implicit none
      integer n0,n1,n2,n3
      double precision cc(n0,4,2),c(0:n1,0:n1,0:n1,4)

      integer ind1,ind2
      integer i,j,k,l

      do 100 l = 1,4
         ind1 = 1
         ind2 = 1
         do 110 i = 0,n2
            do 120 j = 0,n2-i
               do 130 k = 0,n2-i-j
                  cc(ind1,l,1) = c(i,j,k,l)
                  ind1 = ind1 + 1
                  if (i+j+k .le. n3) then
                     cc(ind2,l,2) = c(i,j,k,l)
                     ind2 = ind2 + 1
                  endif
130            continue
120         continue
110      continue
100   continue

      end


* load coefficients for s1 and s2

      subroutine cload
      implicit none

      double precision s1(0:9,0:9,0:9,4),s2(0:9,0:9,0:9,4)
      common /hcoeff/ s1,s2

      integer i,j,k,l,max1,max2

      open (unit=10,file=
     .FILENAME2
     .     ,form='FORMATTED',status='OLD')
      read(10,*) max1
      if (max1 .ne. 8) stop 'wrong data file!'
      do 100 l = 1,4
         do 110 i = 0,8
            do 120 j = 0,8
               do 130 k = 0,8
                  if (i+j+k .le. 8) read (10,*) s1(i,j,k,l)
130            continue
120         continue
110      continue
100   continue

      read(10,*) max2
      if (max2 .ne. 9) stop 'wrong data file!'
      do 200 l = 1,4
         do 210 i = 0,9
            do 220 j = 0,9
               do 230 k = 0,9
                  if (i+j+k .le. 9) read (10,*) s2(i,j,k,l)
230            continue
220         continue
210      continue
200   continue
      close(10)

      end


* generate coefficients for dipole interactions (periodic bc)

      subroutine pbccoef
      implicit none

      double precision s1c(0:9,0:9,0:9,4),s2c(0:9,0:9,0:9,4)
      common /hcoeff/ s1c,s2c

      double precision s1(84,4,2)
      double precision s1x(84,4,2),s1z(84,4,2)
      double precision s1xx(120,4,2),s1zz(120,4,2),s1xy(120,4,2),
     .     s1xz(120,4,2)
      double precision s2xx(84,4,2),s2zz(84,4,2),s2xy(84,4,2),
     .     s2xz(84,4,2)
      double precision s2xxx(84,4,2),s2zzz(84,4,2),s2xxy(84,4,2)
      double precision s2xxz(84,4,2),s2zzx(84,4,2),s2xyz(84,4,2)
      double precision s2xxxx(120,4,2),s2zzzz(120,4,2),s2xxyy(120,4,2)
      double precision s2xxzz(120,4,2),s2xxxy(120,4,2),s2xxxz(120,4,2)
      double precision s2zzzx(120,4,2),s2xxyz(120,4,2),s2zzxy(120,4,2)
      common /pbcc/ s1,s1x,s1z,s1xx,s1zz,s1xy,s1xz,s2xx,s2zz,s2xy,
     .              s2xz,s2xxx,s2zzz,s2xxy,s2xxz,s2zzx,s2xyz,
     .              s2xxxx,s2zzzz,s2xxyy,s2xxzz,s2xxxy,s2xxxz,
     .              s2zzzx,s2xxyz,s2zzxy

      integer i,j
      double precision c(0:9,0:9,0:9,4)

      call compr(s1,s1c,84,9,6,4)

      call cocopy(c,s1c,9)
      call difx(c,9,.false.)
      call compr(s1x,c,84,9,6,4)
      call cocopy(c,s1c,9)
      call difz(c,9,.false.)
      call compr(s1z,c,84,9,6,4)

      call cocopy(c,s1c,9)
      call difx(c,9,.true.)
      call compr(s1xx,c,120,9,7,5)
      call cocopy(c,s1c,9)
      call difz(c,9,.true.)
      call compr(s1zz,c,120,9,7,5)
      call cocopy(c,s1c,9)
      call difx(c,9,.false.)
      call dify(c,9,.false.)
      call compr(s1xy,c,120,9,7,5)
      call cocopy(c,s1c,9)
      call difx(c,9,.false.)
      call difz(c,9,.false.)
      call compr(s1xz,c,120,9,7,5)

      call cocopy(c,s2c,9)
      call difx(c,9,.true.)
      call compr(s2xx,c,84,9,6,4)
      call cocopy(c,s2c,9)
      call difz(c,9,.true.)
      call compr(s2zz,c,84,9,6,4)
      call cocopy(c,s2c,9)
      call difx(c,9,.false.)
      call dify(c,9,.false.)
      call compr(s2xy,c,84,9,6,4)
      call cocopy(c,s2c,9)
      call difx(c,9,.false.)
      call difz(c,9,.false.)
      call compr(s2xz,c,84,9,6,4)

      call cocopy(c,s2c,9)
      call difx(c,9,.true.)
      call difx(c,9,.false.)
      call compr(s2xxx,c,84,9,6,4)
      call cocopy(c,s2c,9)
      call difz(c,9,.true.)
      call difz(c,9,.false.)
      call compr(s2zzz,c,84,9,6,4)
      call cocopy(c,s2c,9)
      call difx(c,9,.true.)
      call dify(c,9,.false.)
      call compr(s2xxy,c,84,9,6,4)
      call cocopy(c,s2c,9)
      call difx(c,9,.true.)
      call difz(c,9,.false.)
      call compr(s2xxz,c,84,9,6,4)
      call cocopy(c,s2c,9)
      call difz(c,9,.true.)
      call difx(c,9,.false.)
      call compr(s2zzx,c,84,9,6,4)
      call cocopy(c,s2c,9)
      call difx(c,9,.false.)
      call dify(c,9,.false.)
      call difz(c,9,.false.)
      call compr(s2xyz,c,84,9,6,4)

      call cocopy(c,s2c,9)
      call difx(c,9,.true.)
      call difx(c,9,.true.)
      call compr(s2xxxx,c,120,9,7,5)
      call cocopy(c,s2c,9)
      call difz(c,9,.true.)
      call difz(c,9,.true.)
      call compr(s2zzzz,c,120,9,7,5)
      call cocopy(c,s2c,9)
      call difx(c,9,.true.)
      call dify(c,9,.true.)
      call compr(s2xxyy,c,120,9,7,5)
      call cocopy(c,s2c,9)
      call difx(c,9,.true.)
      call difz(c,9,.true.)
      call compr(s2xxzz,c,120,9,7,5)
      call cocopy(c,s2c,9)
      call difx(c,9,.true.)
      call difx(c,9,.false.)
      call dify(c,9,.false.)
      call compr(s2xxxy,c,120,9,7,5)
      call cocopy(c,s2c,9)
      call difx(c,9,.true.)
      call difx(c,9,.false.)
      call difz(c,9,.false.)
      call compr(s2xxxz,c,120,9,7,5)
      call cocopy(c,s2c,9)
      call difz(c,9,.true.)
      call difz(c,9,.false.)
      call difx(c,9,.false.)
      call compr(s2zzzx,c,120,9,7,5)
      call cocopy(c,s2c,9)
      call difx(c,9,.true.)
      call dify(c,9,.false.)
      call difz(c,9,.false.)
      call compr(s2xxyz,c,120,9,7,5)
      call cocopy(c,s2c,9)
      call difz(c,9,.true.)
      call difx(c,9,.false.)
      call dify(c,9,.false.)
      call compr(s2zzxy,c,120,9,7,5)

      do 100 i = 1,2
         do 110 j = 1,84
            s1x(j,3,i) = -s1x(j,3,i)
            s1x(j,4,i) = -s1x(j,4,i)
            s1z(j,2,i) = -s1z(j,2,i)
            s1z(j,4,i) = -s1z(j,4,i)
            s2xz(j,2,i) = -s2xz(j,2,i)
            s2xz(j,3,i) = -s2xz(j,3,i)
            s2xxx(j,3,i) = -s2xxx(j,3,i)
            s2xxx(j,4,i) = -s2xxx(j,4,i)
            s2zzz(j,2,i) = -s2zzz(j,2,i)
            s2zzz(j,4,i) = -s2zzz(j,4,i)
            s2xxy(j,3,i) = -s2xxy(j,3,i)
            s2xxy(j,4,i) = -s2xxy(j,4,i)
            s2xxz(j,2,i) = -s2xxz(j,2,i)
            s2xxz(j,4,i) = -s2xxz(j,4,i)
            s2zzx(j,3,i) = -s2zzx(j,3,i)
            s2zzx(j,4,i) = -s2zzx(j,4,i)
            s2xyz(j,2,i) = -s2xyz(j,2,i)
            s2xyz(j,4,i) = -s2xyz(j,4,i)
110      continue
100   continue

      do 120 i = 1,2
         do 130 j = 1,120
            s1xz(j,2,i) = -s1xz(j,2,i)
            s1xz(j,3,i) = -s1xz(j,3,i)
            s2xxxz(j,2,i) = -s2xxxz(j,2,i)
            s2xxxz(j,3,i) = -s2xxxz(j,3,i)
            s2zzzx(j,2,i) = -s2zzzx(j,2,i)
            s2zzzx(j,3,i) = -s2zzzx(j,3,i)
            s2xxyz(j,2,i) = -s2xxyz(j,2,i)
            s2xxyz(j,3,i) = -s2xxyz(j,3,i)
130      continue
120   continue

      end


************************************
* Calculate the interaction matrix *
************************************

* nearest neighbour images

      blockdata images

      integer numimg(4),nnimg(0:2,8,4)
      common /nn/ numimg,nnimg

      data numimg /1,2,4,8/
      data nnimg /0,0,0,  21*0,
     .            0,0,0, 0,0,1, 18*0,
     .            0,0,0, 0,1,0, 1,0,0, 1,1,0, 12*0,
     .            0,0,0, 0,0,1, 0,1,0, 1,0,0,
     .                1,1,0, 1,0,1, 0,1,1, 1,1,1/

      end

* Calculate powers for dipole interactions (periodic bc)

      subroutine dpower(p1,p2,n0,index,xsq,ysq,zsq,n1,n2)
      implicit none
      double precision p1(*),p2(*)
      integer n0,n1,n2,index
      double precision xsq,ysq,zsq

      double precision xp(0:9),yp(0:9),zp(0:9)
      integer i,j,k,n

      xp(0) = 1.d0
      yp(0) = 1.d0
      zp(0) = 1.d0
      do 100 i = 1,n1
         xp(i) = xsq * xp(i-1)
         yp(i) = ysq * yp(i-1)
         zp(i) = zsq * zp(i-1)
100   continue

      if (xsq+ysq+zsq .gt. 0.0625d0) then
         n = n1
         index = 1
      else
         n = n2
         index = 2
      endif
      n0 = 1
      do 200 i = 0,n
         do 210 j = 0,n-i
            do 220 k = 0,n-i-j
               p1(n0) = xp(i)*yp(j)*zp(k)
               p2(n0) = xp(j)*yp(i)*zp(k)
               n0 = n0+1
220         continue
210      continue
200   continue
      n0 = n0-1

      end

* evaluate power series

      double precision function sum(c,p,n)
      implicit none
      integer n
      double precision c(*),p(*)

      integer i

      sum = 0.d0
      do 100 i = 1,n
         sum = sum + c(i)*p(i)
100   continue

      end

* swap two numbers

      subroutine swap(a,b)
      implicit none
      double precision a,b

      double precision temp

      temp = a
      a = b
      b = temp

      end

* calculate interaction matrix with pbc (up to 1/r^3)

      subroutine gpbc(xc,yc,zc)
      implicit none
      double precision xc,yc,zc
 
      double precision s1(84,4,2)
      double precision s1x(84,4,2),s1z(84,4,2)
      double precision s1xx(120,4,2),s1zz(120,4,2),s1xy(120,4,2),
     .     s1xz(120,4,2)
      double precision s2xx(84,4,2),s2zz(84,4,2),s2xy(84,4,2),
     .     s2xz(84,4,2)
      double precision s2xxx(84,4,2),s2zzz(84,4,2),s2xxy(84,4,2)
      double precision s2xxz(84,4,2),s2zzx(84,4,2),s2xyz(84,4,2)
      double precision s2xxxx(120,4,2),s2zzzz(120,4,2),s2xxyy(120,4,2)
      double precision s2xxzz(120,4,2),s2xxxy(120,4,2),s2xxxz(120,4,2)
      double precision s2zzzx(120,4,2),s2xxyz(120,4,2),s2zzxy(120,4,2)
      common /pbcc/ s1,s1x,s1z,s1xx,s1zz,s1xy,s1xz,s2xx,s2zz,s2xy,
     .     s2xz,s2xxx,s2zzz,s2xxy,s2xxz,s2zzx,s2xyz,
     .     s2xxxx,s2zzzz,s2xxyy,s2xxzz,s2xxxy,s2xxxz,
     .     s2zzzx,s2xxyz,s2zzxy
 
      integer numimg(4),nnimg(0:2,8,4)
      common /nn/ numimg,nnimg
 
      double precision box
      common /pbc/ box
 
      double precision g0s0s(3,3),g0s1s(3,5),g0s2s(3,7),g1s1s(5,5)
      double precision g0s0t(3,3),g0s1t(3,5),g0s0p(3,3),s1v
      common /gpbca/ g0s0s,g0s1s,g0s2s,g1s1s,g0s0t,g0s1t,g0s0p,s1v
 
      external sum
      double precision sum
 
      double precision d0s1,d1s1(3),d2s1(6)
      double precision d2s2(6),d3s2(10),d4s2(15)
      double precision p1(120),p2(120)
      double precision x,y,z,xx,yy,zz,xi,yi,zi,xs,ys,zs,xi2,yi2,zi2
      double precision r,r2,r3,box1,box2,box3,f1,f2,f21,f22,f23
      integer i,index,type,terms,ti
 
      x = dabs(xc/box)
      y = dabs(yc/box)
      z = dabs(zc/box)
      type = 0
      if (x.ge.0.25d0 .neqv.  y.ge.0.25d0) then
         if (x.ge.0.25d0 .neqv.  z.ge.0.25d0) then
            call swap(x,z)
            type = 1
         else
            call swap(y,z)
            type = 2
         endif
      endif
      index = 1
      if (x .ge. 0.25d0) index = index+1
      if (y .ge. 0.25d0) index = index+1
      if (z .ge. 0.25d0) index = index+1
      xx = min(x,0.5d0-x)
      yy = min(y,0.5d0-y)
      zz = min(z,0.5d0-z)
 
      call dpower(p1,p2,terms,ti,xx*xx,yy*yy,zz*zz,6,4)
 
      d0s1 = sum(s1(1,index,ti),p1,terms)
 
      d1s1(1) = xx*sum(s1x(1,index,ti),p1,terms)
      d1s1(2) = yy*sum(s1x(1,index,ti),p2,terms)
      d1s1(3) = zz*sum(s1z(1,index,ti),p1,terms)
 
      d2s2(1) = sum(s2xx(1,index,ti),p1,terms)
      d2s2(2) = xx*yy*sum(s2xy(1,index,ti),p1,terms)
      d2s2(3) = xx*zz*sum(s2xz(1,index,ti),p1,terms)
      d2s2(4) = sum(s2xx(1,index,ti),p2,terms)
      d2s2(5) = yy*zz*sum(s2xz(1,index,ti),p2,terms)
      d2s2(6) = sum(s2zz(1,index,ti),p1,terms)
 
      d3s2(1) = xx*sum(s2xxx(1,index,ti),p1,terms)
      d3s2(2) = yy*sum(s2xxy(1,index,ti),p1,terms)
      d3s2(3) = zz*sum(s2xxz(1,index,ti),p1,terms)
      d3s2(4) = xx*sum(s2xxy(1,index,ti),p2,terms)
      d3s2(5) = xx*yy*zz*sum(s2xyz(1,index,ti),p1,terms)
      d3s2(6) = xx*sum(s2zzx(1,index,ti),p1,terms)
      d3s2(7) = yy*sum(s2xxx(1,index,ti),p2,terms)
      d3s2(8) = zz*sum(s2xxz(1,index,ti),p2,terms)
      d3s2(9) = yy*sum(s2zzx(1,index,ti),p2,terms)
      d3s2(10) = zz*sum(s2zzz(1,index,ti),p2,terms)
 
      call dpower(p1,p2,terms,ti,xx*xx,yy*yy,zz*zz,7,5)
 
      d2s1(1) = sum(s1xx(1,index,ti),p1,terms)-4.18879020478639098458d0
      d2s1(2) = xx*yy*sum(s1xy(1,index,ti),p1,terms)
      d2s1(3) = xx*zz*sum(s1xz(1,index,ti),p1,terms)
      d2s1(4) = sum(s1xx(1,index,ti),p2,terms)-4.18879020478639098458d0
      d2s1(5) = yy*zz*sum(s1xz(1,index,ti),p2,terms)
      d2s1(6) = sum(s1zz(1,index,ti),p1,terms)-4.18879020478639098458d0
 
      d4s2(1) = sum(s2xxxx(1,index,ti),p1,terms)
     .     -2.51327412287183459075d0
      d4s2(2) = xx*yy*sum(s2xxxy(1,index,ti),p1,terms)
      d4s2(3) = xx*zz*sum(s2xxxz(1,index,ti),p1,terms)
      d4s2(4) = sum(s2xxyy(1,index,ti),p1,terms)
     .     -0.83775804095727819691d0
      d4s2(5) = yy*zz*sum(s2xxyz(1,index,ti),p1,terms)
      d4s2(6) = sum(s2xxzz(1,index,ti),p1,terms)
     .     -0.83775804095727819691d0
      d4s2(7) = xx*yy*sum(s2xxxy(1,index,ti),p2,terms)
      d4s2(8) = xx*zz*sum(s2xxyz(1,index,ti),p2,terms)
      d4s2(9) = xx*yy*sum(s2zzxy(1,index,ti),p1,terms)
      d4s2(10) = xx*zz*sum(s2zzzx(1,index,ti),p1,terms)
      d4s2(11) = sum(s2xxxx(1,index,ti),p2,terms)
     .     -2.51327412287183459075d0
      d4s2(12) = yy*zz*sum(s2xxxz(1,index,ti),p2,terms)
      d4s2(13) = sum(s2xxzz(1,index,ti),p2,terms)
     .     -0.83775804095727819691d0
      d4s2(14) = yy*zz*sum(s2zzzx(1,index,ti),p2,terms)
      d4s2(15) = sum(s2zzzz(1,index,ti),p1,terms)
     .     -2.51327412287183459075d0

      if (x*x+y*y+z*z .gt. 0.d0) then
         do 100 i = 1,numimg(index)
            xi = x-nnimg(0,i,index)
            yi = y-nnimg(1,i,index)
            zi = z-nnimg(2,i,index)
            r2 = 1.d0/(xi*xi+yi*yi+zi*zi)
            r = r2**0.5d0
            r3 = r*r*r
            xi = xi*r
            yi = yi*r
            zi = zi*r
            xi2 = xi*xi
            yi2 = yi*yi
            zi2 = zi*zi
 
            d0s1 = d0s1 + r
 
            d1s1(1) = d1s1(1) - r2*xi
            d1s1(2) = d1s1(2) - r2*yi
            d1s1(3) = d1s1(3) - r2*zi
 
            d2s1(1) = d2s1(1) + r3*(3.d0*xi2-1.d0)
            d2s1(2) = d2s1(2) + 3.d0*r3*xi*yi
            d2s1(3) = d2s1(3) + 3.d0*r3*xi*zi
            d2s1(4) = d2s1(4) + r3*(3.d0*yi2-1.d0)
            d2s1(5) = d2s1(5) + 3.d0*r3*yi*zi
            d2s1(6) = d2s1(6) + r3*(3.d0*zi2-1.d0)
 
            d2s2(1) = d2s2(1) + 0.5d0*r*(1.d0-xi2)
            d2s2(2) = d2s2(2) - 0.5d0*r*xi*yi
            d2s2(3) = d2s2(3) - 0.5d0*r*xi*zi
            d2s2(4) = d2s2(4) + 0.5d0*r*(1.d0-yi2)
            d2s2(5) = d2s2(5) - 0.5d0*r*yi*zi
            d2s2(6) = d2s2(6) + 0.5d0*r*(1.d0-zi2)
 
            d3s2(1) = d3s2(1) - 1.5d0*r2*xi*(1.d0-xi2)
            d3s2(2) = d3s2(2) - 0.5d0*r2*yi*(1.d0-3.d0*xi2)
            d3s2(3) = d3s2(3) - 0.5d0*r2*zi*(1.d0-3.d0*xi2)
            d3s2(4) = d3s2(4) - 0.5d0*r2*xi*(1.d0-3.d0*yi2)
            d3s2(5) = d3s2(5) + 1.5d0*r2*xi*yi*zi
            d3s2(6) = d3s2(6) - 0.5d0*r2*xi*(1.d0-3.d0*zi2)
            d3s2(7) = d3s2(7) - 1.5d0*r2*yi*(1.d0-yi2)
            d3s2(8) = d3s2(8) - 0.5d0*r2*zi*(1.d0-3.d0*yi2)
            d3s2(9) = d3s2(9) - 0.5d0*r2*yi*(1.d0-3.d0*zi2)
            d3s2(10) = d3s2(10) - 1.5d0*r2*zi*(1.d0-zi2)
 
            d4s2(1) = d4s2(1) + 1.5d0*r3*(xi2*(6.d0-5.d0*xi2)-1.d0)
            d4s2(2) = d4s2(2) + 1.5d0*r3*xi*yi*(3.d0-5.d0*xi2)
            d4s2(3) = d4s2(3) + 1.5d0*r3*xi*zi*(3.d0-5.d0*xi2)
            d4s2(4) = d4s2(4) + 0.5d0*r3*(-15.d0*xi2*yi2+3.d0
     .           *(xi2+yi2)-1.d0)
            d4s2(5) = d4s2(5) + 1.5d0*r3*yi*zi*(1.d0-5.d0*xi2)
            d4s2(6) = d4s2(6) + 0.5d0*r3*(-15.d0*xi2*zi2+3.d0
     .           *(xi2+zi2)-1.d0)
            d4s2(7) = d4s2(7) + 1.5d0*r3*xi*yi*(3.d0-5.d0*yi2)
            d4s2(8) = d4s2(8) + 1.5d0*r3*xi*zi*(1.d0-5.d0*yi2)
            d4s2(9) = d4s2(9) + 1.5d0*r3*xi*yi*(1.d0-5.d0*zi2)
            d4s2(10) = d4s2(10) + 1.5d0*r3*xi*zi*(3.d0-5.d0*zi2)
            d4s2(11) = d4s2(11) + 1.5d0*r3*(yi2*(6.d0-5.d0*yi2)-1.d0)
            d4s2(12) = d4s2(12) + 1.5d0*r3*yi*zi*(3.d0-5.d0*yi2)
            d4s2(13) = d4s2(13) + 0.5d0*r3*(-15.d0*yi2*zi2+3.d0
     .           *(yi2+zi2)-1.d0)
            d4s2(14) = d4s2(14) + 1.5d0*r3*yi*zi*(3.d0-5.d0*zi2)
            d4s2(15) = d4s2(15) + 1.5d0*r3*(zi2*(6.d0-5.d0*zi2)-1.d0)
100      continue
      endif
 
      if (type .eq. 1) then
         call swap(d1s1(1),d1s1(3))
         call swap(d2s1(1),d2s1(6))
         call swap(d2s1(2),d2s1(5))
         call swap(d2s2(1),d2s2(6))
         call swap(d2s2(2),d2s2(5))
         call swap(d3s2(1),d3s2(10))
         call swap(d3s2(2),d3s2(9))
         call swap(d3s2(3),d3s2(6))
         call swap(d3s2(4),d3s2(8))
         call swap(d4s2(1),d4s2(15))
         call swap(d4s2(2),d4s2(14))
         call swap(d4s2(3),d4s2(10))
         call swap(d4s2(4),d4s2(13))
         call swap(d4s2(5),d4s2(9))
         call swap(d4s2(7),d4s2(12))
      endif
      if (type .eq. 2) then
         call swap(d1s1(2),d1s1(3))
         call swap(d2s1(2),d2s1(3))
         call swap(d2s1(4),d2s1(6))
         call swap(d2s2(2),d2s2(3))
         call swap(d2s2(4),d2s2(6))
         call swap(d3s2(2),d3s2(3))
         call swap(d3s2(4),d3s2(6))
         call swap(d3s2(7),d3s2(10))
         call swap(d3s2(8),d3s2(9))
         call swap(d4s2(2),d4s2(3))
         call swap(d4s2(4),d4s2(6))
         call swap(d4s2(7),d4s2(10))
         call swap(d4s2(8),d4s2(9))
         call swap(d4s2(11),d4s2(15))
         call swap(d4s2(12),d4s2(14))
      endif
 
      xs = dsign(1.d0,xc)
      ys = dsign(1.d0,yc)
      zs = dsign(1.d0,zc)
      box1 = 1.d0/box
      box2 = box1*box1
      box3 = box2*box1
 
      d0s1 = d0s1*box1
 
      d1s1(1) = d1s1(1)*xs*box2
      d1s1(2) = d1s1(2)*ys*box2
      d1s1(3) = d1s1(3)*zs*box2
 
      d2s1(1) = d2s1(1)*box3
      d2s1(2) = d2s1(2)*xs*ys*box3
      d2s1(3) = d2s1(3)*xs*zs*box3
      d2s1(4) = d2s1(4)*box3
      d2s1(5) = d2s1(5)*ys*zs*box3
      d2s1(6) = d2s1(6)*box3
 
      d2s2(1) = d2s2(1)*box1
      d2s2(2) = d2s2(2)*xs*ys*box1
      d2s2(3) = d2s2(3)*xs*zs*box1
      d2s2(4) = d2s2(4)*box1
      d2s2(5) = d2s2(5)*ys*zs*box1
      d2s2(6) = d2s2(6)*box1
 
      d3s2(1) = d3s2(1)*xs*box2
      d3s2(2) = d3s2(2)*ys*box2
      d3s2(3) = d3s2(3)*zs*box2
      d3s2(4) = d3s2(4)*xs*box2
      d3s2(5) = d3s2(5)*xs*ys*zs*box2
      d3s2(6) = d3s2(6)*xs*box2
      d3s2(7) = d3s2(7)*ys*box2
      d3s2(8) = d3s2(8)*zs*box2
      d3s2(9) = d3s2(9)*ys*box2
      d3s2(10) = d3s2(10)*zs*box2
 
      d4s2(1) = d4s2(1)*box3
      d4s2(2) = d4s2(2)*xs*ys*box3
      d4s2(3) = d4s2(3)*xs*zs*box3
      d4s2(4) = d4s2(4)*box3
      d4s2(5) = d4s2(5)*ys*zs*box3
      d4s2(6) = d4s2(6)*box3
      d4s2(7) = d4s2(7)*xs*ys*box3
      d4s2(8) = d4s2(8)*xs*zs*box3
      d4s2(9) = d4s2(9)*xs*ys*box3
      d4s2(10) = d4s2(10)*xs*zs*box3
      d4s2(11) = d4s2(11)*box3
      d4s2(12) = d4s2(12)*ys*zs*box3
      d4s2(13) = d4s2(13)*box3
      d4s2(14) = d4s2(14)*ys*zs*box3
      d4s2(15) = d4s2(15)*box3
 
      f1 = 8.37758040957278196917d0*box3
      f2 = 2.51327412287183459075d0*box3
      f21 = 0.5d0*f2
      f22 = -f2/3.d0
      f23 = -2.d0*f22

      g0s0s(1,1) = d0s1-d2s2(6)
      g0s0s(1,2) = -d2s2(5)
      g0s0s(1,3) = -d2s2(3)
      g0s0s(2,1) = -d2s2(5)
      g0s0s(2,2) = d0s1-d2s2(4)
      g0s0s(2,3) = -d2s2(2)
      g0s0s(3,1) = -d2s2(3)
      g0s0s(3,2) = -d2s2(2)
      g0s0s(3,3) = d0s1-d2s2(1)
 
      g0s1s(1,1) = 0.5d0*d1s1(2)-d3s2(9)
      g0s1s(1,2) = -d3s2(8)
      g0s1s(1,3) = 0.5d0*d1s1(1)-d3s2(6)
      g0s1s(1,4) = -d3s2(5)
      g0s1s(1,5) = -d3s2(3)
      g0s1s(2,1) = 0.5d0*d1s1(3)-d3s2(8)
      g0s1s(2,2) = d1s1(2)-d3s2(7)
      g0s1s(2,3) = -d3s2(5)
      g0s1s(2,4) = 0.5d0*d1s1(1)-d3s2(4)
      g0s1s(2,5) = -d3s2(2)
      g0s1s(3,1) = -d3s2(5)
      g0s1s(3,2) = -d3s2(4)
      g0s1s(3,3) = 0.5d0*d1s1(3)-d3s2(3)
      g0s1s(3,4) = 0.5d0*d1s1(2)-d3s2(2)
      g0s1s(3,5) = d1s1(1)-d3s2(1)
 
      g0s2s(1,1) = d2s1(4)/3.d0+d2s1(6)/15.d0-d4s2(13)
      g0s2s(1,2) = 0.2d0*d2s1(5)-d4s2(12)
      g0s2s(1,3) = d2s1(2)/3.d0-d4s2(9)
      g0s2s(1,4) = d2s1(3)/15.d0-d4s2(8)
      g0s2s(1,5) = d2s1(1)/3.d0+d2s1(6)/15.d0-d4s2(6)
      g0s2s(1,6) = d2s1(5)/15.d0-d4s2(5)
      g0s2s(1,7) = 0.2d0*d2s1(3)-d4s2(3)
      g0s2s(2,1) = 11.d0*d2s1(5)/15.d0-d4s2(12)
      g0s2s(2,2) = 1.2d0*d2s1(4)-d4s2(11)
      g0s2s(2,3) = d2s1(3)/3.d0-d4s2(8)
      g0s2s(2,4) = 11.d0*d2s1(2)/15.d0-d4s2(7)
      g0s2s(2,5) = d2s1(5)/15.d0-d4s2(5)
      g0s2s(2,6) = d2s1(1)/3.d0+d2s1(4)/15.d0-d4s2(4)
      g0s2s(2,7) = 0.2d0*d2s1(2)-d4s2(2)
      g0s2s(3,1) = d2s1(3)/15.d0-d4s2(8)
      g0s2s(3,2) = 0.2d0*d2s1(2)-d4s2(7)
      g0s2s(3,3) = d2s1(5)/3.d0-d4s2(5)
      g0s2s(3,4) = d2s1(4)/3.d0+d2s1(1)/15.d0-d4s2(4)
      g0s2s(3,5) = 11.d0*d2s1(3)/15.d0-d4s2(3)
      g0s2s(3,6) = 11.d0*d2s1(2)/15.d0-d4s2(2)
      g0s2s(3,7) = 1.2d0*d2s1(1)-d4s2(1)
 
      g1s1s(1,1) = 0.25d0*(d2s1(4)+d2s1(6))-d4s2(13)+f21
      g1s1s(1,2) = 0.5d0*d2s1(5)-d4s2(12)
      g1s1s(1,3) = 0.25d0*d2s1(2)-d4s2(9)
      g1s1s(1,4) = 0.25d0*d2s1(3)-d4s2(8)
      g1s1s(1,5) = -d4s2(5)
      g1s1s(2,1) = g1s1s(1,2)
      g1s1s(2,2) = d2s1(4)-d4s2(11)+f23
      g1s1s(2,3) = -d4s2(8)
      g1s1s(2,4) = 0.5d0*d2s1(2)-d4s2(7)
      g1s1s(2,5) = -d4s2(4)+f22
      g1s1s(3,1) = g1s1s(1,3)
      g1s1s(3,2) = g1s1s(2,3)
      g1s1s(3,3) = 0.25d0*(d2s1(1)+d2s1(6))-d4s2(6)+f21
      g1s1s(3,4) = 0.25d0*d2s1(5)-d4s2(5)
      g1s1s(3,5) = 0.5d0*d2s1(3)-d4s2(3)
      g1s1s(4,1) = g1s1s(1,4)
      g1s1s(4,2) = g1s1s(2,4)
      g1s1s(4,3) = g1s1s(3,4)
      g1s1s(4,4) = 0.25d0*(d2s1(1)+d2s1(4))-d4s2(4)+f21
      g1s1s(4,5) = 0.5d0*d2s1(2)-d4s2(2)
      g1s1s(5,1) = g1s1s(1,5)
      g1s1s(5,2) = g1s1s(2,5)
      g1s1s(5,3) = g1s1s(3,5)
      g1s1s(5,4) = g1s1s(4,5)
      g1s1s(5,5) = d2s1(1)-d4s2(1)+f23
 
      g0s0t(1,1) = 0.d0
      g0s0t(1,2) = d1s1(1)
      g0s0t(1,3) = -d1s1(2)
      g0s0t(2,1) = -g0s0t(1,2)
      g0s0t(2,2) = 0.d0
      g0s0t(2,3) = d1s1(3)
      g0s0t(3,1) = -g0s0t(1,3)
      g0s0t(3,2) = -g0s0t(2,3)
      g0s0t(3,3) = 0.d0
 
      g0s1t(1,1) = 0.5d0*d2s1(3)
      g0s1t(1,2) = d2s1(2)
      g0s1t(1,3) = -0.5d0*d2s1(5)
      g0s1t(1,4) = 0.5d0*(d2s1(1)-d2s1(4))
      g0s1t(1,5) = -d2s1(2)
      g0s1t(2,1) = -0.5d0*d2s1(2)
      g0s1t(2,2) = 0.d0
      g0s1t(2,3) = 0.5d0*(d2s1(6)-d2s1(1))
      g0s1t(2,4) = 0.5d0*d2s1(5)
      g0s1t(2,5) = d2s1(3)
      g0s1t(3,1) = 0.5d0*(d2s1(4)-d2s1(6))
      g0s1t(3,2) = -d2s1(5)
      g0s1t(3,3) = 0.5d0*d2s1(2)
      g0s1t(3,4) = -0.5d0*d2s1(3)
      g0s1t(3,5) = 0.d0
 
      g0s0p(1,1) = -d2s1(6)+f1
      g0s0p(1,2) = -d2s1(5)
      g0s0p(1,3) = -d2s1(3)
      g0s0p(2,1) = -d2s1(5)
      g0s0p(2,2) = -d2s1(4)+f1
      g0s0p(2,3) = -d2s1(2)
      g0s0p(3,1) = -d2s1(3)
      g0s0p(3,2) = -d2s1(2)
      g0s0p(3,3) = -d2s1(1)+f1
 
      s1v = d0s1
 
      end
