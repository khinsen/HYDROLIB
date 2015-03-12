************************************
* Friction and mobility            *
* for a hard-sphere configuration  *
* Module LUBRICATION               *
*                                  *
* K. Hinsen                        *
* Last revision: July 4, 1994      *
************************************


* Load tables for short-range contributions
      
      subroutine luload(fname)
      implicit none
      character*150 fname

      integer npoints
      parameter (npoints = 1000)
      double precision rv(0:npoints)
      double precision x11a(0:npoints),x12a(0:npoints)
      double precision y11a(0:npoints),y12a(0:npoints)
      double precision y11b(0:npoints),y12b(0:npoints)
      double precision x11c(0:npoints),x12c(0:npoints)
      double precision y11c(0:npoints),y12c(0:npoints)
      common /tab/ rv,x11a,x12a,y11a,y12a,y11b,y12b,x11c,x12c,y11c,y12c
      
      integer i
      
      open (unit = 20, file = fname, status = 'UNKNOWN')
      read (20,*) i
      if (i .ne. npoints+1) stop 'Wrong table file!'
      do 300 i = 0,npoints
         read (20,*) rv(i)
         rv(i) = 2.d0*rv(i)
300   continue
      do 301 i = 0,npoints
         read (20,*) x11a(i),x12a(i)
301   continue
      do 302 i = 0,npoints
         read (20,*) y11a(i),y12a(i)
302   continue
      do 303 i = 0,npoints
         read (20,*) x11c(i),x12c(i)
303   continue
      do 304 i = 0,npoints
         read (20,*) y11c(i),y12c(i)
304   continue
      do 305 i = 0,npoints
         read (20,*) y11b(i),y12b(i)
305   continue
      close(20)

      end


* Calculate short-range contribution
      
      subroutine lub(c,cnct,fr,np,nfixed,rigid,lm)
      implicit none
      integer np,nfixed,lm
      double precision c(0:2,np),fr(6*np,6*np)
      character*1 cnct(np,np)
      logical rigid

#ifdef PERIODIC      
      common /pbc/ box
      double precision box
#endif

      integer npoints
      parameter (npoints = 1000)
      double precision rv(0:npoints)
      double precision x11a(0:npoints),x12a(0:npoints)
      double precision y11a(0:npoints),y12a(0:npoints)
      double precision y11b(0:npoints),y12b(0:npoints)
      double precision x11c(0:npoints),x12c(0:npoints)
      double precision y11c(0:npoints),y12c(0:npoints)
      common /tab/ rv,x11a,x12a,y11a,y12a,y11b,y12b,x11c,x12c,y11c,y12c
      
      double precision dr,r(3),ra,f,x1,x2,lx1,lx2,f1,f2
      double precision xa(2),ya(2),yb(2),xc(2),yc(2)
      double precision ftt(3,3,2),frr(3,3,2),frt(3,3,2)
      integer i,j,k,l,m
      logical ndivp

      double precision cutoff(0:3)
      data cutoff /4.d0,3.8d0,3.d0,2.6d0/
      
      dr = rv(1)-rv(0)
      do 100 i = 0,np-1
         do 110 j = i+1,np-1
            r(1) = c(0,j+1)-c(0,i+1)
            r(2) = c(1,j+1)-c(1,i+1)
            r(3) = c(2,j+1)-c(2,i+1)
#ifdef PERIODIC
            if (r(1) .gt. 0.5d0*box) r(1) = r(1)-box
            if (r(1) .lt. -0.5d0*box) r(1) = r(1)+box
            if (r(2) .gt. 0.5d0*box) r(2) = r(2)-box
            if (r(2) .lt. -0.5d0*box) r(2) = r(2)+box
            if (r(3) .gt. 0.5d0*box) r(3) = r(3)-box
            if (r(3) .lt. -0.5d0*box) r(3) = r(3)+box
#endif
            ndivp = .false.
            if (rigid) ndivp = cnct(i+1,j+1) .ne. ' '
            if (nfixed .gt. 0) ndivp = ndivp .or.
     .           ((i+1 .gt. np-nfixed) .and. (j+1 .gt. np-nfixed))
            ra = dsqrt(r(1)*r(1)+r(2)*r(2)+r(3)*r(3))
            if (ra .le. cutoff(lm)) then
               r(1) = r(1)/ra
               r(2) = r(2)/ra
               r(3) = r(3)/ra
               k = idint((ra-rv(0))/dr)
               f = (ra-rv(k))/dr
               if (ndivp) then
                  f1 = 0.d0
               else
                  x1 = 1.d0-2.d0/ra
                  lx1 = log(x1)
                  f1 = 0.125d0/x1 - 9.d0*lx1/40.d0 - 9.d0*x1*lx1/168.d0
               endif
               x2 = 1.d0+2.d0/ra
               lx2 = log(x2)
               f2 = 0.125d0/x2 - 9.d0*lx2/40.d0 - 9.d0*x2*lx2/168.d0
               xa(1) = f2+f1 + x11a(k) + f*(x11a(k+1)-x11a(k))
               xa(2) = f2-f1 + x12a(k) + f*(x12a(k+1)-x12a(k))
               if (ndivp) then
                  f1 = 0.d0
               else
                  f1 = -lx1/6.d0
               endif
               f2 = -lx2/6.d0
               ya(1) = f2+f1 + y11a(k) + f*(y11a(k+1)-y11a(k))
               ya(2) = f2-f1 + y12a(k) + f*(y12a(k+1)-y12a(k))
               if (ndivp) then
                  f1 = 0.d0
               else
                  f1 = 0.25d0*lx1*(1.d0+x1)
               endif
               f2 = -0.25d0*lx2*(1.d0+x2)
               yb(1) = f2+f1 + y11b(k) + f*(y11b(k+1)-y11b(k))
               yb(2) = f2-f1 + y12b(k) + f*(y12b(k+1)-y12b(k))
               if (ndivp) then
                  f1 = 0.d0
               else
                  f1 = 0.25d0*x1*lx1
               endif
               f2 = 0.25d0*x2*lx2
               xc(1) = f2+f1 + x11c(k) + f*(x11c(k+1)-x11c(k))
               xc(2) = f2-f1 + x12c(k) + f*(x12c(k+1)-x12c(k))
               if (ndivp) then
                  yc(1) = -0.2d0*lx2-47.d0*x2*lx2/125.d0
     .                 + y11c(k) + f*(y11c(k+1)-y11c(k))
                  yc(2) = 0.05d0*lx2+31.d0*x2*lx2/250.d0
     .                 + y12c(k) + f*(y12c(k+1)-y12c(k))
               else
                  yc(1) = -0.2d0*(lx1+lx2)-47.d0*(x1*lx1+x2*lx2)/125.d0
     .                 + y11c(k) + f*(y11c(k+1)-y11c(k))
                  yc(2) = -0.05d0*(lx1-lx2)-31.d0*(x1*lx1-x2*lx2)/250.d0
     .                 + y12c(k) + f*(y12c(k+1)-y12c(k))
               endif
               do 120 k = 1,2
                  do 121 l = 1,3
                     do 122 m = 1,3
                        ftt(l,m,k) = 1.5d0*(xa(k)-ya(k))*r(l)*r(m)
                        frr(l,m,k) = 2.d0*(xc(k)-yc(k))*r(l)*r(m)
122                  continue
                     ftt(l,l,k) = ftt(l,l,k)+1.5d0*ya(k)
                     frr(l,l,k) = frr(l,l,k)+2.d0*yc(k)
                     frt(l,l,k) = 0.d0
121               continue
                  frt(1,2,k) = yb(k)*r(3)
                  frt(2,3,k) = yb(k)*r(1)
                  frt(3,1,k) = yb(k)*r(2)
                  frt(2,1,k) = -frt(1,2,k)
                  frt(3,2,k) = -frt(2,3,k)
                  frt(1,3,k) = -frt(3,1,k)
120            continue
               do 130 l = 1,3
                  do 131 m = 1,3
                     fr(i*6+l,j*6+m) = fr(i*6+l,j*6+m)
     .                    + ftt(l,m,2)
                     fr(j*6+m,i*6+l) = fr(j*6+m,i*6+l)
     .                    + ftt(l,m,2)
                     fr(i*6+l,i*6+m) = fr(i*6+l,i*6+m)
     .                    + ftt(l,m,1)
                     fr(j*6+m,j*6+l) = fr(j*6+m,j*6+l)
     .                    + ftt(l,m,1)
                     fr(i*6+l+3,j*6+m+3) = fr(i*6+l+3,j*6+m+3)
     .                    + frr(l,m,2)
                     fr(j*6+m+3,i*6+l+3) = fr(j*6+m+3,i*6+l+3)
     .                    + frr(l,m,2)
                     fr(i*6+l+3,i*6+m+3) = fr(i*6+l+3,i*6+m+3)
     .                    + frr(l,m,1)
                     fr(j*6+m+3,j*6+l+3) = fr(j*6+m+3,j*6+l+3)
     .                    + frr(l,m,1)

                     fr(i*6+l+3,j*6+m) = fr(i*6+l+3,j*6+m)
     .                    + frt(l,m,2)
                     fr(j*6+l+3,i*6+m) = fr(j*6+l+3,i*6+m)
     .                    - frt(l,m,2)
                     fr(i*6+l+3,i*6+m) = fr(i*6+l+3,i*6+m)
     .                    + frt(l,m,1)
                     fr(j*6+l+3,j*6+m) = fr(j*6+l+3,j*6+m)
     .                    - frt(l,m,1)

                     fr(j*6+m,i*6+l+3) = fr(j*6+m,i*6+l+3)
     .                    - frt(m,l,2)
                     fr(i*6+m,j*6+l+3) = fr(i*6+m,j*6+l+3)
     .                    + frt(m,l,2)
                     fr(i*6+m,i*6+l+3) = fr(i*6+m,i*6+l+3)
     .                    - frt(m,l,1)
                     fr(j*6+m,j*6+l+3) = fr(j*6+m,j*6+l+3)
     .                    + frt(m,l,1)
131               continue
130            continue
            endif
110      continue
100   continue
      
      end
