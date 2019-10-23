#ifdef DPSOURCE
#include <gcalc_dp.f>
#else
************************************
* Friction and mobility            *
* for a hard-sphere configuration  *
* Module GCALC                     *
*                                  *
* K. Hinsen                        *
* Last revision: December 7, 1993  *
************************************

#include <filenames.h>

#ifdef PERIODIC
#include <pbc.f>
#endif

* Load coefficients for non-periodic interactions

      subroutine ocload
      implicit none

      integer lm
      parameter (lm = _LM_)

      integer llen
      parameter (llen = 3000)
      integer g1s(0:lm,0:lm,2*lm+3,2*lm+3),g1e(0:lm,0:lm,2*lm+3,2*lm+3)
      integer g2s(0:lm,0:lm,2*lm+3,2*lm+3),g2e(0:lm,0:lm,2*lm+3,2*lm+3)
      integer g3s(2:2*lm+2,4*lm+5),g3e(2:2*lm+2,4*lm+5)
      real cf(llen)
      integer xp(llen),yp(llen),zp(llen)
      common /oc/ g1s,g1e,g2s,g2e,g3s,g3e,cf,xp,yp,zp

      integer ocf
      parameter (ocf = 10)

      integer l
      integer i,j,n,n1,n2,m

      open (unit=ocf,file=
     .FILENAME1
     .     ,form='FORMATTED',status='OLD')
      read (ocf,*) l
      if (l .lt. lm) stop 'Wrong coefficient file!'

      do 100 i = 0,l
         do 110 j = i,l
            read (ocf,*) n1,n2
            if (n1 .ne. i .or. n2 .ne. j) stop
            if (i.le.lm .and. j.le.lm) then
               do 200 n = 1,2*i+3
                  read (ocf,*) (g1s(i,j,n,m), m = 1,2*j+3)
200            continue
               do 210 n = 1,2*i+3
                  read (ocf,*) (g1e(i,j,n,m), m = 1,2*j+3)
210            continue
               do 220 n = 1,2*i+3
                  read (ocf,*) (g2s(i,j,n,m), m = 1,2*j+3)
220            continue
               do 230 n = 1,2*i+3
                  read (ocf,*) (g2e(i,j,n,m), m = 1,2*j+3)
230            continue
            else
               do 300 n = 1,2*i+3
                  read (ocf,*)
                  read (ocf,*)
                  read (ocf,*)
                  read (ocf,*)
300            continue
            endif
110      continue
100   continue

      n = 0
      do 400 i = 2,2*l+2
         read (ocf,*) n1
         if (n1 .ne. i) stop
         if (i .le. 2*lm+2) then
            do 410 j = 1,2*n1+1
               read (ocf,*) g3s(i,j)
410         continue
            do 420 j = 1,2*n1+1
               read (ocf,*) g3e(i,j)
               n = max(n,g3e(i,j))
420         continue
         else
            do 430 j = 1,2*n1+1
               read (ocf,*)
               read (ocf,*)
430         continue
         endif
400   continue

      read (ocf,*) n1
      if (n1 .gt. llen) stop 'Array size too small!'
      if (n .gt. n1) stop
      do 500 i = 1,n
         read (ocf,*) cf(i),xp(i),yp(i),zp(i)
500   continue

      close (ocf)

      end


* Calculate double factorial (!!)

      real function dfact(n)
      implicit none
      integer n

      integer i

      dfact = 1.
      do 100 i = 3,n,2
         dfact = dfact*i
100   continue

      end


* Determine correction matrix to compensate for eliminated components

      subroutine subst(trmc,ncc,v1,v2,dxc,dyc,dzc,ccnt)
      implicit none
      integer ncc,ccnt,v1,v2,dxc(ccnt),dyc(ccnt),dzc(ccnt)
      real trmc(ncc,ncc)

      integer v11,v12,v21,v22
      integer i

      v11 = 0
      v12 = 0
      v21 = 0
      v22 = 0
      do 100 i = 1,ccnt
         if (dzc(i) .eq. dzc(v1)-2) then
            if (dxc(i).eq.dxc(v1)+2 .and. dyc(i).eq.dyc(v1)) v11 = i
            if (dxc(i).eq.dxc(v1) .and. dyc(i).eq.dyc(v1)+2) v12 = i
         endif
         if (dzc(i) .eq. dzc(v2)-2) then
            if (dxc(i).eq.dxc(v2)+2 .and. dyc(i).eq.dyc(v2)) v21 = i
            if (dxc(i).eq.dxc(v2) .and. dyc(i).eq.dyc(v2)+2) v22 = i
         endif
100   continue
      if (v11.eq.0 .or. v12.eq.0 .or. v21.eq.0 .or. v22.eq.0)
     .         stop 'Error 1'
      trmc(v11,v21) = trmc(v11,v21) + trmc(v1,v2)
      trmc(v11,v22) = trmc(v11,v22) + trmc(v1,v2)
      trmc(v12,v21) = trmc(v12,v21) + trmc(v1,v2)
      trmc(v12,v22) = trmc(v12,v22) + trmc(v1,v2)

      end


      subroutine elim(trmc,ncc,var,dxc,dyc,dzc,ccnt)
      implicit none
      integer ncc,var,ccnt,dxc(ccnt),dyc(ccnt),dzc(ccnt)
      real trmc(ncc,ncc)

      integer i,j

      do 100 i = var,ccnt-1
         do 110 j = 1,ccnt
            trmc(i,j) = trmc(i+1,j)
110      continue
100   continue
      do 120 i = var,ccnt-1
         do 130 j = 1,ccnt
            trmc(j,i) = trmc(j,i+1)
130      continue
         dxc(i) = dxc(i+1)
         dyc(i) = dyc(i+1)
         dzc(i) = dzc(i+1)
120   continue
      ccnt = ccnt - 1

      end


      subroutine msqrt(m,nm,dm)
      implicit none
      integer nm,dm
      real m(dm,dm)

      integer nmax
      parameter (nmax = 11)

      integer i,j,k
#ifdef IMSL
      real eval(nmax),evec(nmax,nmax),ievec(nmax,nmax)
#endif
#ifdef LAPACK
      integer lwork
      parameter(lwork = 3*nmax)
      integer ipiv(nmax)
      real eval(nmax),ievec(nmax,nmax),evec(nmax,nmax)
      real work(lwork),info
#endif

      if (nm .gt. nmax) stop 'Error in msqrt!'

#ifdef IMSL
#ifdef DP
      call devcsf(nm,m,dm,eval,evec,nmax)
      call dlinrg(nm,evec,nmax,ievec,nmax)
#else
      call evcsf(nm,m,dm,eval,evec,nmax)
      call linrg(nm,evec,nmax,ievec,nmax)
#endif
#endif
#ifdef LAPACK
#ifdef DP
      call dsyev('V','U',nm,m,dm,eval,work,lwork,info)
#else
      call ssyev('V','U',nm,m,dm,eval,work,lwork,info)
#endif
      if (info .ne. 0) stop 'Error in eigenvalue calculation!'
      do 10 i = 1,nm
         do 11 j = 1,nm
            evec(i,j) = m(i,j)
11       continue
10    continue
#ifdef DP
      call dgetrf(nm,nm,m,dm,ipiv,info)
#else
      call sgetrf(nm,nm,m,dm,ipiv,info)
#endif
      if (info .ne. 0) stop 'Error in factorization!'
#ifdef DP
      call dgetri(nm,m,dm,ipiv,work,lwork,info)
#else
      call sgetri(nm,m,dm,ipiv,work,lwork,info)
#endif
      if (info .ne. 0) stop 'Error in inversion!'
      do 20 i = 1,nm
         do 21 j = 1,nm
            ievec(i,j) = m(i,j)
21       continue
20    continue
#endif

      do 100 i = 1,nm
         eval(i) = sqrt(eval(i))
100   continue

      do 120 i = 1,nm
         do 121 j = 1,nm
            m(i,j) = 0.
            do 122 k = 1,nm
               m(i,j) = m(i,j) + evec(i,k)*eval(k)*ievec(k,j)
122         continue
121      continue
120   continue

      end


      subroutine comp(xi)
      implicit none
      real xi

      real dfact
      external dfact

      integer lm
      parameter (lm = _LM_)

      integer ttype(3*lm+3),trank(3*lm+3),tindex(3*lm+3)
      common /tens/ ttype,trank,tindex

      integer im
      parameter (im = lm+1)
      integer nc
      parameter (nc = 2*im+1)
      integer dx(2*im,4*im+1),dy(2*im,4*im+1),dz(2*im,4*im+1)
      real trm(im,nc,nc),trms(im,nc,nc)
      common /tdf/ dx,dy,dz,trm,trms
 
      integer ncc
      parameter (ncc = (im+1)*(im+2)/2)
      integer dxc(ncc),dyc(ncc),dzc(ncc)
      real trmc(ncc,ncc)

      real z(0:lm,4),zs(0:lm,4)
      common /zm/ z,zs

      integer i,j,k,l,lc,n,idx,ccnt,max
      real fi,fj,fk,fl
      real m(2,2)
      logical ef(ncc)

      idx = 0
      do 10 i = 0,lm
         do 20 j = 0,2
            ttype(3*i+j+1) = j
            trank(3*i+j+1) = i+1
            tindex(3*i+j+1) = idx
            idx = idx + 2*i+3
20       continue
10    continue

      fl = 1
      do 100 l = 1,im
         lc = (l+1)*(l+2)/2
         do 101 i = 1,lc
            do 102 j = 1,lc
               trmc(i,j) = 0
102         continue
101      continue
         idx = 1
         fl = l*fl
         fi = 1
         do 110 i = 0,l
            if (i .gt. 1) fi = i*fi
            fj = 1
            do 120 j = 0,l-i
               if (j .gt. 1) fj = j*fj
               k = l-i-j
               fk = 1.
               do 130 n = 2,k
                  fk = n*fk
130            continue
               dxc(idx) = i
               dyc(idx) = j
               dzc(idx) = k
               trmc(idx,idx) = fl/(fi*fj*fk)
               idx = idx + 1
120         continue
110      continue
         if (idx-1 .ne. lc) stop 'Error 2'
         ccnt = lc
300      max = 1
         do 310 i = 2,ccnt
            if (dzc(i) .gt. dzc(max)) max = i
310      continue
         if (dzc(max) .lt. 2) goto 399
         do 320 i = 1,ccnt
            ef(i) = trmc(i,max) .ne. 0
320      continue
         do 330 i = 1,ccnt
            if (ef(i)) then
               do 340 j = 1,ccnt
                  if (ef(j)) call subst(trmc,ncc,i,j,dxc,dyc,dzc,ccnt)
340            continue
            endif
330      continue
         do 350 i = ccnt,1,-1
            if (ef(i)) call elim(trmc,ncc,i,dxc,dyc,dzc,ccnt)
350      continue
         goto 300
399      if (ccnt .ne. 2*l+1) stop 'Error 3'
         do 400 i = 1,nc
            dx(l,i) = dxc(i)
            dy(l,i) = dyc(i)
            dz(l,i) = dzc(i)
            do 410 j = 1,2*l+1
               trm(l,i,j) = trmc(i,j)
410         continue
400      continue
         call msqrt(trmc,2*l+1,ncc)
         do 420 i = 1,nc
            do 430 j = 1,2*l+1
               trms(l,i,j) = trmc(i,j)
430         continue
420      continue
100   continue

      do 500 l = im,2*im
         idx = 1
         do 510 i = 0,l
            do 520 j = 0,l-i
               k = l-i-j
               if (k .lt. 2) then
                  dx(l,idx) = i
                  dy(l,idx) = j
                  dz(l,idx) = k
                  idx = idx+1
               endif
520         continue
510      continue
500   continue

      fl = 1.
      do 600 l = 0,lm
         if (l .gt. 0) fl = l*fl
         z(l,1) = -(2.*l+3.)*(1.-xi)/
     .        (fl*(1.+2.*l*xi)*(l+2.)*dfact(2*l-1))
         z(l,2) = -(1.-(l+3.)*xi)/(fl*(1.+l*xi)*(l+2.)*dfact(2*l+1))
         z(l,3) = -(2.*l+3.)**2*(1.-5.*xi)/
     .        (4.*fl*(1.+2.*l*xi)*(l+2.)*dfact(2*l+5))
         z(l,4) = -(1.-3.*xi)/(2.*fl*(l+2.)*(1.+2.*l*xi)*dfact(2*l-1))
600   continue

      do 610 l = 0,lm
         m(1,1) = -z(l,1)
         m(1,2) = -z(l,4)
         m(2,1) = -z(l,4)
         m(2,2) = -z(l,3)
         call msqrt(m,2,2)
         zs(l,1) = m(1,1)
         zs(l,3) = m(2,2)
         zs(l,4) = m(1,2)
         zs(l,2) = sqrt(-z(l,2))
610   continue

      end


************************************
* Calculate the interaction matrix *
************************************

* Calculate powers of r and x/r

      subroutine powers(x,y,z)
      implicit none
      real x,y,z

      integer lm
      parameter (lm = _LM_)

      real px(3,0:2*lm+2),pr(0:2*lm+3)
      common /po/ px,pr

      integer i
      real r,xh,yh,zh

      r = sqrt(x**2+y**2+z**2)
      pr(0) = 1.
      do 100 i = 1,2*lm+3
         pr(i) = r*pr(i-1)
100   continue
      xh = x/r
      yh = y/r
      zh = z/r
      px(1,0) = 1.
      px(2,0) = 1.
      px(3,0) = 1.
      do 110 i = 1,2*lm+2
         px(1,i) = xh*px(1,i-1)
         px(2,i) = yh*px(2,i-1)
         px(3,i) = zh*px(3,i-1)
110   continue

      end


* Evaluate terms for one component

      real function ev(si,ei)
      implicit none
      integer si,ei

      integer lm
      parameter (lm = _LM_)

      integer llen
      parameter (llen = 3000)
      integer g1s(0:lm,0:lm,2*lm+3,2*lm+3),g1e(0:lm,0:lm,2*lm+3,2*lm+3)
      integer g2s(0:lm,0:lm,2*lm+3,2*lm+3),g2e(0:lm,0:lm,2*lm+3,2*lm+3)
      integer g3s(2:2*lm+2,4*lm+5),g3e(2:2*lm+2,4*lm+5)
      real cf(llen)
      integer xp(llen),yp(llen),zp(llen)
      common /oc/ g1s,g1e,g2s,g2e,g3s,g3e,cf,xp,yp,zp

      real px(3,0:2*lm+2),pr(0:2*lm+3)
      common /po/ px,pr

      integer i

      ev = 0.
      do 100 i = si,ei
         ev = ev + cf(i)*px(1,xp(i))*px(2,yp(i))*px(3,zp(i))
100   continue

      end


* Find coefficients for a term

      integer function find(x,y,z,n1,n2)
      implicit none
      integer x,y,z,n1,n2

      integer lm
      parameter (lm = _LM_)

      integer im
      parameter (im = lm+1)
      integer nc
      parameter (nc = 2*im+1)
      integer dx(2*im,4*im+1),dy(2*im,4*im+1),dz(2*im,4*im+1)
      real trm(im,nc,nc),trms(im,nc,nc)
      common /tdf/ dx,dy,dz,trm,trms

      integer i

      do 100 i = 1,2*(n1+n2)+1
         if (x.eq.dx(n1+n2,i) .and. y.eq.dy(n1+n2,i)
     .               .and. z.eq.dz(n1+n2,i)) then
            find = i
            return
         endif
100   continue
      stop 'Not found!'

      end


* Calculate submatrices of G

      subroutine gss(g,ng,n1,n2,first)
      implicit none
      integer ng,n1,n2
      real g(ng,ng)
      logical first

      integer lm
      parameter (lm = _LM_)

      integer llen
      parameter (llen = 3000)
      integer g1s(0:lm,0:lm,2*lm+3,2*lm+3),g1e(0:lm,0:lm,2*lm+3,2*lm+3)
      integer g2s(0:lm,0:lm,2*lm+3,2*lm+3),g2e(0:lm,0:lm,2*lm+3,2*lm+3)
      integer g3s(2:2*lm+2,4*lm+5),g3e(2:2*lm+2,4*lm+5)
      real cf(llen)
      integer xp(llen),yp(llen),zp(llen)
      common /oc/ g1s,g1e,g2s,g2e,g3s,g3e,cf,xp,yp,zp

      integer im
      parameter (im = lm+1)
      integer nc
      parameter (nc = 2*im+1)
      integer dx(2*im,4*im+1),dy(2*im,4*im+1),dz(2*im,4*im+1)
      real trm(im,nc,nc),trms(im,nc,nc)
      common /tdf/ dx,dy,dz,trm,trms

      real px(3,0:2*lm+2),pr(0:2*lm+3)
      common /po/ px,pr

#ifdef PERIODIC
      real g0s0s(3,3),g0s1s(3,5),g0s2s(3,7),g1s1s(5,5)
      real g0s0t(3,3),g0s1t(3,5),g0s0p(3,3),s1v
      common /gpbca/ g0s0s,g0s1s,g0s2s,g1s1s,g0s0t,g0s1t,g0s0p,s1v
#endif

      real dfact,ev
      external dfact,ev

      integer i,j
      real f

#ifdef PERIODIC
      if (n1+n2 .le. 2) then
         if (n1 .eq. 0 .and. n2 .eq. 0) then
            do 10 i = 1,3
               do 20 j = 1,3
                  if (first) then
                     g(i,j) = g0s0s(i,j)
                  else
                     g(i,j) = 0.
                  endif
20             continue
10          continue
         else if (n1 .eq. 0 .and. n2 .eq. 1) then
            do 30 i = 1,3
               do 40 j = 1,5
                  if (first) then
                     g(i,j) = g0s1s(i,j)
                  else
                     g(i,j) = 0.
                  endif
40             continue
30          continue
         else if (n1 .eq. 0 .and. n2 .eq. 2) then
            do 50 i = 1,3
               do 60 j = 1,7
                  if (first) then
                     g(i,j) = g0s2s(i,j)
                  else
                     g(i,j) = 0.
                  endif
60             continue
50          continue
         else if (n1 .eq. 1 .and. n2 .eq. 1) then
            do 70 i = 1,5
               do 80 j = 1,5
                  if (first) then
                     g(i,j) = g1s1s(i,j)
                  else
                     g(i,j) = 0.
                  endif
80             continue
70          continue
         endif
      else
#endif
         f = -0.5*pr(2)/(2*(n1+n2)+3)
         do 90 i = 1,2*n1+3
            do 91 j = 1,2*n2+3
               g(i,j) = f*g(i,j)
91          continue
90       continue

         f = (n1+n2-n1*n2+2)*dfact(2*(n1+n2)-1)
     .                     /((2*(n1+n2)+3)*pr(n1+n2+1))
         if (mod(n1+n2,2) .eq. 1) f = -f
         if (n1 .eq. n2) then
            do 100 i = 1,2*n1+3
               do 110 j = i,2*n2+3
                  g(i,j) = g(i,j) + f*ev(g1s(n1,n2,i,j),g1e(n1,n2,i,j))
110            continue
100         continue
            do 120 i = 2,2*n1+3
               do 130 j = 1,i-1
                  g(i,j) = g(j,i)
130            continue
120         continue
         else
            do 200 i = 1,2*n1+3
               do 210 j = 1,2*n2+3
                  if (n2 .ge. n1) then
                     g(i,j) = g(i,j)
     .                     + f*ev(g1s(n1,n2,i,j),g1e(n1,n2,i,j))
                  else
                     g(i,j) = g(i,j)
     .                     + f*ev(g1s(n2,n1,j,i),g1e(n2,n1,j,i))
                  endif
210            continue
200         continue
         endif
#ifdef PERIODIC
      endif
#endif

      end


      subroutine gst(g,ng,n1,n2,first)
      implicit none
      integer ng,n1,n2
      real g(ng,ng)
      logical first

      integer lm
      parameter (lm = _LM_)

      integer llen
      parameter (llen = 3000)
      integer g1s(0:lm,0:lm,2*lm+3,2*lm+3),g1e(0:lm,0:lm,2*lm+3,2*lm+3)
      integer g2s(0:lm,0:lm,2*lm+3,2*lm+3),g2e(0:lm,0:lm,2*lm+3,2*lm+3)
      integer g3s(2:2*lm+2,4*lm+5),g3e(2:2*lm+2,4*lm+5)
      real cf(llen)
      integer xp(llen),yp(llen),zp(llen)
      common /oc/ g1s,g1e,g2s,g2e,g3s,g3e,cf,xp,yp,zp

      integer im
      parameter (im = lm+1)
      integer nc
      parameter (nc = 2*im+1)
      integer dx(2*im,4*im+1),dy(2*im,4*im+1),dz(2*im,4*im+1)
      real trm(im,nc,nc),trms(im,nc,nc)
      common /tdf/ dx,dy,dz,trm,trms

      real px(3,0:2*lm+2),pr(0:2*lm+3)
      common /po/ px,pr

#ifdef PERIODIC
      real g0s0s(3,3),g0s1s(3,5),g0s2s(3,7),g1s1s(5,5)
      real g0s0t(3,3),g0s1t(3,5),g0s0p(3,3),s1v
      common /gpbca/ g0s0s,g0s1s,g0s2s,g1s1s,g0s0t,g0s1t,g0s0p,s1v
#endif

      real dfact,ev
      external dfact,ev

      integer i,j
      real f

#ifdef PERIODIC
      if (n1 .eq. 0 .and. n2 .eq. 0) then
         do 10 i = 1,3
            do 20 j = 1,3
               if (first) then
                  g(i,j) = g0s0t(i,j)
               else
                  g(i,j) = 0.
               endif
20          continue
10       continue
      else if (n1 .eq. 0 .and. n2 .eq. 1) then
         do 30 i = 1,3
            do 40 j = 1,5
               if (first) then
                  g(i,j) = g0s1t(i,j)
               else
                  g(i,j) = 0.
               endif
40          continue
30       continue
      else
#endif
         f = dfact(2*(n1+n2)+1)/pr(n1+n2+2)
         if (mod(n1+n2,2) .eq. 1) f = -f

         if (n1 .eq. n2) then
            do 100 i = 1,2*n1+3
               g(i,i) = 0.
               do 110 j = i+1,2*n2+3
                  g(i,j) = f*ev(g2s(n1,n2,i,j),g2e(n1,n2,i,j))
110            continue
100         continue
            do 120 i = 2,2*n1+3
               do 130 j = 1,i-1
                  g(i,j) = -g(j,i)
130            continue
120         continue
         else
            do 200 i = 1,2*n1+3
               do 210 j = 1,2*n2+3
                  if (n2 .ge. n1) then
                     g(i,j) = f*ev(g2s(n1,n2,i,j),g2e(n1,n2,i,j))
                  else
                     g(i,j) = -f*ev(g2s(n2,n1,j,i),g2e(n2,n1,j,i))
                  endif
210            continue
200         continue
         endif
#ifdef PERIODIC
      endif
#endif

      end


      subroutine gsp(g,ng,n1,n2,first)
      implicit none
      integer ng,n1,n2
      real g(ng,ng)
      logical first

      integer lm
      parameter (lm = _LM_)

      integer llen
      parameter (llen = 3000)
      integer g1s(0:lm,0:lm,2*lm+3,2*lm+3),g1e(0:lm,0:lm,2*lm+3,2*lm+3)
      integer g2s(0:lm,0:lm,2*lm+3,2*lm+3),g2e(0:lm,0:lm,2*lm+3,2*lm+3)
      integer g3s(2:2*lm+2,4*lm+5),g3e(2:2*lm+2,4*lm+5)
      real cf(llen)
      integer xp(llen),yp(llen),zp(llen)
      common /oc/ g1s,g1e,g2s,g2e,g3s,g3e,cf,xp,yp,zp

      integer im
      parameter (im = lm+1)
      integer nc
      parameter (nc = 2*im+1)
      integer dx(2*im,4*im+1),dy(2*im,4*im+1),dz(2*im,4*im+1)
      real trm(im,nc,nc),trms(im,nc,nc)
      common /tdf/ dx,dy,dz,trm,trms

      real px(3,0:2*lm+2),pr(0:2*lm+3)
      common /po/ px,pr

#ifdef PERIODIC
      real g0s0s(3,3),g0s1s(3,5),g0s2s(3,7),g1s1s(5,5)
      real g0s0t(3,3),g0s1t(3,5),g0s0p(3,3),s1v
      common /gpbca/ g0s0s,g0s1s,g0s2s,g1s1s,g0s0t,g0s1t,g0s0p,s1v
#endif

      integer find
      real dfact,ev
      external find,dfact,ev

      integer i,j,k,k1,k2
      integer ix,iy,iz
      real f

#ifdef PERIODIC
      if (n1 .eq. 0 .and. n2 .eq. 0) then
         do 10 i = 1,3
            do 20 j = 1,3
               if (first) then
                  g(i,j) = g0s0p(i,j)
               else
                  g(i,j) = 0.
               endif
20          continue
10       continue
      else
#endif
         f = dfact(2*(n1+n2)+3)/pr(n1+n2+3)
         if (mod(n1+n2,2) .eq. 0) f = -f

         if (n1 .eq. n2) then
            do 100 i = 1,2*n1+3
               do 110 j = i,2*n2+3
                  ix = dx(n1+1,i)+dx(n2+1,j)
                  iy = dy(n1+1,i)+dy(n2+1,j)
                  iz = dz(n1+1,i)+dz(n2+1,j)
                  if (iz .eq. 2) then
                     k1 = find(ix+2,iy,0,n1+1,n2+1)
                     k2 = find(ix,iy+2,0,n1+1,n2+1)
                     g(i,j) = -ev(g3s(n1+n2+2,k1),g3e(n1+n2+2,k1))
     .                             -ev(g3s(n1+n2+2,k2),g3e(n1+n2+2,k2))
                else
                     k = find(ix,iy,iz,n1+1,n2+1)
                     g(i,j) = ev(g3s(n1+n2+2,k),g3e(n1+n2+2,k))
                  endif
                  g(i,j) = f*g(i,j)
110            continue
100         continue
            do 130 i = 2,2*n1+3
               do 140 j = 1,i-1
                  g(i,j) = g(j,i)
140            continue
130         continue
         else
            do 200 i = 1,2*n1+3
               do 210 j = 1,2*n2+3
                  ix = dx(n1+1,i)+dx(n2+1,j)
                  iy = dy(n1+1,i)+dy(n2+1,j)
                  iz = dz(n1+1,i)+dz(n2+1,j)
                  if (iz .eq. 2) then
                     k1 = find(ix+2,iy,0,n1+1,n2+1)
                     k2 = find(ix,iy+2,0,n1+1,n2+1)
                     g(i,j) = -ev(g3s(n1+n2+2,k1),g3e(n1+n2+2,k1))
     .                             -ev(g3s(n1+n2+2,k2),g3e(n1+n2+2,k2))
                  else
                     k = find(ix,iy,iz,n1+1,n2+1)
                     g(i,j) = ev(g3s(n1+n2+2,k),g3e(n1+n2+2,k))
                  endif
                  g(i,j) = f*g(i,j)
210            continue
200         continue
         endif
#ifdef PERIODIC
      endif
#endif

      end


* Insert a contribution to G into the final matrix

      subroutine ins(g,dim,gx,ngx,n1,n2,t1,t2,f,tr)
      implicit none
      integer dim,ngx,n1,n2,t1,t2
      logical tr
      real g(dim,dim),gx(ngx,ngx),f

      integer lm
      parameter (lm = _LM_)

      integer ttype(3*lm+3),trank(3*lm+3),tindex(3*lm+3)
      common /tens/ ttype,trank,tindex

      integer im
      parameter (im = lm+1)
      integer nc
      parameter (nc = 2*im+1)
      integer dx(2*im,4*im+1),dy(2*im,4*im+1),dz(2*im,4*im+1)
      real trm(im,nc,nc),trms(im,nc,nc)
      common /tdf/ dx,dy,dz,trm,trms

      integer i,j,m,n

      do 100 i = 1,3*lm+3
         if (ttype(i) .eq. t1 .and. trank(i) .eq. n1+1) goto 110
100   continue
      return
110   do 120 j = 1,3*lm+3
         if (ttype(j) .eq. t2 .and. trank(j) .eq. n2+1) goto 130
120   continue
      return
130   continue

      do 200 m = 1,2*trank(i)+1
         do 210 n = 1,2*trank(j)+1
            if (tr) then
               g(tindex(i)+m,tindex(j)+n) =
     .               g(tindex(i)+m,tindex(j)+n) + f*gx(n,m)
            else
               g(tindex(i)+m,tindex(j)+n) =
     .               g(tindex(i)+m,tindex(j)+n) + f*gx(m,n)
            endif
210      continue
200   continue

      end


* Calculate complete symmetrized G matrix for one pair of particles

      subroutine gcalcs(xc,yc,zc)
      implicit none
      real xc,yc,zc

      real cutoff
#ifdef PERIODIC
      parameter (cutoff = 6.)
#else
      parameter (cutoff = 1.e6)
#endif

      integer lm
      parameter (lm = _LM_)

      integer im
      parameter (im = lm+1)
      integer nc
      parameter (nc = 2*im+1)
      integer dx(2*im,4*im+1),dy(2*im,4*im+1),dz(2*im,4*im+1)
      real trm(im,nc,nc),trms(im,nc,nc)
      common /tdf/ dx,dy,dz,trm,trms

      integer ttype(3*lm+3),trank(3*lm+3),tindex(3*lm+3)
      common /tens/ ttype,trank,tindex

      integer dim
      parameter (dim = 3*(lm+1)*(lm+3))
      real g(dim,dim),gt(dim,dim)
      common /mat/ g,gt

#ifdef PERIODIC
      common /pbc/ box
      real box
#endif

      real z(0:lm,4),zs(0:lm,4)
      common /zm/ z,zs

      integer i,j,k,m
#ifdef PERIODIC
      integer ns,nx,ny,nz
#endif
      logical ho,first
      real xl,yl,zl,s,st,r
      real gx(4*lm+5,4*lm+5)
      real temp(dim,2*lm+3)

      do 100 k = 1,dim
         do 110 m = 1,dim
            g(k,m) = 0.
110      continue
100   continue

      first = .true.

#ifdef PERIODIC
      call gpbc(xc,yc,zc)
      ns = int((2*cutoff+box-0.1)/box)-1
      do 150 nx = -ns,ns
      xl = xc + nx*box
      do 151 ny = -ns,ns
      yl = yc + ny*box
      do 152 nz = -ns,ns
      zl = zc + nz*box
#else
      xl = xc
      yl = yc
      zl = zc
#endif
      r = xl*xl+yl*yl+zl*zl
      ho = (r .le. cutoff*cutoff) .and. (r .gt. 0.)
      if (ho) call powers(xl,yl,zl)

      do 200 m = 0,lm
         st = 1.-2.*mod(m,2)
         do 210 k = m,lm
            s = 1.-2.*mod(k,2)

            if (ho .or. (k .eq. 0 .and. first)) then

               call gsp(gx,4*lm+5,m,k,first)

               call ins(g,dim,gx,4*lm+5,m,k,0,0,
     .                  s*(zs(m,1)*zs(k,4)+zs(m,4)*zs(k,1)),.false.)
               call ins(g,dim,gx,4*lm+5,m,k,0,2,
     .                  s*(zs(m,1)*zs(k,3)+zs(m,4)*zs(k,4)),.false.)
               call ins(g,dim,gx,4*lm+5,m,k,2,0,
     .                  s*(zs(m,3)*zs(k,1)+zs(m,4)*zs(k,4)),.false.)
               call ins(g,dim,gx,4*lm+5,m,k,2,2,
     .                  s*(zs(m,3)*zs(k,4)+zs(m,4)*zs(k,3)),.false.)
               call ins(g,dim,gx,4*lm+5,m,k,1,1,
     .                  -s*zs(m,2)*zs(k,2),.false.)
               if (m .ne. k) then
                  call ins(g,dim,gx,4*lm+5,k,m,0,0,
     .                  st*(zs(k,1)*zs(m,4)+zs(k,4)*zs(m,1)),.true.)
                  call ins(g,dim,gx,4*lm+5,k,m,0,2,
     .                  st*(zs(k,1)*zs(m,3)+zs(k,4)*zs(m,4)),.true.)
                  call ins(g,dim,gx,4*lm+5,k,m,2,0,
     .                  st*(zs(k,3)*zs(m,1)+zs(k,4)*zs(m,4)),.true.)
                  call ins(g,dim,gx,4*lm+5,k,m,2,2,
     .                  st*(zs(k,3)*zs(m,4)+zs(k,4)*zs(m,3)),.true.)
                  call ins(g,dim,gx,4*lm+5,k,m,1,1,
     .                  -st*zs(k,2)*zs(m,2),.true.)
               endif

            endif
            if (ho .or. (m+k .le. 2 .and. first)) then

               call gss(gx,4*lm+5,m,k,first)

               call ins(g,dim,gx,4*lm+5,m,k,0,0,
     .                              s*zs(m,1)*zs(k,1),.false.)
               call ins(g,dim,gx,4*lm+5,m,k,0,2,
     .                              s*zs(m,1)*zs(k,4),.false.)
               call ins(g,dim,gx,4*lm+5,m,k,2,0,
     .                              s*zs(m,4)*zs(k,1),.false.)
               call ins(g,dim,gx,4*lm+5,m,k,2,2,
     .                              s*zs(m,4)*zs(k,4),.false.)
               if (m .ne. k) then
                  call ins(g,dim,gx,4*lm+5,k,m,0,0,
     .                              st*zs(k,1)*zs(m,1),.true.)
                  call ins(g,dim,gx,4*lm+5,k,m,0,2,
     .                              st*zs(k,1)*zs(m,4),.true.)
                  call ins(g,dim,gx,4*lm+5,k,m,2,0,
     .                              st*zs(k,4)*zs(m,1),.true.)
                  call ins(g,dim,gx,4*lm+5,k,m,2,2,
     .                              st*zs(k,4)*zs(m,4),.true.)
               endif

            endif
            if (ho .or. (m+k .le. 1 .and. first)) then

               call gst(gx,4*lm+5,m,k,first)

               call ins(g,dim,gx,4*lm+5,m,k,0,1,
     .                              s*zs(m,1)*zs(k,2),.false.)
               call ins(g,dim,gx,4*lm+5,m,k,1,0,
     .                              s*zs(m,2)*zs(k,1),.false.)
               call ins(g,dim,gx,4*lm+5,m,k,1,2,
     .                              s*zs(m,2)*zs(k,4),.false.)
               call ins(g,dim,gx,4*lm+5,m,k,2,1,
     .                              s*zs(m,4)*zs(k,2),.false.)
               if (m .ne. k) then
                  call ins(g,dim,gx,4*lm+5,k,m,0,1,
     .                              -st*zs(k,1)*zs(m,2),.true.)
                  call ins(g,dim,gx,4*lm+5,k,m,1,0,
     .                              -st*zs(k,2)*zs(m,1),.true.)
                  call ins(g,dim,gx,4*lm+5,k,m,1,2,
     .                              -st*zs(k,2)*zs(m,4),.true.)
                  call ins(g,dim,gx,4*lm+5,k,m,2,1,
     .                              -st*zs(k,4)*zs(m,2),.true.)
               endif

            endif

210      continue
200   continue

      first = .false.

#ifdef PERIODIC
152   continue
151   continue
150   continue
#endif

      do 300 k = 1,3*lm+3
         do 310 i = 1,dim
            do 311 j = 1,2*trank(k)+1
               temp(i,j) = 0.
               do 312 m = 1,2*trank(k)+1
                  temp(i,j) = temp(i,j)
     .                  + g(i,tindex(k)+m)*trms(trank(k),m,j)
312            continue
311         continue
310      continue
         do 320 i = 1,dim
            do 321 j = 1,2*trank(k)+1
              g(i,tindex(k)+j) = temp(i,j)
321         continue
320      continue
300   continue

      do 400 k = 1,3*lm+3
         do 410 i = 1,dim
            do 411 j = 1,2*trank(k)+1
               temp(i,j) = 0.
               do 412 m = 1,2*trank(k)+1
                  temp(i,j) = temp(i,j)
     .                  + g(tindex(k)+m,i)*trms(trank(k),m,j)
412            continue
411         continue
410      continue
         do 420 i = 1,dim
            do 421 j = 1,2*trank(k)+1
              g(tindex(k)+j,i) = temp(i,j)
421         continue
420      continue
400   continue

      end


* Calculate complete G matrix for one pair of particles
      
      subroutine gcalc(xc,yc,zc)
      implicit none
      real xc,yc,zc

      real cutoff
#ifdef PERIODIC
      parameter (cutoff = 6.)
#else
      parameter (cutoff = 1.e6)
#endif
      
      integer lm
      parameter (lm = _LM_)
      
      integer im
      parameter (im = lm+1)
      integer nc
      parameter (nc = 2*im+1)
      integer dx(2*im,4*im+1),dy(2*im,4*im+1),dz(2*im,4*im+1)
      real trm(im,nc,nc),trms(im,nc,nc)
      common /tdf/ dx,dy,dz,trm,trms
      
      integer ttype(3*lm+3),trank(3*lm+3),tindex(3*lm+3)
      common /tens/ ttype,trank,tindex
      
      integer dim
      parameter (dim = 3*(lm+1)*(lm+3))
      real g(dim,dim),gt(dim,dim)
      common /mat/ g,gt
      
#ifdef PERIODIC
      common /pbc/ box
      real box
#endif

      real z(0:lm,4),zs(0:lm,4)
      common /zm/ z,zs
      
      integer i,j,k,m
#ifdef PERIODIC
      integer ns,nx,ny,nz
#endif
      logical ho,first
      real xl,yl,zl,s,st,r
      real gx(4*lm+5,4*lm+5)
      real temp(dim,2*lm+3)
      
      do 100 k = 1,dim
         do 110 m = 1,dim
            g(k,m) = 0.
            gt(k,m) = 0.
110      continue
100   continue
      
      first = .true.

#ifdef PERIODIC
      call gpbc(xc,yc,zc)
      ns = int((2*cutoff+box-0.1)/box)-1
      do 150 nx = -ns,ns
      xl = xc + nx*box
      do 151 ny = -ns,ns
      yl = yc + ny*box
      do 152 nz = -ns,ns
      zl = zc + nz*box
#else
      xl = xc
      yl = yc
      zl = zc
#endif
      r = xl*xl+yl*yl+zl*zl
      ho = (r .le. cutoff*cutoff) .and. (r .gt. 0.)
      if (ho) call powers(xl,yl,zl)
      
      do 200 m = 0,lm
         st = 1.-2.*mod(m,2)
         do 210 k = m,lm
            s = 1.-2.*mod(k,2)
            
            if (ho .or. (k .eq. 0 .and. first)) then

               call gsp(gx,4*lm+5,m,k,first)
            
               call ins(g,dim,gx,4*lm+5,m,k,0,2,s*z(m,1),.false.)
               call ins(g,dim,gx,4*lm+5,m,k,2,2,s*z(m,4),.false.)
               call ins(g,dim,gx,4*lm+5,m,k,2,0,s*z(m,3),.false.)
               call ins(g,dim,gx,4*lm+5,m,k,0,0,s*z(m,4),.false.)
               call ins(g,dim,gx,4*lm+5,m,k,1,1,-s*z(m,2),.false.)
               if (m .ne. k) then
                  call ins(g,dim,gx,4*lm+5,k,m,0,2,st*z(k,1),.true.)
                  call ins(g,dim,gx,4*lm+5,k,m,2,2,st*z(k,4),.true.)
                  call ins(g,dim,gx,4*lm+5,k,m,2,0,st*z(k,3),.true.)
                  call ins(g,dim,gx,4*lm+5,k,m,0,0,st*z(k,4),.true.)
                  call ins(g,dim,gx,4*lm+5,k,m,1,1,-st*z(k,2),.true.)
               endif

               call ins(gt,dim,gx,4*lm+5,k,m,0,2,s*z(k,1),.true.)
               call ins(gt,dim,gx,4*lm+5,k,m,2,2,s*z(k,4),.true.)
               call ins(gt,dim,gx,4*lm+5,k,m,2,0,s*z(k,3),.true.)
               call ins(gt,dim,gx,4*lm+5,k,m,0,0,s*z(k,4),.true.)
               call ins(gt,dim,gx,4*lm+5,k,m,1,1,-s*z(k,2),.true.)
               if (m .ne. k) then
                  call ins(gt,dim,gx,4*lm+5,m,k,0,2,st*z(m,1),.false.)
                  call ins(gt,dim,gx,4*lm+5,m,k,2,2,st*z(m,4),.false.)
                  call ins(gt,dim,gx,4*lm+5,m,k,2,0,st*z(m,3),.false.)
                  call ins(gt,dim,gx,4*lm+5,m,k,0,0,st*z(m,4),.false.)
                  call ins(gt,dim,gx,4*lm+5,m,k,1,1,-st*z(m,2),.false.)
               endif

            endif
            if (ho .or. (m+k .le. 2 .and. first)) then
            
               call gss(gx,4*lm+5,m,k,first)
            
               call ins(g,dim,gx,4*lm+5,m,k,0,0,s*z(m,1),.false.)
               call ins(g,dim,gx,4*lm+5,m,k,2,0,s*z(m,4),.false.)
               if (m .ne. k) then
                  call ins(g,dim,gx,4*lm+5,k,m,0,0,st*z(k,1),.true.)
                  call ins(g,dim,gx,4*lm+5,k,m,2,0,st*z(k,4),.true.)
               endif
            
               call ins(gt,dim,gx,4*lm+5,k,m,0,0,s*z(k,1),.true.)
               call ins(gt,dim,gx,4*lm+5,k,m,2,0,s*z(k,4),.true.)
               if (m .ne. k) then
                  call ins(gt,dim,gx,4*lm+5,m,k,0,0,st*z(m,1),.false.)
                  call ins(gt,dim,gx,4*lm+5,m,k,2,0,st*z(m,4),.false.)
               endif
            
            endif
            if (ho .or. (m+k .le. 1 .and. first)) then
            
               call gst(gx,4*lm+5,m,k,first)
            
               call ins(g,dim,gx,4*lm+5,m,k,0,1,s*z(m,1),.false.)
               call ins(g,dim,gx,4*lm+5,m,k,2,1,s*z(m,4),.false.)
               call ins(g,dim,gx,4*lm+5,m,k,1,0,s*z(m,2),.false.)
               if (m .ne. k) then
                  call ins(g,dim,gx,4*lm+5,k,m,0,1,-st*z(k,1),.true.)
                  call ins(g,dim,gx,4*lm+5,k,m,2,1,-st*z(k,4),.true.)
                  call ins(g,dim,gx,4*lm+5,k,m,1,0,-st*z(k,2),.true.)
               endif
            
               call ins(gt,dim,gx,4*lm+5,k,m,0,1,s*z(k,1),.true.)
               call ins(gt,dim,gx,4*lm+5,k,m,2,1,s*z(k,4),.true.)
               call ins(gt,dim,gx,4*lm+5,k,m,1,0,s*z(k,2),.true.)
               if (m .ne. k) then
                  call ins(gt,dim,gx,4*lm+5,m,k,0,1,-st*z(m,1),.false.)
                  call ins(gt,dim,gx,4*lm+5,m,k,2,1,-st*z(m,4),.false.)
                  call ins(gt,dim,gx,4*lm+5,m,k,1,0,-st*z(m,2),.false.)
               endif
            
            endif

210      continue
200   continue
      
      first = .false.

#ifdef PERIODIC
152   continue
151   continue
150   continue
#endif

      do 300 k = 1,3*lm+3
         do 310 i = 1,dim
            do 311 j = 1,2*trank(k)+1
               temp(i,j) = 0.
               do 312 m = 1,2*trank(k)+1
                  temp(i,j) = temp(i,j)
     .                 + g(i,tindex(k)+m)*trm(trank(k),m,j)
312            continue
311         continue
310      continue
         do 320 i = 1,dim
            do 321 j = 1,2*trank(k)+1
               g(i,tindex(k)+j) = temp(i,j)
321         continue
320      continue
         do 330 i = 1,dim
            do 331 j = 1,2*trank(k)+1
               temp(i,j) = 0.
               do 332 m = 1,2*trank(k)+1
                  temp(i,j) = temp(i,j)
     .                 + gt(i,tindex(k)+m)*trm(trank(k),m,j)
332            continue
331         continue
330      continue
         do 340 i = 1,dim
            do 341 j = 1,2*trank(k)+1
               gt(i,tindex(k)+j) = temp(i,j)
341         continue
340      continue
300   continue
      
      end
#endif
