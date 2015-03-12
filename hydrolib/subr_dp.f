*************************************
* Friction and mobility             *
* for a hard-sphere configuration   *
* Module SUBR                       *
*                                   *
* K. Hinsen                         *
* Last revision: September 29, 1994 *
*************************************

#ifndef _LMS_
#define _LMS_ _LM_
#endif

#ifdef LUBRICATION
#include <filenames.h>
#endif

#ifdef INCLUDE
#include <gcalc.f>
#ifdef LUBRICATION
#include <lubrication.f>
#endif
#endif

#ifdef RB_VELOCITIES
#ifndef VELOCITIES
#define VELOCITIES
#endif
#endif

#ifdef RIGID

#ifdef MOBILITY
#ifndef FRICTION
#define FRICTION
#define MOBILITY_DIRECT
#endif
#undef MOBILITY
#define MOBILITY_RB
#endif

#ifdef VELOCITIES
#ifndef FRICTION
#define FRICTION
#define MOBILITY_DIRECT
#endif
#undef VELOCITIES
#define VELOCITIES_RB
#endif

#ifdef NORMALIZE
#undef NORMALIZE
#define NORMALIZE_RB
#endif

#include <rigid.f>

#endif

#ifdef FIXED_FAST
#define FIXED
#endif

#ifdef FIXED

#ifdef MOBILITY
#ifndef FRICTION
#define FRICTION
#define MOBILITY_DIRECT
#endif
#endif

#ifdef VELOCITIES
#ifndef FRICTION
#define FRICTION
#define MOBILITY_DIRECT
#endif
#endif

#ifdef NORMALIZE
#undef NORMALIZE
#define NORMALIZE_RB
#endif

#endif

#ifndef MOBILITY_DIRECT
#ifdef MOBILITY
#define MOBILITY_DIRECT
#endif
#ifdef VELOCITIES
#define MOBILITY_DIRECT
#endif
#ifdef FRICTION
#undef MOBILITY_DIRECT
#endif
#endif

******************
* initialization *      
******************

      subroutine init
      implicit none

      double precision xi
      parameter (xi = 0.d0)
#ifdef LUBRICATION
      character*150 lname
#endif
      integer lm
      parameter (lm = _LM_)
      integer lms
      parameter (lms = _LMS_)

#if _LM_ < _LMS_
      stop 'Inconsistent multipole order specifications!'
#endif
      call comp(xi)
      call ocload
#ifdef LUBRICATION
      lname = 'IMPOSSIBLE'
      if (lm .eq. lms) then
         if (lm .eq. 0)
     .         lname =
     .FILENAME3
         if (lm .eq. 1)
     .         lname =
     .FILENAME4
         if (lm .eq. 2)
     .         lname =
     .FILENAME5
         if (lm .eq. 3)
     .         lname =
     .FILENAME6
      else
         if (lm .eq. 2 .and. lms .eq. 1)
     .         lname =
     .FILENAME7
      endif
      call luload(lname)
#endif
#ifdef PERIODIC
      call cload
      call pbccoef
#endif

      end



* Calculate the full symmetric interaction matrix

#ifdef FRICTION
      subroutine initts(c,np)
      implicit none
      integer np
      double precision c(0:2,np)

      integer npart
      parameter (npart = _NP_)

      integer lm
      parameter (lm = _LM_)
      integer lms
      parameter (lms = _LMS_)

#ifdef PERIODIC      
      common /pbc/ box
      double precision box

      double precision g0s0s(3,3),g0s1s(3,5),g0s2s(3,7),g1s1s(5,5)
      double precision g0s0t(3,3),g0s1t(3,5),g0s0p(3,3),s1v
      common /gpbca/ g0s0s,g0s1s,g0s2s,g1s1s,g0s0t,g0s1t,g0s0p,s1v
#endif

      integer dim
      parameter (dim = 3*(lm+1)*(lm+3))
      double precision g(dim,dim),gt(dim,dim)
      common /mat/ g,gt

      integer dims
      parameter (dims = 3*(lms+1)*(lms+3))

#ifdef PACKED
      double precision t(npart*dims*(npart*dims+1)/2)
#else
      double precision t(npart*dims,npart*dims)
#endif
      common t

#ifdef FIXED_FAST
#ifdef PACKED
      double precision tb(npart*dims*(npart*dims+1)/2)
#else
      double precision tb(npart*dims,npart*dims)
#endif
      save tb

      integer nfixed
      logical full
      common /fix/ nfixed,full
#endif

#if _LM_ != _LMS_
      integer ddim
      parameter (ddim = dim-dims)
      double precision b(npart*dims,ddim)
#endif

      integer i,j,k,m
      double precision rx,ry,rz

      if (np .ne. npart) stop 'Incompatible parameters!'

#ifdef PACKED
      do 50 i = 1,npart*dims*(npart*dims+1)/2
         t(i) = 0.d0
50    continue
#else
      do 50 i = 1,npart*dims
         do 51 j = 1,npart*dims
            t(i,j) = 0.d0
51       continue
50    continue
#endif

#if _LM_ == _LMS_
      do 100 j = 0,np-2
         do 110 i = j+1,np-1
#ifdef FIXED_FAST
           if (.not. full .and. j+1 .gt. np-nfixed) then
            do 120 k = 1,dims
               do 121 m = 1,dims
#ifdef PACKED
                  t(j*dims+k+(i*dims+m-1)*(i*dims+m)/2) =
     .               tb(j*dims+k+(i*dims+m-1)*(i*dims+m)/2)
#else
                  t(i*dims+k,j*dims+m) = tb(i*dims+k,j*dims+m)
                  t(j*dims+k,i*dims+m) = tb(j*dims+k,i*dims+m)
#endif
121            continue
120         continue
           else
#endif
            rx = c(0,i+1)-c(0,j+1)
            ry = c(1,i+1)-c(1,j+1)
            rz = c(2,i+1)-c(2,j+1)
#ifdef PERIODIC
            if (rx .gt. 0.5d0*box) rx = rx-box
            if (rx .lt. -0.5d0*box) rx = rx+box
            if (ry .gt. 0.5d0*box) ry = ry-box
            if (ry .lt. -0.5d0*box) ry = ry+box
            if (rz .gt. 0.5d0*box) rz = rz-box
            if (rz .lt. -0.5d0*box) rz = rz+box
#endif
            call gcalcs(rx,ry,rz)
            do 130 k = 1,dims
               do 131 m = 1,dims
#ifdef PACKED
                  t(j*dims+k+(i*dims+m-1)*(i*dims+m)/2) = g(m,k)
#else
                  t(i*dims+k,j*dims+m) = g(k,m)
                  t(j*dims+k,i*dims+m) = g(m,k)
#endif
131            continue
130         continue
#ifdef FIXED_FAST
           endif
#endif
110      continue
100   continue
#else
      do 100 j = 0,np-1
         do 110 i = 0,np-1
            if (i .ne. j) then
               rx = c(0,i+1)-c(0,j+1)
               ry = c(1,i+1)-c(1,j+1)
               rz = c(2,i+1)-c(2,j+1)
#ifdef PERIODIC
               if (rx .gt. 0.5d0*box) rx = rx-box
               if (rx .lt. -0.5d0*box) rx = rx+box
               if (ry .gt. 0.5d0*box) ry = ry-box
               if (ry .lt. -0.5d0*box) ry = ry+box
               if (rz .gt. 0.5d0*box) rz = rz-box
               if (rz .lt. -0.5d0*box) rz = rz+box
#endif
               call gcalcs(rx,ry,rz)
               do 130 k = 1,dims
                  do 131 m = 1,ddim
                     b(i*dims+k,m) = g(k,dims+m)
131               continue
130            continue
            endif
            if (i .gt. j) then
               do 120 k = 1,dims
                  do 121 m = 1,dims
#ifdef PACKED
                     t(j*dims+k+(i*dims+m-1)*(i*dims+m)/2) =
     .                    t(j*dims+k+(i*dims+m-1)*(i*dims+m)/2) + g(m,k)
#else
                     t(i*dims+k,j*dims+m) =
     .                    t(i*dims+k,j*dims+m) + g(k,m)
                     t(j*dims+k,i*dims+m) =
     .                    t(j*dims+k,i*dims+m) + g(m,k)
#endif
121               continue
120            continue
            endif
110      continue
         do 140 k = 1,dims
            do 141 m = 1,ddim
               b(j*dims+k,m) = 0
141         continue
140      continue
#ifdef PACKED
         do 150 i = 1,dims*np
            do 151 k = i,dims*np
               do 152 m = 1,ddim
                  t(i+(k-1)*k/2) = t(i+(k-1)*k/2) - b(i,m)*b(k,m)
152            continue
151         continue
150      continue
#else
         do 150 i = 1,dims*np
            do 151 k = 1,dims*np
               do 152 m = 1,ddim
                  t(i,k) = t(i,k) - b(i,m)*b(k,m)
152            continue
151         continue
150      continue
#endif
100   continue
#endif

#ifdef FIXED_FAST
      if (full) then
#ifdef PACKED
         do 170 i = 1,npart*dims*(npart*dims+1)/2
            tb(i) = t(i)
170      continue
#else
         do 170 i = 1,npart*dims
            do 171 j = 1,npart*dims
               tb(i,j) = t(i,j)
171         continue
170      continue
#endif
      endif
#endif

#ifdef PERIODIC
      call gcalcs(0.d0,0.d0,0.d0)
      do 200 i = 0,np-1
#ifdef PACKED
         do 210 k = 1,dims
            do 220 m = k,dims
               t(i*dims+k+(i*dims+m-1)*(i*dims+m)/2) =
     .              t(i*dims+k+(i*dims+m-1)*(i*dims+m)/2) + g(k,m)
220         continue
            t((i*dims+k+1)*(i*dims+k)/2) =
     .           1.d0 + t((i*dims+k+1)*(i*dims+k)/2)
210      continue
#else
         do 210 k = 1,dims
            do 220 m = 1,dims
               t(i*dims+k,i*dims+m) = t(i*dims+k,i*dims+m) + g(k,m)
220         continue
            t(i*dims+k,i*dims+k) = 1.d0 + t(i*dims+k,i*dims+k)
210      continue
#endif
200   continue
#else
      do 200 i = 0,np-1
#ifdef PACKED
         do 210 k = 1,dims
            t((i*dims+k+1)*(i*dims+k)/2) =
     .           1.d0 + t((i*dims+k+1)*(i*dims+k)/2)
210      continue
#else
         do 210 k = 1,dims
            t(i*dims+k,i*dims+k) = 1.d0 + t(i*dims+k,i*dims+k)
210      continue
#endif
200   continue
#endif

      end
#endif
      
      
* Intialize velocities to calculate friction matrix

#ifdef FRICTION      
      subroutine vfmat(i,np)
      implicit none
      integer i,np

      integer npart
      parameter (npart = _NP_)
      
      integer lm
      parameter (lm = _LM_)
      integer lms
      parameter (lms = _LMS_)
      
      integer dim
      parameter (dim = 3*(lm+1)*(lm+3))
      integer dims
      parameter (dims = 3*(lms+1)*(lms+3))
      double precision v(dims*npart,6)
      common /vec/ v
      
      integer ttype(3*lm+3),trank(3*lm+3),tindex(3*lm+3)
      common /tens/ ttype,trank,tindex

      double precision z(0:lm,4),zs(0:lm,4)
      common /zm/ z,zs
      
      integer j
      
      if (np .ne. npart) stop 'Incompatible parameters!'

      do 100 j = 1,dims*npart
         v(j,1) = 0.d0
         v(j,2) = 0.d0
         v(j,3) = 0.d0
         v(j,4) = 0.d0
         v(j,5) = 0.d0
         v(j,6) = 0.d0
100   continue

      v((i-1)*dims+3,1) = -zs(0,1)
      v((i-1)*dims+2,2) = -zs(0,1)
      v((i-1)*dims+1,3) = -zs(0,1)
      v((i-1)*dims+9,1) = -zs(0,4)
      v((i-1)*dims+8,2) = -zs(0,4)
      v((i-1)*dims+7,3) = -zs(0,4)
      v((i-1)*dims+4,6) = -2.d0*zs(0,2)
      v((i-1)*dims+5,5) = -2.d0*zs(0,2)
      v((i-1)*dims+6,4) = -2.d0*zs(0,2)

      end
#endif
      
* Calculate friction and mobility matrices for a single configuration
      
#ifdef FRICTION
      subroutine evalf
      implicit none

      integer npart
      parameter (npart = _NP_)
      
      integer lm
      parameter (lm = _LM_)
      integer lms
      parameter (lms = _LMS_)
      
      integer dim
      parameter (dim = 3*(lm+1)*(lm+3))
      integer dims
      parameter (dims = 3*(lms+1)*(lms+3))
      double precision v(dims*npart,6)
      common /vec/ v
      
#ifdef PACKED
      double precision t(npart*dims*(npart*dims+1)/2)
#else
      double precision t(npart*dims,npart*dims)
#endif
      common t
      
      double precision c(0:2,npart)
#ifdef RIGID
      character*1 cnct(npart,npart)
      common /conf/ c,cnct
#else
      common /conf/ c
#endif
#ifdef FIXED
      integer nfixed
      logical full
      common /fix/ nfixed,full
#else
      integer nfixed
      parameter (nfixed = 0)
#endif

      double precision fr(6*npart,6*npart)
      common /frict/ fr

      double precision z(0:lm,4),zs(0:lm,4)
      common /zm/ z,zs
      
#ifndef POSDEF
      integer ipvt(dims*npart)
#ifdef LAPACK
#ifndef PACKED
      integer lwork
      parameter(lwork = 10*dims*npart)
      double precision work(lwork)
#endif
#endif
#endif

      integer i,j,k

#ifdef IMSL
#ifdef PACKED
      stop 'IMSL does not support packed arrays!'
#endif
#endif

      call initts(c,npart)
#ifdef POSDEF
#ifdef IMSL
#ifdef DP
      call dlftds(dims*npart,t,dims*npart,t,dims*npart)
#else
      call lftds(dims*npart,t,dims*npart,t,dims*npart)
#endif
#endif
#ifdef LAPACK
#ifdef PACKED
#ifdef DP
      call dpptrf('U',dims*npart,t,i)
#else
      call spptrf('U',dims*npart,t,i)
#endif
#else
#ifdef DP
      call dpotrf('U',dims*npart,t,dims*npart,i)
#else
      call spotrf('U',dims*npart,t,dims*npart,i)
#endif
#endif
      if (i .ne. 0) print *,'Error during factorization 1! (',i,')'
#endif
#else /* not positive definite */
#ifdef IMSL
#ifdef DP
      call dlftsf(dims*npart,t,dims*npart,t,dims*npart,ipvt)
#else
      call lftsf(dims*npart,t,dims*npart,t,dims*npart,ipvt)
#endif
#endif
#ifdef LAPACK
#ifdef PACKED
#ifdef DP
      call dsptrf('U',dims*npart,t,ipvt,i)
#else
      call ssptrf('U',dims*npart,t,ipvt,i)
#endif
#else
#ifdef DP
      call dsytrf('U',dims*npart,t,dims*npart,ipvt,work,lwork,i)
#else
      call ssytrf('U',dims*npart,t,dims*npart,ipvt,work,lwork,i)
#endif
#endif
      if (i .ne. 0) print *,'Error during factorization 1! (',i,')'
#endif
#endif

      do 100 i = 0,npart-nfixed-1
         call vfmat(i+1,npart)
#ifdef IMSL
         do 101 j = 1,6
#ifdef POSDEF
#ifdef DP
            call dlfsds(dims*npart,t,dims*npart,v(1,j),v(1,j))
#else
            call lfsds(dims*npart,t,dims*npart,v(1,j),v(1,j))
#endif
#else /* not positive definite */
#ifdef DP
            call dlfssf(dims*npart,t,dims*npart,ipvt,v(1,j),v(1,j))
#else
            call lfssf(dims*npart,t,dims*npart,ipvt,v(1,j),v(1,j))
#endif
#endif
101      continue
#endif
#ifdef LAPACK
#ifdef POSDEF
#ifdef PACKED
#ifdef DP
         call dpptrs('U',dims*npart,6,t,v,dims*npart,j)
#else
         call spptrs('U',dims*npart,6,t,v,dims*npart,j)
#endif
#else
#ifdef DP
         call dpotrs('U',dims*npart,6,t,dims*npart,v,dims*npart,j)
#else
         call spotrs('U',dims*npart,6,t,dims*npart,v,dims*npart,j)
#endif
#endif
#else /* not positive definite */
#ifdef PACKED
#ifdef DP
         call dsptrs('U',dims*npart,6,t,ipvt,v,dims*npart,j)
#else
         call ssptrs('U',dims*npart,6,t,ipvt,v,dims*npart,j)
#endif
#else
#ifdef DP
         call dsytrs('U',dims*npart,6,t,dims*npart,ipvt,v,dims*npart,j)
#else
         call ssytrs('U',dims*npart,6,t,dims*npart,ipvt,v,dims*npart,j)
#endif
#endif
#endif
         if (j .ne. 0) print *,'Error during solution 1! (',j,')'
#endif

         do 120 j = 0,npart-nfixed-1
            do 130 k = 1,6
               fr(6*j+1,6*i+k) =
     .              -zs(0,1)*v(j*dims+3,k)-zs(0,4)*v(j*dims+9,k)
               fr(6*j+2,6*i+k) =
     .              -zs(0,1)*v(j*dims+2,k)-zs(0,4)*v(j*dims+8,k)
               fr(6*j+3,6*i+k) =
     .              -zs(0,1)*v(j*dims+1,k)-zs(0,4)*v(j*dims+7,k)
               fr(6*j+4,6*i+k) = -2.d0*zs(0,2)*v(j*dims+6,k)
               fr(6*j+5,6*i+k) = -2.d0*zs(0,2)*v(j*dims+5,k)
               fr(6*j+6,6*i+k) = -2.d0*zs(0,2)*v(j*dims+4,k)
130         continue
120      continue
100   continue

#ifdef LUBRICATION
#ifdef RIGID
      call lub(c,cnct,fr,npart,nfixed,.true.,lms)
#else
      call lub(c,c,fr,npart,nfixed,.false.,lms)
#endif
#endif

#ifdef VELOCITIES
      call solvev
#endif

#ifdef MOBILITY
      call inv(npart-nfixed)
#endif

#ifdef NORMALIZE
      call norm(npart-nfixed)
#endif

      end
#endif


* Call appropriate evaluation routine

      subroutine eval
      implicit none

#ifdef FIXED
      integer nfixed
      logical full
      common /fix/ nfixed,full
#else
      integer nfixed
      parameter (nfixed = 0)
#endif

#ifdef RIGID
      call evalfr
#else
#ifdef MOBILITY_DIRECT
#ifdef FIXED
      if (nfixed .eq. 0) then
         call evalm
      else
         call evalf
      endif
#else
      call evalm
#endif
#else
      call evalf
#endif
#endif
      end


#ifdef MOBILITY_RB
#define MOBILITY
#endif

#ifdef VELOCITIES_RB
#define VELOCITIES
#endif

#ifdef NORMALIZE_RB
#define NORMALIZE
#endif


* Calculate the full interaction matrix for mobility evaluations
 
#ifdef MOBILITY_DIRECT
      subroutine initzg(c,np)
      implicit none
      integer np
      double precision c(0:2,np)

      integer npart
      parameter (npart = _NP_)
 
      integer lm
      parameter (lm = _LM_)
      integer lms
      parameter (lms = _LMS_)

#if _LM_ != _LMS_
      stop 'Truncated scheme not yet implemented!'
#endif
#ifdef PERIODIC      
      common /pbc/ box
      double precision box

      double precision g0s0s(3,3),g0s1s(3,5),g0s2s(3,7),g1s1s(5,5)
      double precision g0s0t(3,3),g0s1t(3,5),g0s0p(3,3),s1v
      common /gpbca/ g0s0s,g0s1s,g0s2s,g1s1s,g0s0t,g0s1t,g0s0p,s1v
#endif

      integer dim
      parameter (dim = 3*(lm+1)*(lm+3))
      double precision g(dim,dim),gt(dim,dim)
      common /mat/ g,gt
 
      double precision z(0:lm,4),zs(0:lm,4)
      common /zm/ z,zs
 
      double precision zg(npart*dim,npart*dim)
      common zg
 
      integer i,j,k,m
      double precision rx,ry,rz
 
      if (np .ne. npart) stop 'Incompatible parameters!'
 
      do 100 j = 0,np-2
         do 110 i = j+1,np-1
            rx = c(0,i+1)-c(0,j+1)
            ry = c(1,i+1)-c(1,j+1)
            rz = c(2,i+1)-c(2,j+1)
#ifdef PERIODIC
            if (rx .gt. 0.5d0*box) rx = rx-box
            if (rx .lt. -0.5d0*box) rx = rx+box
            if (ry .gt. 0.5d0*box) ry = ry-box
            if (ry .lt. -0.5d0*box) ry = ry+box
            if (rz .gt. 0.5d0*box) rz = rz-box
            if (rz .lt. -0.5d0*box) rz = rz+box
#endif
            call gcalc(rx,ry,rz)
            do 120 k = 1,dim
               do 130 m = 1,dim
                  zg(i*dim+k,j*dim+m) = -g(k,m)
                  zg(j*dim+k,i*dim+m) = -gt(k,m)
130            continue
120         continue
110      continue
100   continue
 
#ifdef PERIODIC
      call gcalc(0.d0,0.d0,0.d0)
      do 200 i = 0,np-1
         do 210 k = 1,dim
            do 220 m = 1,dim
               zg(i*dim+k,i*dim+m) = -g(k,m)
220         continue
210      continue
200   continue
#else
      do 200 i = 0,np-1
         do 210 k = 1,dim
            do 220 m = 1,dim
               zg(i*dim+k,i*dim+m) = 0.d0
220         continue
210      continue
200   continue
#endif
 
      end
#endif     
      
* Initialize "velocities"

#ifdef MOBILITY_DIRECT
      subroutine initmp
      implicit none

      integer npart
      parameter (npart = _NP_)

      integer lm
      parameter (lm = _LM_)
      integer lms
      parameter (lms = _LMS_)

      integer dim
      parameter (dim = 3*(lm+1)*(lm+3))
      integer dims
      parameter (dims = 3*(lms+1)*(lms+3))

#ifdef MOBILITY
      double precision mp(dim*npart,6*npart)
      common /vecm/ mp
#endif
#ifdef VELOCITIES
      double precision mpf(dim*npart)
      common /vecf/ mpf
#endif

#ifdef VELOCITIES
      double precision f(6*npart),v(6*npart)
      common /fv/ f,v
#endif

      double precision zg(npart*dim,npart*dim)
      common zg

      integer ttype(3*lm+3),trank(3*lm+3),tindex(3*lm+3)
      common /tens/ ttype,trank,tindex
 
      double precision z(0:lm,4),zs(0:lm,4)
      common /zm/ z,zs
      
#ifdef VELOCITIES
      integer ti(6)
      data ti /3,2,1,6,5,4/
#endif

      integer i,j,k
      
#if _LM_ != _LMS_
      stop 'Truncated scheme not yet implemented!'
#endif
      do 100 i = 1,dim*npart
#ifdef MOBILITY
         do 101 j = 0,npart-1
            do 102 k = 1,6
               mp(i,6*j+k) = zg(i,dim*j+k)
102         continue
101      continue
#endif
#ifdef VELOCITIES
         mpf(i) = 0.d0
         do 105 j = 0,npart-1
            do 106 k = 1,3
               mpf(i) = mpf(i) + zg(i,dim*j+k)*f(6*j+ti(k))
106         continue
            do 107 k = 4,6
               mpf(i) = mpf(i) + zg(i,dim*j+k)*f(6*j+ti(k))
     .              /(4.d0*dabs(z(0,2)))
107         continue
105      continue
#endif
100   continue

      end
#endif

* Calculate short-range contribution for mobility

#ifdef MOBILITY_DIRECT
#ifdef LUBRICATION
      subroutine lubmo
      implicit none

      integer npart
      parameter (npart = _NP_)
      
      integer lm
      parameter (lm = _LM_)
      integer lms
      parameter (lms = _LMS_)

      integer dim
      parameter (dim = 3*(lm+1)*(lm+3))
      integer dims
      parameter (dims = 3*(lms+1)*(lms+3))

      double precision zg(npart*dim,npart*dim)
      common zg
      
      double precision dz(dim*npart,6*npart)
      common /lz/ dz
      
      double precision c(0:2,npart)
#ifdef RIGID
      character*1 cnct(npart,npart)
      common /conf/ c,cnct
#else
      common /conf/ c
#endif

#ifdef FIXED
      integer nfixed
      logical full
      common /fix/ nfixed,full
#endif

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
      double precision frij(6,6),fri(6,6),frj(6,6)
      integer i,j,k,l,m

      double precision z1(6),fc(6)
      logical ndivp
      integer ti(6)
      data z1 /3*1.5d0,3*1.d0/
      data fc /3*1.d0,3*2.d0/
      data ti /3,2,1,6,5,4/

      double precision cutoff(0:3)
      data cutoff /4.d0,3.8d0,3.d0,2.6d0/
      
#if _LM_ != _LMS_
      stop 'Truncated scheme not yet implemented!'
#endif

      do 50 i = 1,dim*npart
         do 51 j = 1,6*npart
            dz(i,j) = 0.d0
51       continue
50    continue

      dr = rv(1)-rv(0)
      do 100 i = 1,npart
         do 110 j = i+1,npart
            r(1) = c(0,j)-c(0,i)
            r(2) = c(1,j)-c(1,i)
            r(3) = c(2,j)-c(2,i)
#ifdef PERIODIC
            if (r(1) .gt. 0.5d0*box) r(1) = r(1)-box
            if (r(1) .lt. -0.5d0*box) r(1) = r(1)+box
            if (r(2) .gt. 0.5d0*box) r(2) = r(2)-box
            if (r(2) .lt. -0.5d0*box) r(2) = r(2)+box
            if (r(3) .gt. 0.5d0*box) r(3) = r(3)-box
            if (r(3) .lt. -0.5d0*box) r(3) = r(3)+box
#endif
#ifdef RIGID
            ndivp = cnct(i,j) .ne. ' '
#else
            ndivp = .false.
#endif
#ifdef FIXED
            ndivp = ndivp .or.
     .         ((i .gt. npart-nfixed) .and. (j .gt. npart-nfixed))
#endif
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
                     frij(l,m) =ftt(l,m,2)
                     frij(l+3,m+3) = frr(l,m,2)
                     frij(l+3,m) = frt(l,m,2)
                     frij(l,m+3) = frt(l,m,2)

                     fri(l,m) = ftt(l,m,1)
                     fri(l+3,m+3) = frr(l,m,1)
                     fri(l+3,m) = frt(l,m,1)
                     fri(l,m+3) = frt(m,l,1)

                     frj(l,m) = ftt(m,l,1)
                     frj(l+3,m+3) = frr(m,l,1)
                     frj(l+3,m) = frt(m,l,1)
                     frj(l,m+3) = frt(l,m,1)
131               continue
130            continue
               do 250 k = 1,dim*npart
                  do 202 l = 1,6
                     do 211 m = 1,6
                        dz(k,6*(i-1)+ti(l)) = dz(k,6*(i-1)+ti(l))
     .                       + (zg(k,dim*(j-1)+ti(m))*frij(l,m)
     .                          + zg(k,dim*(i-1)+ti(m))*fri(m,l))
     .                       /(z1(l)*fc(m))
                        dz(k,6*(j-1)+ti(l)) = dz(k,6*(j-1)+ti(l))
     .                       + (zg(k,dim*(i-1)+ti(m))*frij(m,l)
     .                          + zg(k,dim*(j-1)+ti(m))*frj(m,l))
     .                       /(z1(l)*fc(m))
211                  continue
202               continue
250            continue
            endif
110      continue
100   continue
      
      end
#endif
#endif


* Calculate mobility matrix for a single configuration
      
#ifdef MOBILITY_DIRECT
      subroutine evalm      
      implicit none

      integer npart
      parameter (npart = _NP_)
      
      integer lm
      parameter (lm = _LM_)
      integer lms
      parameter (lms = _LMS_)
      
      integer dim
      parameter (dim = 3*(lm+1)*(lm+3))
      integer dims
      parameter (dims = 3*(lms+1)*(lms+3))

#ifdef MOBILITY
      double precision mp(dim*npart,6*npart)
      common /vecm/ mp
#endif
#ifdef VELOCITIES
      double precision mpf(dim*npart)
      common /vecf/ mpf
#endif
      
      double precision zg(npart*dim,npart*dim)
      common zg
      
#ifdef LUBRICATION
      double precision dz(dim*npart,6*npart)
      common /lz/ dz
#endif
      
      double precision c(0:2,npart)
#ifdef RIGID
      character*1 cnct(npart,npart)
      common /conf/ c,cnct
#else
      common /conf/ c
#endif

#ifdef VELOCITIES
      double precision f(6*npart),v(6*npart)
      common /fv/ f,v
#endif

#ifdef MOBILITY
      double precision mo(6*npart,6*npart)
      common /mobil/ mo
#endif

      double precision z(0:lm,4),zs(0:lm,4)
      common /zm/ z,zs
      
      integer i,j,k,l,m,indx(dim*npart)
      integer ti(6)
      double precision fa
      data ti /3,2,1,6,5,4/

#if _LM_ != _LMS_
      stop 'Truncated scheme not yet implemented!'
#endif

      call initzg(c,npart)
      do 100 i = 1,dim*npart
         zg(i,i) = 1.d0+zg(i,i)
100   continue

      call initmp

#ifdef LUBRICATION
      call lubmo
#endif

      fa = -z(0,4)/z(0,1)
      do 200 i = 1,dim*npart
         zg(i,i) = zg(i,i)-1.d0
200   continue
      do 210 i = 1,dim*npart
         do 211 j = 0,npart-1
            do 212 m = 1,3
               zg(i,j*dim+m) = fa*zg(i,j*dim+m+6)
               zg(i,j*dim+m+3) = 0.d0
212         continue
211      continue
210   continue
      do 220 i = 1,dim*npart
         zg(i,i) = zg(i,i)+1.d0
220   continue

#ifdef LUBRICATION
      do 230 i = 1,dim*npart
         do 231 j = 0,npart-1
            do 232 k = 1,6
               zg(i,dim*j+k) = zg(i,dim*j+k) + dz(i,6*j+k)
232         continue
231      continue
230   continue
#endif

#ifdef IMSL
      stop 'Not implemented yet!'
#endif
#ifdef LAPACK
#ifdef DP
      call dgetrf(dim*npart,dim*npart,zg,dim*npart,indx,i)
#else
      call sgetrf(dim*npart,dim*npart,zg,dim*npart,indx,i)
#endif
      if (i .ne. 0) print *,'Error during factorization 1! (',i,')'
#endif

#ifdef MOBILITY
#ifdef IMSL
      stop 'Not implemented yet!'
#endif
#ifdef LAPACK
#ifdef DP
      call dgetrs('N',dim*npart,6*npart,zg,dim*npart,indx,mp,
     .     dim*npart,j)
#else
      call sgetrs('N',dim*npart,6*npart,zg,dim*npart,indx,mp,
     .     dim*npart,j)
#endif
      if (j .ne. 0) print *,'Error during solution 1! (',j,')'
#endif
      do 300 i = 0,npart-1
         do 310 j = 0,npart-1
            do 320 k = 1,6
               do 330 l = 1,6
                  mo(6*j+l,6*i+k) = mp(j*dim+ti(l),6*i+ti(k))
330            continue
320         continue
310      continue
300   continue

      do 400 i = 0,npart-1
         do 401 j = 0,npart-1
            do 402 k = 1,3
               do 403 l = 1,3
                  mo(6*i+l,6*j+k) = mo(6*i+l,6*j+k)/dabs(z(0,1))
                  mo(6*i+l+3,6*j+k+3) = mo(6*i+l+3,6*j+k+3)
     .                 /(4.d0*dabs(z(0,2)))
                  mo(6*i+l,6*j+k+3) = mo(6*i+l,6*j+k+3)
     .                 /(4.d0*dabs(z(0,2)*z(0,1)))
403            continue
402         continue
401      continue
400   continue
#endif

#ifdef VELOCITIES
#ifdef IMSL
      stop 'Not implemented yet!'
#endif
#ifdef LAPACK
#ifdef DP
      call dgetrs('N',dim*npart,1,zg,dim*npart,indx,mpf,dim*npart,j)
#else
      call sgetrs('N',dim*npart,1,zg,dim*npart,indx,mpf,dim*npart,j)
#endif
      if (j .ne. 0) print *,'Error during solution 1! (',j,')'
#endif
      do 350 i = 0,npart-1
         do 351 k = 1,3
            v(6*i+k) = mpf(dim*i+ti(k))/dabs(z(0,1))
351      continue
         do 352 k = 4,6
            v(6*i+k) = mpf(dim*i+ti(k))
352      continue
350   continue
#endif

#ifdef NORMALIZE
      call norm(npart)
#endif
 
      end
#endif


* Invert friction matrix to get mobility matrix

#ifdef FRICTION
#ifdef MOBILITY
      subroutine inv(np)
      implicit none
      integer np

      integer npart
      parameter (npart = _NP_)
      
      double precision fr(6*npart,6*npart)
      common /frict/ fr
      
      double precision mo(6*npart,6*npart)
      common /mobil/ mo

      integer i,j

      do 300 i = 1,6*np
         do 301 j = 1,6*np
            mo(i,j) = fr(i,j)
301      continue
300   continue

#ifdef IMSL
#ifdef DP
      call dlinds(6*np,mo,6*npart,mo,6*npart)
#else
      call linds(6*np,mo,6*npart,mo,6*npart)
#endif
#endif
#ifdef LAPACK
#ifdef DP
      call dpotrf('U',6*np,mo,6*npart,i)
      if (i .ne. 0) print *,'Error during factorization 2! (',i,')'
      call dpotri('U',6*np,mo,6*npart,i)
      if (i .ne. 0) print *,'Error during inversion 2! (',i,')'
#else
      call spotrf('U',6*np,mo,6*npart,i)
      if (i .ne. 0) print *,'Error during factorization 2! (',i,')'
      call spotri('U',6*np,mo,6*npart,i)
      if (i .ne. 0) print *,'Error during inversion 2! (',i,')'
#endif
#endif

      do 320 i = 1,6*np
         do 321 j = 1,i-1
            mo(i,j) = mo(j,i)
321      continue
320   continue

      end
#endif
#endif


#ifdef NORMALIZE
      subroutine norm(np)
      implicit none
      integer np

      integer npart
      parameter (npart = _NP_)

#ifdef FRICTION
      double precision fr(6*npart,6*npart)
      common /frict/ fr
#endif
      
#ifdef MOBILITY
      double precision mo(6*npart,6*npart)
      common /mobil/ mo
#endif

      integer lm
      parameter (lm = _LM_)

      double precision z(0:lm,4),zs(0:lm,4)
      common /zm/ z,zs

      integer i,j,k,l

      do 330 i = 0,np-1
         do 331 j = 0,np-1
            do 332 k = 1,3
               do 333 l = 1,3
#ifdef FRICTION
                  fr(6*i+l,6*j+k) = fr(6*i+l,6*j+k)/dabs(z(0,1))
                  fr(6*i+l+3,6*j+k+3) = fr(6*i+l+3,6*j+k+3)/
     .                 (4.d0*dabs(z(0,2)))
#endif
#ifdef MOBILITY
                  mo(6*i+l,6*j+k) = mo(6*i+l,6*j+k)*dabs(z(0,1))
                  mo(6*i+l+3,6*j+k+3) = mo(6*i+l+3,6*j+k+3)*
     .                 (4.d0*dabs(z(0,2)))
#endif
333            continue
332         continue
331      continue
330   continue

      end
#endif
      
      
* Solve for velocities given the friction matrix and acting forces
      
#ifdef FRICTION
#ifdef VELOCITIES
      subroutine solvev
      implicit none
      
      integer npart
      parameter (npart = _NP_)
      
      double precision fr(6*npart,6*npart)
      common /frict/ fr
      
      double precision f(6*npart),v(6*npart)
      common /fv/ f,v
      
#ifdef RIGID
      double precision vrb(6*npart,6*npart)
      integer nrb,irb(npart)
      common /rb/ vrb,nrb,irb
#endif
      
      double precision fmat(6*npart,6*npart)
      common fmat
      
      double precision fv(6*npart)
      integer np
      integer i,j
      
#ifdef RIGID
      np = nrb
#else
      np = npart
#endif
      
      do 100 i = 1,6*np
         fv(i) = f(i)
100   continue
#ifdef RIGID
#ifndef RB_VELOCITIES
      if (nrb .ne. npart) then
         do 110 i = 1,6*nrb
            fv(i) = 0.d0
            do 111 j = 1,6*npart
               fv(i) = fv(i) + vrb(j,i)*f(j)
111         continue
110      continue
      endif
#endif
#endif
      
      do 120 i = 1,6*np
         do 121 j = 1,6*np
            fmat(i,j) = fr(i,j)
121      continue
120   continue
      
#ifdef IMSL
#ifdef DP
      call dlftds(6*np,fmat,6*npart,fmat,6*npart)
      call dlfsds(6*np,fmat,6*np,fv,fv)
#else
      call lftds(6*np,fmat,6*npart,fmat,6*npart)
      call lfsds(6*np,fmat,6*np,fv,fv)
#endif
#endif
#ifdef LAPACK
#ifdef DP
      call dpotrf('U',6*np,fmat,6*npart,i)
      if (i .ne. 0) print *,'Error during factorization 3! (',i,')'
      call dpotrs('U',6*np,1,fmat,6*npart,fv,6*npart,i)
#else
      call spotrf('U',6*np,fmat,6*npart,i)
      if (i .ne. 0) print *,'Error during factorization 3! (',i,')'
      call spotrs('U',6*np,1,fmat,6*npart,fv,6*npart,i)
#endif
      if (i .ne. 0) print *,'Error during solution 3! (',i,')'
#endif
      
      do 130 i = 1,6*np
         v(i) = fv(i)
130   continue
#ifdef RIGID
#ifndef RB_VELOCITIES
      if (nrb .ne. npart) then
         do 140 i = 1,6*npart
            v(i) = 0.d0
            do 141 j = 1,6*nrb
               v(i) = v(i) + vrb(i,j)*fv(j)
141         continue
140      continue
      endif
#endif
#endif
      
      end
#endif
#endif
