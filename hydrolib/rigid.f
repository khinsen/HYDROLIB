#ifdef DPSOURCE
#include <rigid_dp.f>
#else
************************************
* Friction and mobility            *
* for a hard-sphere configuration  *
* Module RIGID                     *
*                                  *
* K. Hinsen                        *
* Last revision: June 16, 1994     *
************************************


* Calculate connection matrix from rigid body list

      subroutine cnctm(cl)
      implicit none
      integer cl(*)

      integer npart
      parameter (npart = _NP_)
      
      real c(0:2,npart)
      character*1 cnct(npart,npart)
      common /conf/ c,cnct

      integer i,j

      do 100 i = 1,npart
         do 101 j = 1,npart
            cnct(i,j) = ' '
            if (cl(i) .eq. cl(j)) cnct(i,j) = '+'
101      continue
100   continue
      
      end


* construct list of rigid bodies from connection matrix

      subroutine rgdbd(cl,ncl)
      implicit none
      integer cl(*),ncl

      integer npart
      parameter (npart = _NP_)
      
      real c(0:2,npart)
      character*1 cnct(npart,npart)
      common /conf/ c,cnct

      integer i,j,k,no,nn

      do 100 i = 1,npart
         cl(i) = 0
100   continue
      ncl = 0
      do 110 i = 1,npart
         if (cl(i) .eq. 0) then
            cl(i) = i
            ncl = ncl + 1
         endif
         do 111 j = i+1,npart
            if (cnct(i,j) .ne. ' ' .and. cl(i) .ne. cl(j)) then
               if (cl(j) .eq. 0) then
                  cl(j) = cl(i)
               else
                  no = cl(j)
                  nn = cl(i)
                  do 112 k = 1,npart
                     if (cl(k) .eq. no) cl(k) = nn
112               continue
                  ncl = ncl - 1
               endif
            endif
111      continue
110   continue
      
      call cnctm(cl)

      end
      
      
* calculate rigid body velocity matrix
      
      subroutine velmat
      implicit none
      
      integer npart
      parameter (npart = _NP_)
      
      real c(0:2,npart)
      character*1 cnct(npart,npart)
      common /conf/ c,cnct
      
      real vrb(6*npart,6*npart)
      integer nrb,irb(npart)
      common /rb/ vrb,nrb,irb
      
      integer cl(npart),np(npart)
      real re(0:2,npart)
      
      integer i,j,k
      real x,y,z
      
      call rgdbd(cl,nrb)
      do 100 i = 1,npart
         np(i) = 0
         re(0,i) = 0.
         re(1,i) = 0.
         re(2,i) = 0.
100   continue
      do 110 i = 1,npart
         j = cl(i)
         np(j) = np(j) + 1
         re(0,j) = re(0,j) + c(0,i)
         re(1,j) = re(1,j) + c(1,i)
         re(2,j) = re(2,j) + c(2,i)
110   continue
      do 120 i = 1,npart
         j = cl(i)
         re(0,i) = re(0,j)
         re(1,i) = re(1,j)
         re(2,i) = re(2,j)
120   continue
      do 130 i = 1,npart
         j = cl(i)
         re(0,i) = re(0,i)/np(j)
         re(1,i) = re(1,i)/np(j)
         re(2,i) = re(2,i)/np(j)
130   continue
      
      j = 1
      do 140 i = 1,npart
         if (np(i) .gt. 0) then
            np(i) = j
            j = j + 1
         endif
140   continue
      do 150 i = 1,npart
         irb(i) = np(cl(i))
150   continue
      
      do 200 i = 1,6*npart
         do 201 j = 1,6*nrb
            vrb(i,j) = 0.
201      continue
200   continue
      do 210 i = 0,npart-1
         j = irb(i+1)-1
         do 211 k = 1,6
            vrb(6*i+k,6*j+k) = 1.
211      continue
         x = c(0,i+1)-re(0,i+1)
         y = c(1,i+1)-re(1,i+1)
         z = c(2,i+1)-re(2,i+1)
         vrb(6*i+2,6*j+6) = x
         vrb(6*i+3,6*j+5) = -x
         vrb(6*i+3,6*j+4) = y
         vrb(6*i+1,6*j+6) = -y
         vrb(6*i+1,6*j+5) = z
         vrb(6*i+2,6*j+4) = -z
210   continue
      
      end
      
      
* calculate reduced friction and mobility matrices
      
      subroutine evalfr
      implicit none
      
      integer npart
      parameter (npart = _NP_)

      real vrb(6*npart,6*npart)
      integer nrb,irb(npart)
      common /rb/ vrb,nrb,irb
      
      real fr(6*npart,6*npart)
      common /frict/ fr
      
      real frrb(6*npart,6*npart)
      common frrb
      
      integer i,j,k,l
      
      call velmat

#ifdef MOBILITY_DIRECT
      if (nrb .eq. npart) then
         call evalm
      else
#endif

      call evalf
      
      if (nrb .ne. npart) then
         do 100 i = 1,6*nrb
            do 101 j = 1,6*nrb
               frrb(i,j) = 0.
               do 102 k = 1,6*npart
                  do 103 l = 1,6*npart
                     frrb(i,j) = frrb(i,j) + vrb(k,i)*fr(k,l)*vrb(l,j)
103               continue
102            continue
101         continue
100      continue
         
         do 110 i = 1,6*nrb
            do 111 j = 1,6*nrb
               fr(i,j) = frrb(i,j)
111         continue
110      continue
         
      endif
      
#ifdef VELOCITIES_RB
      call solvev
#endif

#ifdef MOBILITY_RB
      call inv(nrb)
#endif

#ifdef NORMALIZE_RB
      call norm(nrb)
#endif

#ifdef MOBILITY_DIRECT
      endif
#endif

      end
#endif
