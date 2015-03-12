********************************************
* Friction coefficients for a linear chain *
* of touching spheres.                     *
* Example application for HYDROLIB         *
*                                          *
* K. Hinsen                                *
********************************************


* Include the configuration file
#include "example.h"

* Include HYDROLIB routines
#include <subr.f>

      
      program chain
      
      integer npart
      parameter (npart = _NP_)
      
* The arrays that stores the particle configuration.
* c contains the coordinates,
* cnct defines the rigid units.
#ifdef DPSOURCE
      double precision c(0:2,npart)
#else
      real c(0:2,npart)
#endif
      character*1 cnct(npart,npart)
      common /conf/ c,cnct
      
* The array for the friction matrix.
#ifdef DPSOURCE
      double precision fr(6*npart,6*npart)
#else
      real fr(6*npart,6*npart)
#endif
      common /frict/ fr

* Other variables
      integer i,j
#ifdef DPSOURCE
      double precision zpar,zpp,zraxis,zrpp
#else
      real zpar,zpp,zraxis,zrpp
#endif
      

* Library initialization
      call init


* Set up configuration. The particles are lined
* up along the x-axis.

      do 110 j = 1,npart
         c(0,j) = 2.*(j-1)
         c(1,j) = 0.
         c(2,j) = 0.
110   continue

* Definition of rigid connections. It would be sufficient to connect
* every particle to its neighbours, but specifying connections between
* all pairs is equivalent for this case and a bit easier.
      do 130 i = 1,npart
         do 131 j = 1,npart
            cnct(i,j) = '*'
131      continue
130   continue

* Library call: evaluate friction matrix
      call eval


* Pick out the four distinct friction coefficients for a linear chain.

* Translational motion along the axis
      zpar = fr(1,1)
* Translational motion perpendicular to the axis
      zpp = fr(2,2)
* Rotation around the axis
      zraxis = fr(4,4)
* Rotation perpendicular to the axis
      zrpp = fr(5,5)

* Print results
      print *,'N=',npart
      print *,zpar,zpp
      print *,zraxis,zrpp
      
      end
