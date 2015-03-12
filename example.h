/* Configuration */


/* numerical parameters */

#define _NP_ 6        /* number of particles in the system */
#define _LM_ 2        /* multipole order to be used */


/* compile-time options */

#undef  PERIODIC      /* periodic boundary conditions */
#define RIGID         /* enable treatment of rigid clusters */
#define LUBRICATION   /* short-range corrections */

#define FRICTION      /* calculate friction matrix */
#undef  MOBILITY      /* calculate mobility matrix */
#undef  VELOCITIES    /* calculate velocities for given forces */
#undef  RB_VELOCITIES /* calculate rigid-body velocities */

#undef  NORMALIZE     /* normalize friction and mobility matrices */

#define POSDEF        /* interaction matrix is positive definite */
#define PACKED        /* use packed storage mode (LAPACK only) */

#define DP            /* call double-precision library routines */


/* The following definitions must agree with those used
 * while compiling the library files
 */

#define LAPACK        /* use LAPACK for linear algebra */
