
/* The return codes from the geom library. */

#define GEOM_COMPUTE_ENERGY   1 /* Compute an energy and give it to libgeom. */
#define GEOM_COMPUTE_GRADIENT 2 /* Compute a gradient and give it to libgeom. */
#define GEOM_COMPUTE_HESSIAN  3 /* Compute a hessian and give it to libgeom. */
#define GEOM_ABORT            4 /* Abort the optimization. */
#define GEOM_DONE             5 /* The geometry is optimized. */
#define GEOM_NOTDONE          6 /* The geometry is not optimized. */

