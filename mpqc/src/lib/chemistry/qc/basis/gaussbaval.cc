
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <util/keyval/keyval.h>

#include <math/topology/point.h>
#include "gaussbas.h"
#include "gaussshell.h"

int GaussianBasisSet::values(cart_point& r, double* basis_values) const
{
  return grad_values(r, 0, basis_values);
}

int GaussianBasisSet::grad_values(cart_point& r,
                             double* g_values,
                             double* basis_values) const
{
    cart_point r_diff;
    int ishell = 0;
    int ibasis = 0;
    int nreturns;

    // calculate the value of each basis
    for (int icenter=0; icenter < ncenter_; icenter++) 
      {
        int nshell = center_to_nshell_[icenter];

	// Calculate r_diff
	r_diff.x()=r.x()-center_to_r_(icenter,0);
	r_diff.y()=r.y()-center_to_r_(icenter,1);
	r_diff.z()=r.z()-center_to_r_(icenter,2);

#ifdef EXTRA_PRINT
        static int iflag=0;
	if (iflag)
	{
	    iflag--;
	    printf("Center %d, (%lf,%lf,%lf)\n",icenter,r_center(center,0),
		   r_center(center,1),r_center(center,2));
	}
#endif

	for (int ish=0; ish < nshell; ish++) {
            if (basis_values && g_values)
              {
#ifdef __GNUC__
	        nreturns=operator()(ishell).grad_values(r_diff,
#else
	        nreturns=((GaussianShell)operator()(ishell)).grad_values(r_diff,
#endif
                                                        &g_values[ibasis*3],
                                                        &basis_values[ibasis]);
              }
            else if (g_values)
              {
#ifdef __GNUC__
	        nreturns=operator()(ishell).grad_values(r_diff,
#else
	        nreturns=((GaussianShell)operator()(ishell)).grad_values(r_diff,
#endif
                                                        &g_values[ibasis*3],
                                                        &basis_values[ibasis]);
              }
            else if (basis_values)
              {
#ifdef __GNUC__
	        nreturns=operator()(ishell).grad_values(r_diff,
#else
	        nreturns=((GaussianShell)operator()(ishell)).grad_values(r_diff,
#endif
                                                        0,
                                                        &basis_values[ibasis]);
              }
            ibasis += nreturns;
	    ishell++;
	}
    }

    return ibasis;
}
