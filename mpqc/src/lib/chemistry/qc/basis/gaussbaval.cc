
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
    int icenter = 0;
    int ishell = 0;
    int ibasis = 0;
    int nreturns;

    // calculate the value of each basis
    for (Pix center=first_center(); center != 0;next_center(center)) 
    {
	// Calculate r_diff
	r_diff.x()=r.x()-r_center(center,0);
	r_diff.y()=r.y()-r_center(center,1);
	r_diff.z()=r.z()-r_center(center,2);

#ifdef EXTRA_PRINT
        static int iflag=0;
	if (iflag)
	{
	    iflag--;
	    printf("Center %d, (%lf,%lf,%lf)\n",icenter,r_center(center,0),
		   r_center(center,1),r_center(center,2));
	}
#endif

	for (Pix shell=first_shell_on_center(center);
	     shell != 0;
	     next_shell_on_center(center,shell)) 
	{
            if (basis_values && g_values)
            {
	        nreturns=operator[](shell).grad_values(r_diff,
                                                  &g_values[ibasis*3],
                                                  &basis_values[ibasis]);
            }
            else if (g_values)
            {
	        nreturns=operator[](shell).grad_values(r_diff,
                                                  &g_values[ibasis*3],
                                                  &basis_values[ibasis]);
            }
            else if (basis_values)
            {
	        nreturns=operator[](shell).grad_values(r_diff,
                                                  0,
                                                  &basis_values[ibasis]);
            }
            ibasis += nreturns;
	    ishell++;
	}
	icenter++;
    }

    return ibasis;
}
