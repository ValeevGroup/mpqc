
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <util/keyval/keyval.h>

#include <math/topology/point.h>
#include <chemistry/qc/basis/gaussbas.h>
#include <chemistry/qc/basis/gaussshell.h>

int
GaussianBasisSet::values(const RefIntegral& ints,
                         const SCVector3& r, double* basis_values) const
{
  return grad_values(ints, r, 0, basis_values);
}

int
GaussianBasisSet::grad_values(const RefIntegral& ints,
                              const SCVector3& r,
                              double* g_values,
                              double* basis_values) const
{
    SCVector3 r_diff;
    int ishell = 0;
    int ibasis = 0;
    int nreturns;

    // for convenience
    const GaussianBasisSet& gbs = *this;

    // calculate the value of each basis
    for (int icenter=0; icenter < ncenter_; icenter++) 
      {
        int nshell = center_to_nshell_[icenter];

	// Calculate r_diff
	r_diff.x()=r.x()-GaussianBasisSet::r(icenter,0);
	r_diff.y()=r.y()-GaussianBasisSet::r(icenter,1);
	r_diff.z()=r.z()-GaussianBasisSet::r(icenter,2);

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
	        nreturns=gbs(ishell).grad_values(ints, r_diff,
                                                 &g_values[ibasis*3],
                                                 &basis_values[ibasis]);
              }
            else if (g_values)
              {
	        nreturns=gbs(ishell).grad_values(ints, r_diff,
                                                 &g_values[ibasis*3],
                                                 &basis_values[ibasis]);
              }
            else if (basis_values)
              {
	        nreturns=gbs(ishell).grad_values(ints, r_diff, 0,
                                                 &basis_values[ibasis]);
              }
            ibasis += nreturns;
	    ishell++;
	}
    }

    return ibasis;
}
