//
// gaussbaval.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#include <stdlib.h>
#include <math.h>

#include <util/misc/formio.h>
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
	    cout << node0 << indent
                 << scprintf("Center %d, (%lf,%lf,%lf)\n",
                             icenter,r_center(center,0),
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

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
// End:
