//
// readvect.cc
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

#include <chemistry/qc/mpqc/mpqc.h>
#include <math/newmat7/newmat.h>

void MPSCF::read_vector(char *fname,int n_basis, RefSCMatrix &scf_vector)
{
    // Try to open file containing scf vector
    FILE *fd= fopen(fname,"r");
    if(fd == NULL)
    {
        fprintf(stderr,"Couldn't open %s to read \"%s\"\n",
                fname);
        exit(1);
    }

    // Read in vector
    double *scf_tmp=new double[n_basis];
    if (!scf_tmp)
    {
        fprintf(stderr," Could not allocate scf_tmp in read_vector\n");
        exit(1);
    }

    for (int i=0; i<n_basis; i++)
    {
        int readval=0; // Tmp for holding # elements read in

        readval= fread(scf_tmp, sizeof(double), n_basis, fd);
        if (readval != n_basis)
        {
            fprintf(stderr," Error in reading vector nbasis=%d, readval=%d",
                    n_basis, readval);
            exit(1);
        }

        // Copy this row to scf_vector matrix
        //printf(" MO #%d\n",i);
        for (int j=0; j< n_basis; j++)
        {
            scf_vector.set_element(j,i,scf_tmp[j]);
            //printf(" %d %lf\n",j,scf_tmp[j]);
        }
    }
    delete[] scf_tmp;

    _have_eigenvectors = 1;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
// End:
