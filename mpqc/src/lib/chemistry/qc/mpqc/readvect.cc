
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

