
#include <stdio.h>
#include <math/array/math_lib.h>
#include <chemistry/qc/intv2/int_libv2.h>
#include <chemistry/molecule/molecule.h>
#include <math/topology/pointbag.h>
#include "integral.h"
#include "integralv2.h"

OneBodyIntv2::OneBodyIntv2(const GaussianBasisSet*bs,const Molecule*mol):
  OneBodyInt(bs)
{
  c = bs->convert_to_centers_t(mol);
  int_normalize_centers(c);
  int_initialize_1e(0,0,c,c);
  int_initialize_offsets1(c,c);
}

OneBodyIntv2::~OneBodyIntv2()
{
  int_done_offsets1(c,c);
  int_done_1e();
  free_centers(c);
  free(c);
}

GaussianOverlapIntv2::GaussianOverlapIntv2(const GaussianBasisSet*bs_):
  OneBodyIntv2(bs_,0)
{
}

GaussianOverlapIntv2::~GaussianOverlapIntv2()
{
}

void GaussianOverlapIntv2::compute_shell(int i,int j,double*buf)
{
  int_shell_overlap(c,c,buf,i,j);
}

GaussianKineticIntv2::GaussianKineticIntv2(const GaussianBasisSet*bs_,
                                           const Molecule*m):
  OneBodyIntv2(bs_,m)
{
}

GaussianKineticIntv2::~GaussianKineticIntv2()
{
}

void GaussianKineticIntv2::compute_shell(int i,int j,double*buf)
{
  int_shell_kinetic(c,c,buf,i,j);
}

GaussianPointChargeIntv2::GaussianPointChargeIntv2(PointBag_double*charges,
						   const GaussianBasisSet*bs_,
						   const Molecule*m):
  OneBodyIntv2(bs_,m)
{
  ncharge = charges->length();
  
  if (ncharge) {
      position = new double*[ncharge];
      charge = new double[ncharge];
    }
  
  int i = 0;
  for (Pix pix= charges->first(); pix!=0; charges->next(pix)) {
      position[i] = new double[3];
      charge[i] = charges->get(pix);
      for (int j=0; j<3; j++) {
	  position[i][j] = charges->point(pix)[j];
	}
      i++;
    }
  delete charges;
}

GaussianPointChargeIntv2::~GaussianPointChargeIntv2()
{
  for (int i=0; i<ncharge; i++) {
      delete position[i];
    }
  if (ncharge) {
      delete[] charge;
      delete[] position;
    }
}

void GaussianPointChargeIntv2::compute_shell(int i,int j,double*buf)
{
  int_shell_point_charge(c,c,buf,i,j,ncharge,charge,position);
}

GaussianNuclearIntv2::GaussianNuclearIntv2(const GaussianBasisSet*bs_,
				           const Molecule*m):
GaussianPointChargeIntv2(m->charges(),bs_,m)
{
}

GaussianNuclearIntv2::~GaussianNuclearIntv2()
{
}
