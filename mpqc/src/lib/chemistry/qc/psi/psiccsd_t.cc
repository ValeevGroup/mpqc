
#include "psicc.h"
#include "file11.h"

extern "C" {
#include <stdlib.h>
#include <math.h>
}

#include <util/keyval/keyval.h>
#include <util/misc/libmisc.h>
#include <math/scmat/matrix.h>
#include <math/symmetry/pointgrp.h>
#include <chemistry/molecule/molecule.h>


#define CLASSNAME PSI_CCSD_T
#define PARENTS public Wavefunction
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>

void *
PSI_CCSD_T::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = Wavefunction::_castdown(cd);
  return do_castdowns(casts,cd);
}

PSI_CCSD_T::PSI_CCSD_T(const RefKeyVal&keyval):
  Wavefunction(keyval), psi_in(keyval)
{
}

PSI_CCSD_T::~PSI_CCSD_T()
{
}

PSI_CCSD_T::PSI_CCSD_T(StateIn&s):
  Wavefunction(s)
  maybe_SavableState(s)
{
  abort();
}

void
PSI_CCSD_T::save_data_state(StateOut&s)
{
  Wavefunction::save_data_state(s);
  abort();
}

void
PSI_CCSD_T::print(ostream&o)
{
  Wavefunction::print(o);
}

void
PSI_CCSD_T::compute()
{  
  int i;
  int ni = molecule()->point_group()->char_table().nirrep();

  //molecule()->transform_to_principal_axes();

  if (!psi_in.test()){
    psi_in.write_input_file(do_gradient() ? "FIRST" : "NONE", "CCSD_T",
      (int)-log10(desired_value_accuracy()), "input.dat");

    system("inputth");
    system("psi");
 
  }

  // read output
  if (do_gradient()) {
    int i, j, ii;
    int *reorder;
    double tol = 1e-6;
    reorder = new int[molecule()->natom()];
    FILE11 file11(0);
    for(i=0; i<molecule()->natom(); i++)
      for(j=0; j<molecule()->natom(); j++)
        if( fabs(file11.coordinate(0,i) - molecule()->r(j,0)) < tol &&
            fabs(file11.coordinate(1,i) - molecule()->r(j,1)) < tol &&
            fabs(file11.coordinate(2,i) - molecule()->r(j,2)) < tol){
          reorder[i] = j;
          break;
          }
    RefSCVector g(moldim(),matrixkit());

    for (ii=0, i=0; i<molecule()->natom(); i++) {
      for (j=0; j<3; j++, ii++) {
        g(ii) = file11.gradient(j,reorder[i]);
        }
      }
    set_gradient(g);
    set_actual_gradient_accuracy(desired_gradient_accuracy());
    set_energy(file11.energy());
    }
  else {
      // read the energy from the output file
      FILE *in;
      system("rm -f psitmp.energy");
      if(!psi_in.test()){
        system("grep \"CCT \" energy.dat > psitmp.energy");
        }
      in = fopen("psitmp.energy","r");
      if (!in) {
          fprintf(stderr,"PSI_CCSD_T::compute(): cannot open psitmp.energy\n");
          abort();
        }
      double r;
      while(fscanf(in,"%*s %lf", &r)>0){ }
      set_energy(r);
      fclose(in);
    }
}

double
PSI_CCSD_T::density(cart_point&d)
{
  abort();
}

RefSymmSCMatrix
PSI_CCSD_T::density()
{
  abort();
}

int
PSI_CCSD_T::spin_polarized()
{
  return 1;
}

int
PSI_CCSD_T::nelectron()
{
  abort();
  return 0;
}

int
PSI_CCSD_T::value_implemented()
{
  return 1;
}

int
PSI_CCSD_T::gradient_implemented()
{
  return 1;
}
