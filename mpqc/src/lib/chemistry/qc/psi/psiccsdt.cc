
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


#define CLASSNAME PSI_CCSDT
#define PARENTS public Wavefunction
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>

void *
PSI_CCSDT::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = Wavefunction::_castdown(cd);
  return do_castdowns(casts,cd);
}

PSI_CCSDT::PSI_CCSDT(const RefKeyVal&keyval):
  Wavefunction(keyval), psi_in(keyval)
{
}

PSI_CCSDT::~PSI_CCSDT()
{
}

PSI_CCSDT::PSI_CCSDT(StateIn&s):
  Wavefunction(s)
  maybe_SavableState(s)
{
  abort();
}

void
PSI_CCSDT::save_data_state(StateOut&s)
{
  Wavefunction::save_data_state(s);
  abort();
}

void
PSI_CCSDT::print(ostream&o)
{
  Wavefunction::print(o);
}

void
PSI_CCSDT::compute()
{  
  int i;
  int ni = molecule()->point_group()->char_table().nirrep();

  //molecule()->transform_to_principal_axes();

  double energy_acc = desired_value_accuracy();
  double grad_acc = desired_gradient_accuracy();
  if (energy_acc > 1.0e-6) energy_acc = 1.0e-6;
  if (grad_acc > 1.0e-7) grad_acc = 1.0e-7;
  if (do_gradient() && energy_acc > grad_acc/10.0) energy_acc = grad_acc/10.0;

  if (!psi_in.test()){
    psi_in.write_input_file(do_gradient() ? "FIRST" : "NONE", "CCSDT",
      (int)-log10(energy_acc), "input.dat");

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
    set_actual_gradient_accuracy(grad_acc);
    set_energy(file11.energy());
    }
  else {
      // read the energy from the output file
      FILE *in;
      system("rm -f psitmp.energy");
      if(!psi_in.test()){
        system("grep \"FSDT\" energy.dat > psitmp.energy");
        }
      in = fopen("psitmp.energy","r");
      if (!in) {
          fprintf(stderr,"PSI_CCSDT::compute(): cannot open psitmp.energy\n");
          abort();
        }
      double r;
      while(fscanf(in,"%*s %lf", &r)>0){ }
      set_energy(r);
      fclose(in);
    }
  set_actual_value_accuracy(energy_acc);
}

double
PSI_CCSDT::density(cart_point&d)
{
  abort();
}

RefSymmSCMatrix
PSI_CCSDT::density()
{
  abort();
}

int
PSI_CCSDT::spin_polarized()
{
  return 1;
}

int
PSI_CCSDT::nelectron()
{
  abort();
  return 0;
}

int
PSI_CCSDT::value_implemented()
{
  return 1;
}

int
PSI_CCSDT::gradient_implemented()
{
  return 1;
}
