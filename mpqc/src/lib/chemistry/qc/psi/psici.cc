
#include "psici.h"
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


#define CLASSNAME PSI_CI
#define PARENTS public Wavefunction
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>

void *
PSI_CI::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = Wavefunction::_castdown(cd);
  return do_castdowns(casts,cd);
}

PSI_CI::PSI_CI(const RefKeyVal&keyval):
  Wavefunction(keyval), psi_in(keyval)
{
}

PSI_CI::~PSI_CI()
{
}

PSI_CI::PSI_CI(StateIn&s):
  Wavefunction(s)
  maybe_SavableState(s)
{
  abort();
}

void
PSI_CI::save_data_state(StateOut&s)
{
  Wavefunction::save_data_state(s);
  abort();
}

void
PSI_CI::print(ostream&o)
{
  Wavefunction::print(o);
}

void
PSI_CI::compute()
{  
  int i;
  int ni = _mol->point_group().char_table().nirrep();

  mol_transform_to_principal_axes(_mol);

  if (!psi_in.test()){
    psi_in.write_input_file(_gradient.needed() ? "FIRST" : "NONE", "CI",
      (int)-log10(desired_value_accuracy()), "input.dat");

    system("inputth");
    system("psi");
 
  }

  // read output
  if (_gradient.needed()) {
    int i, j, ii;
    int *reorder;
    double tol = 1e-6;
    reorder = new int[_mol->natom()];
    FILE11 file11(0);
    for(i=0; i<_mol->natom(); i++)
      for(j=0; j<_mol->natom(); j++)
        if( fabs(file11.coordinate(0,i) - _mol->atom(j)[0]) < tol &&
            fabs(file11.coordinate(1,i) - _mol->atom(j)[1]) < tol &&
            fabs(file11.coordinate(2,i) - _mol->atom(j)[2]) < tol){
          reorder[i] = j;
          break;
          }
    RefSCVector g(_moldim);

    for (ii=0, i=0; i<_mol->natom(); i++) {
      for (j=0; j<3; j++, ii++) {
        g(ii) = file11.gradient(j,reorder[i]);
        }
      }
    set_gradient(g);
    _gradient.set_actual_accuracy(desired_gradient_accuracy());
    set_energy(file11.energy());
    }
  else {
      // read the energy from the output file
      FILE *in;
      system("rm -f psitmp.energy");
      if(!psi_in.test()){
        system("grep \"1 ECI \" output.dat > psitmp.energy");
        }
      in = fopen("psitmp.energy","r");
      if (!in) {
          fprintf(stderr,"PSI_CI::compute(): cannot open psitmp.energy\n");
          abort();
        }
      double r;
      fscanf(in,"%*s %*s %*s %*s %lf", &r);
      set_energy(r);
      fclose(in);
    }
}

double
PSI_CI::density(cart_point&d)
{
  abort();
}

RefSymmSCMatrix
PSI_CI::density()
{
  abort();
}
