
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


#define CLASSNAME PSI_CCSD
#define PARENTS public Wavefunction
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>

void *
PSI_CCSD::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = Wavefunction::_castdown(cd);
  return do_castdowns(casts,cd);
}

PSI_CCSD::PSI_CCSD(KeyVal&keyval):
  Wavefunction(keyval), psi_in(keyval)
{
}

PSI_CCSD::~PSI_CCSD()
{
}

PSI_CCSD::PSI_CCSD(StateIn&s):
  SavableState(s,class_desc_),
  Wavefunction(s)
{
  abort();
}

void
PSI_CCSD::save_data_state(StateOut&s)
{
  Wavefunction::save_data_state(s);
  abort();
}

void
PSI_CCSD::print(SCostream&o)
{
  Wavefunction::print(o);
}

void
PSI_CCSD::compute()
{  
  int i;
  int ni = _mol->point_group().char_table().nirrep();

  mol_transform_to_principal_axes(*_mol.pointer());

  if (!psi_in.test()){
    psi_in.write_input_file(_gradient.needed() ? "FIRST" : "NONE", "CCSD",
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
        system("grep \"CCSD \" energy.dat > psitmp.energy");
        }
      in = fopen("psitmp.energy","r");
      if (!in) {
          fprintf(stderr,"PSI_CCSD::compute(): cannot open psitmp.energy\n");
          abort();
        }
      double r;
      while(fscanf(in,"%*s %lf", &r)>0){ }
      set_energy(r);
      fclose(in);
    }
}

double
PSI_CCSD::density(cart_point&d)
{
  abort();
}

RefSymmSCMatrix
PSI_CCSD::density()
{
  abort();
}
