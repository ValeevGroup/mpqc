
#ifdef __GNUC__
#pragma implementation
#endif

#include <iostream.h>
#include <iomanip.h>

#include <util/misc/formio.h>

#include <math/scmat/local.h>
#include <math/scmat/repl.h>
#include <math/scmat/dist.h>

#include <chemistry/qc/scf/scf.h>

///////////////////////////////////////////////////////////////////////////
// SCF

#define CLASSNAME SCF
#define PARENTS public OneBodyWavefunction
#include <util/class/classia.h>
void *
SCF::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = OneBodyWavefunction::_castdown(cd);
  return do_castdowns(casts,cd);
}

SCF::SCF(StateIn& s) :
  OneBodyWavefunction(s)
  maybe_SavableState(s)
{
  s.get(maxiter_);
  s.get(int_store_);
  s.get(dens_reset_freq_);
  s.get(level_shift_);

  integral()->set_storage(int_store_);
  scf_grp_ = basis()->matrixkit()->messagegrp();
}

SCF::SCF(const RefKeyVal& keyval) :
  OneBodyWavefunction(keyval),
  maxiter_(40),
  int_store_(0),
  level_shift_(0),
  dens_reset_freq_(10)
{
  if (keyval->exists("maxiter"))
    maxiter_ = keyval->intvalue("maxiter");

  if (keyval->exists("density_reset_frequency"))
    dens_reset_freq_ = keyval->intvalue("density_reset_frequency");

  if (keyval->exists("integral_storage"))
    int_store_ = keyval->intvalue("integral_storage");

  if (keyval->exists("level_shift"))
    level_shift_ = keyval->doublevalue("level_shift");

  integral()->set_storage(int_store_);

  scf_grp_ = basis()->matrixkit()->messagegrp();
}

SCF::~SCF()
{
}

void
SCF::save_data_state(StateOut& s)
{
  OneBodyWavefunction::save_data_state(s);
  s.put(maxiter_);
  s.put(int_store_);
  s.put(dens_reset_freq_);
  s.put(level_shift_);
}

RefSCMatrix
SCF::eigenvectors()
{
  return eigenvectors_.result();
}

void
SCF::print(ostream&o)
{
  OneBodyWavefunction::print(o);
  if (scf_grp_->me()==0) {
    o << indent << "SCF Parameters:\n" << incindent;
    o << indent << "maxiter = " << maxiter_ << endl;
    o << indent << "integral_storage = " << int_store_ << endl;
    o << indent << "density_reset_freq = " << dens_reset_freq_ << endl;
    o << indent << "level_shift = " << level_shift_ << endl;
    o << decindent << endl;
  }
}

void
SCF::compute()
{
  int me=scf_grp_->me();
  
  if (hessian_needed())
    set_desired_gradient_accuracy(desired_hessian_accuracy()/100.0);

  if (gradient_needed())
    set_desired_value_accuracy(desired_gradient_accuracy()/100.0);

  if (value_needed()) {
    if (me==0)
      cout << endl << indent << "SCF::compute: energy accuracy = " <<
        desired_value_accuracy() << endl << endl;

    double eelec;
    compute_vector(eelec);
      
    // this will be done elsewhere eventually
    double nucrep = molecule()->nuclear_repulsion_energy();
    if (me==0)
      cout << endl << indent << "total scf energy = " <<
        setw(20) << setprecision(15) << eelec+nucrep << endl;

    set_energy(eelec+nucrep);
    set_actual_value_accuracy(desired_value_accuracy());
  }

  if (gradient_needed()) {
    RefSCVector gradient = matrixkit()->vector(moldim());

    if (me==0)
      cout << endl << indent << "SCF::compute: gradient accuracy = " <<
        desired_gradient_accuracy() << endl << endl;

    compute_gradient(gradient);
    gradient.print("cartesian gradient");
    set_gradient(gradient);

    set_actual_gradient_accuracy(desired_gradient_accuracy());
  }
  
  if (hessian_needed()) {
    RefSymmSCMatrix hessian = matrixkit()->symmmatrix(moldim());
    
    if (me==0)
      cout << endl << indent << "SCF::compute: hessian accuracy = " <<
        desired_hessian_accuracy() << endl << endl;

    compute_hessian(hessian);
    set_hessian(hessian);

    set_actual_hessian_accuracy(desired_hessian_accuracy());
  }
}
