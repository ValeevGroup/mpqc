
#ifdef __GNUC__
#pragma implementation
#endif

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

  integral()->set_storage(int_store_);
}

SCF::SCF(const RefKeyVal& keyval) :
  OneBodyWavefunction(keyval),
  maxiter_(40),
  int_store_(0),
  dens_reset_freq_(10)
{
  if (keyval->exists("maxiter"))
    maxiter_ = keyval->intvalue("maxiter");

  if (keyval->exists("density_reset_frequency"))
    dens_reset_freq_ = keyval->intvalue("density_reset_frequency");

  if (keyval->exists("integral_storage"))
    int_store_ = keyval->intvalue("integral_storage");

  integral()->set_storage(int_store_);
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
}

void
SCF::compute()
{
  if (hessian_needed())
    set_desired_gradient_accuracy(desired_hessian_accuracy()/100.0);

  if (gradient_needed())
    set_desired_value_accuracy(desired_gradient_accuracy()/100.0);

  if (value_needed()) {
    printf("\n  SCF::compute: energy accuracy = %g\n\n",
           desired_value_accuracy());

    double eelec;
    compute_vector(eelec);
      
    // this will be done elsewhere eventually
    double nucrep = molecule()->nuclear_repulsion_energy();
    printf("\n  total scf energy = %20.15f\n",eelec+nucrep);

    set_energy(eelec+nucrep);
    set_actual_value_accuracy(desired_value_accuracy());
  }

  if (gradient_needed()) {
    RefSCVector gradient = matrixkit()->vector(moldim());

    printf("\n  SCF::compute: gradient accuracy = %g\n\n",
           desired_gradient_accuracy());

    compute_gradient(gradient);
    gradient.print("cartesian gradient");
    set_gradient(gradient);

    set_actual_gradient_accuracy(desired_gradient_accuracy());
  }
  
  if (hessian_needed()) {
    RefSymmSCMatrix hessian = matrixkit()->symmmatrix(moldim());
    
    printf("\n  SCF::compute: hessian accuracy = %g\n\n",
           desired_hessian_accuracy());

    compute_hessian(hessian);
    set_hessian(hessian);

    set_actual_hessian_accuracy(desired_hessian_accuracy());
  }
}
