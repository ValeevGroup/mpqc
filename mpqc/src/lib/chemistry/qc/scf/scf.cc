
#ifdef __GNUC__
#pragma implementation
#endif

#include <math.h>
#include <sys/stat.h>

#include <util/misc/formio.h>

#include <math/scmat/local.h>
#include <math/scmat/offset.h>

#include <math/optimize/diis.h>

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
  s.get(dens_reset_freq_);
  s.get(reset_occ_);
  s.get(local_dens_);
  s.get(storage_);
  s.get(level_shift_);

  extrap_.restore_state(s);

  scf_grp_ = basis()->matrixkit()->messagegrp();
}

SCF::SCF(const RefKeyVal& keyval) :
  OneBodyWavefunction(keyval),
  maxiter_(40),
  level_shift_(0),
  reset_occ_(0),
  local_dens_(1),
  storage_(0),
  dens_reset_freq_(10)
{
  if (keyval->exists("maxiter"))
    maxiter_ = keyval->intvalue("maxiter");

  if (keyval->exists("density_reset_frequency"))
    dens_reset_freq_ = keyval->intvalue("density_reset_frequency");

  if (keyval->exists("reset_occupations"))
    reset_occ_ = keyval->booleanvalue("reset_occupations");

  if (keyval->exists("level_shift"))
    level_shift_ = keyval->doublevalue("level_shift");

  extrap_ = keyval->describedclassvalue("extrap");
  if (extrap_.null())
    extrap_ = new DIIS;
  
  storage_ = keyval->intvalue("memory");
  
  if (keyval->exists("local_density"))
    local_dens_ = keyval->booleanvalue("local_density");
    
  // first see if guess_wavefunction is a wavefunction, then check to
  // see if it's a string.
  if (keyval->exists("guess_wavefunction")) {
    guess_wfn_ = keyval->describedclassvalue("guess_wavefunction");
    if (guess_wfn_.null()) {
      char *path = keyval->pcharvalue("guess_wavefunction");
      struct stat sb;
      if (path && stat(path, &sb)==0 && sb.st_size) {
        StateInBinXDR s(path);

        // reset the default matrixkit so that the matrices in the guess
        // wavefunction will match those in this wavefunction
        RefSCMatrixKit oldkit = SCMatrixKit::default_matrixkit();
        SCMatrixKit::set_default_matrixkit(basis()->matrixkit());

        guess_wfn_.restore_state(s);

        // go back to the original default matrixkit
        SCMatrixKit::set_default_matrixkit(oldkit);
        delete[] path;
      }
    }
  }
  
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
  s.put(dens_reset_freq_);
  s.put(reset_occ_);
  s.put(local_dens_);
  s.put(storage_);
  s.put(level_shift_);
  extrap_.save_state(s);
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
  o << node0 << indent << "SCF Parameters:\n" << incindent
    << indent << "maxiter = " << maxiter_ << endl
    << indent << "density_reset_freq = " << dens_reset_freq_ << endl
    << indent << scprintf("level_shift = %f\n",level_shift_)
    << decindent << endl;
}

//////////////////////////////////////////////////////////////////////////////

void
SCF::compute()
{
  int me=scf_grp_->me();
  
  local_ = (LocalSCMatrixKit::castdown(basis()->matrixkit().pointer())) ?1:0;
  
  if (hessian_needed())
    set_desired_gradient_accuracy(desired_hessian_accuracy()/100.0);

  if (gradient_needed())
    set_desired_value_accuracy(desired_gradient_accuracy()/100.0);

  if (value_needed()) {
    cout << node0 << endl << indent
         << scprintf("SCF::compute: energy accuracy = %10.7e\n\n",
                     desired_value_accuracy());

    double eelec;
    compute_vector(eelec);
      
    // this will be done elsewhere eventually
    double nucrep = molecule()->nuclear_repulsion_energy();
    cout << node0 << endl << indent
         << scprintf("total scf energy = %20.15f\n", eelec+nucrep);

    set_energy(eelec+nucrep);
    set_actual_value_accuracy(desired_value_accuracy());
  }

  if (gradient_needed()) {
    RefSCVector gradient = matrixkit()->vector(moldim());

    cout << node0 << endl << indent
         << scprintf("SCF::compute: gradient accuracy = %10.7e\n\n",
                     desired_gradient_accuracy());

    compute_gradient(gradient);
    gradient.print("cartesian gradient");
    set_gradient(gradient);

    set_actual_gradient_accuracy(desired_gradient_accuracy());
  }
  
  if (hessian_needed()) {
    RefSymmSCMatrix hessian = matrixkit()->symmmatrix(moldim());
    
    cout << node0 << endl << indent
         << scprintf("SCF::compute: hessian accuracy = %10.7e\n\n",
                     desired_hessian_accuracy());

    compute_hessian(hessian);
    set_hessian(hessian);

    set_actual_hessian_accuracy(desired_hessian_accuracy());
  }
}

//////////////////////////////////////////////////////////////////////////////

char *
SCF::init_pmax(double *pmat_data)
{
  double l2inv = 1.0/log(2.0);
  double tol = pow(2.0,-126.0);
  
  GaussianBasisSet& gbs = *basis().pointer();
  
  char * pmax = new char[i_offset(gbs.nshell())];

  int ish, jsh, ij;
  for (ish=ij=0; ish < gbs.nshell(); ish++) {
    int istart = gbs.shell_to_function(ish);
    int iend = istart + gbs(ish).nfunction();
    
    for (jsh=0; jsh <= ish; jsh++,ij++) {
      int jstart = gbs.shell_to_function(jsh);
      int jend = jstart + gbs(jsh).nfunction();
      
      double maxp=0, tmp;

      for (int i=istart; i < iend; i++) {
        int ijoff = i_offset(i);
        for (int j=jstart; j < ((ish==jsh) ? i+1 : jend); j++,ijoff++)
          if ((tmp=fabs(pmat_data[ijoff])) > maxp)
            maxp=tmp;
      }

      if (maxp <= tol)
        maxp=tol;

      pmax[ij] = (signed char) (log(maxp)*l2inv);
    }
  }

  return pmax;
}

//////////////////////////////////////////////////////////////////////////////

RefSymmSCMatrix
SCF::get_local_data(const RefSymmSCMatrix& m, double*& p, Access access)
{
  RefSymmSCMatrix l = m;
  
  if (!LocalSymmSCMatrix::castdown(l.pointer())) {
    RefSCMatrixKit k = new LocalSCMatrixKit;
    l = k->symmmatrix(m.dim());
    l->convert(m);

    if (access == Accum)
      l->assign(0.0);
  } else if (scf_grp_->n() > 1 && access==Accum) {
    l = m.clone();
    l.assign(0.0);
  }

  p = LocalSymmSCMatrix::castdown(l.pointer())->get_data();
  return l;
}

//////////////////////////////////////////////////////////////////////////////

void
SCF::initial_vector()
{
  // if guess_wfn_ is non-null then try to get a guess vector from it.
  // First check that the same basis is used...if not, then project the
  // guess vector into the present basis.
  // right now the check is crude...there should be an equiv member in
  // GaussianBasisSet
  if (guess_wfn_.nonnull()) {
    if (guess_wfn_->basis()->nbasis() == basis()->nbasis()) {
      cout << node0 << indent
           << "Using guess wavefunction as starting vector\n";

      // indent output of eigenvectors() call if there is any
      cout << incindent << incindent;
      eigenvectors_ = guess_wfn_->eigenvectors();
      cout << decindent << decindent;
    } else {
      cout << node0 << indent
           << "Projecting guess wavefunction into the present basis set\n";

      // indent output of projected_eigenvectors() call if there is any
      cout << incindent << incindent;
      eigenvectors_ = projected_eigenvectors(guess_wfn_);
      cout << decindent << decindent;
    }

    // we should only have to do this once, so free up memory used
    // for the old wavefunction
    guess_wfn_=0;

    cout << node0 << endl;
      
  } else {
    cout << node0 << indent << "Starting from core Hamiltonian guess\n\n";
    eigenvectors_ = hcore_guess();
  }
}

//////////////////////////////////////////////////////////////////////////////

void
SCF::init_mem(int nm)
{
  // if local_den_ is already 0, then that means it was set to zero by
  // the user.
  if (!local_dens_) {
    integral()->set_storage(storage_);
    return;
  }
  
  int nmem = i_offset(basis()->nbasis())*nm*sizeof(double);

  // if we're actually using local matrices, then there's no choice
  if (LocalSCMatrixKit::castdown(basis()->matrixkit().pointer())) {
    if (nmem > storage_)
      return;
  } else {
    if (nmem > storage_) {
      local_dens_=0;
      integral()->set_storage(storage_);
      return;
    }
  }

  integral()->set_storage(storage_-nmem);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "ETS")
