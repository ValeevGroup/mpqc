
#ifdef __GNUC__
#pragma implementation
#endif

#include <iostream.h>
#include <iomanip.h>
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
  s.get(int_store_);
  s.get(dens_reset_freq_);
  s.get(reset_occ_);
  s.get(level_shift_);

  extrap_.restore_state(s);

  integral()->set_storage(int_store_);
  scf_grp_ = basis()->matrixkit()->messagegrp();
}

SCF::SCF(const RefKeyVal& keyval) :
  OneBodyWavefunction(keyval),
  maxiter_(40),
  int_store_(0),
  level_shift_(0),
  reset_occ_(0),
  dens_reset_freq_(10)
{
  if (keyval->exists("maxiter"))
    maxiter_ = keyval->intvalue("maxiter");

  if (keyval->exists("density_reset_frequency"))
    dens_reset_freq_ = keyval->intvalue("density_reset_frequency");

  if (keyval->exists("integral_storage"))
    int_store_ = keyval->intvalue("integral_storage");

  if (keyval->exists("reset_occupations"))
    reset_occ_ = keyval->booleanvalue("reset_occupations");

  if (keyval->exists("level_shift"))
    level_shift_ = keyval->doublevalue("level_shift");

  extrap_ = keyval->describedclassvalue("extrap");
  if (extrap_.null())
    extrap_ = new DIIS;
  
  // first see if guess_wavefunction is a wavefunction, then check to
  // see if it's a string.
  if (keyval->exists("guess_wavefunction")) {
    guess_wfn_ = keyval->describedclassvalue("guess_wavefunction");
    if (guess_wfn_.null()) {
      char *path = keyval->pcharvalue("guess_wavefunction");
      struct stat sb;
      if (path && stat(path, &sb)==0 && sb.st_size) {
        StateInBinXDR s(path);
        guess_wfn_.restore_state(s);
        delete[] path;
      }
    }
  }
  
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
  s.put(reset_occ_);
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
    if (me==0) {
      cout.unsetf(ios::fixed);
      cout << endl << indent << "SCF::compute: energy accuracy = "
           << setw(10) << setprecision(7) << desired_value_accuracy()
           << endl << endl;
    }

    double eelec;
    compute_vector(eelec);
      
    // this will be done elsewhere eventually
    double nucrep = molecule()->nuclear_repulsion_energy();
    if (me==0) {
      cout.setf(ios::fixed);
      cout << endl << indent << "total scf energy = " <<
        setw(20) << setprecision(15) << eelec+nucrep << endl;
    }

    set_energy(eelec+nucrep);
    set_actual_value_accuracy(desired_value_accuracy());
  }

  if (gradient_needed()) {
    RefSCVector gradient = matrixkit()->vector(moldim());

    if (me==0) {
      cout.unsetf(ios::fixed);
      cout << endl << indent << "SCF::compute: gradient accuracy = "
           << setw(10) << setprecision(7) << desired_gradient_accuracy()
           << endl << endl;
    }

    compute_gradient(gradient);
    gradient.print("cartesian gradient");
    set_gradient(gradient);

    set_actual_gradient_accuracy(desired_gradient_accuracy());
  }
  
  if (hessian_needed()) {
    RefSymmSCMatrix hessian = matrixkit()->symmmatrix(moldim());
    
    if (me==0) {
      cout.unsetf(ios::fixed);
      cout << endl << indent << "SCF::compute: hessian accuracy = "
           << setw(10) << setprecision(7) << desired_hessian_accuracy()
           << endl << endl;
    }

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
  
  if (!LocalSymmSCMatrix::castdown(l)) {
    RefSCMatrixKit k = new LocalSCMatrixKit;
    l = k->symmmatrix(m.dim());
    l->convert(m);

    if (access == Accum)
      l->assign(0.0);
  } else if (scf_grp_->n() > 1 && access==Accum) {
    l = m.clone();
    l.assign(0.0);
  }

  p = LocalSymmSCMatrix::castdown(l)->get_data();
  return l;
}

void
SCF::initial_vector()
{
  int me=scf_grp_->me();
  
  // if guess_wfn_ is non-null then try to get a guess vector from it.
  // First check that the same basis is used...if not, then project the
  // guess vector into the present basis.
  // right now the check is crude...there should be an equiv member in
  // GaussianBasisSet
  if (guess_wfn_.nonnull()) {
    if (guess_wfn_->basis()->nbasis() == basis()->nbasis()) {
      if (me==0) {
        cout << indent
             << "Using guess_wavefunction as starting vector\n";
      }
      cout << incindent << incindent;
      eigenvectors_ = guess_wfn_->eigenvectors();
      cout << decindent << decindent;
    } else {
      if (me==0) {
        cout << indent
             << "Projecting guess_wavefunction into the present basis set\n";
      }
      cout << incindent << incindent;
      eigenvectors_ = projected_eigenvectors(guess_wfn_);
      cout << decindent << decindent;
    }

    // we should only have to do this once, so free up memory used
    // for the old wavefunction
    guess_wfn_=0;

    if (me==0)
      cout << endl;
      
  } else {
    if (me==0) {
      cout << indent
           << "Starting from core Hamiltonian guess\n\n";
    }
    eigenvectors_ = hcore_guess();
  }
}
