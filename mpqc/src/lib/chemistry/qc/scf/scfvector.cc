
#include <iostream.h>
#include <iomanip.h>

#include <util/misc/timer.h>
#include <util/misc/formio.h>

#include <math/scmat/offset.h>
#include <math/scmat/blocked.h>

#include <math/optimize/diis.h>
#include <math/optimize/scextrapmat.h>

#include <chemistry/qc/basis/symmint.h>

#include <chemistry/qc/scf/scf.h>
#include <chemistry/qc/scf/scfden.h>
#include <chemistry/qc/scf/scflocal.h>

///////////////////////////////////////////////////////////////////////////

void
SCF::compute_vector(double& eelec)
{
  tim_enter("vector");
  int i;

  // one day this should be in the input
  RefSelfConsistentExtrapolation extrap = new DIIS;
  
  // create level shifter
  LevelShift *level_shift = new LevelShift(this);
  level_shift->reference();
  
  // set up subclass for vector calculation
  init_vector();
  
  // calculate the nuclear repulsion energy
  double nucrep = molecule()->nuclear_repulsion_energy();

  for (int iter=0; iter < maxiter_; iter++) {
    // form the density from the current vector 
    tim_enter("density");
    double delta = new_density();
    tim_exit("density");
    
    // check convergence
    if (delta < desired_value_accuracy())
      break;

    // reset the density from time to time
    if (iter && !(iter%dens_reset_freq_))
      reset_density();
      
    // form the AO basis fock matrix
    tim_enter("fock");
    ao_fock();
    tim_exit("fock");

    // calculate the electronic energy
    eelec = scf_energy();
    if (scf_grp_->me()==0)
      cout << indent << "iter " << setw(5) << iter+1 <<
        " energy = " << setw(20) << setprecision(15) << eelec+nucrep <<
        " delta = " << setw(20) << setprecision(15) << delta << endl;

    // now extrapolate the fock matrix
    tim_enter("extrap");
    RefSCExtrapData data = extrap_data();
    RefSCExtrapError error = extrap_error();
    extrap->extrapolate(data,error);
    data=0;
    error=0;
    tim_exit("extrap");

    // diagonalize effective MO fock to get MO vector
    tim_enter("evals");
    RefSCMatrix nvector = scf_vector_.clone();
    RefDiagSCMatrix evals(basis_dimension(), basis_matrixkit());
  
    RefSymmSCMatrix eff = effective_fock();

    // level shift effective fock
    level_shift->set_shift(level_shift_);
    eff.element_op(level_shift);
    
    eff.diagonalize(evals,nvector);
    tim_exit("evals");

    // now un-level shift eigenvalues
    level_shift->set_shift(-level_shift_);
    evals.element_op(level_shift);
    
    set_occupations(evals);

    eff=0;
    evals=0;
    
    // transform MO vector to AO basis
    scf_vector_ = scf_vector_ * nvector;
    nvector=0;
    
    // and orthogonalize vector
    tim_enter("schmidt");
    scf_vector_->schmidt_orthog(overlap().pointer(),basis()->nbasis());
    tim_exit("schmidt");
  }
      
  eigenvectors_ = scf_vector_;
  eigenvectors_.computed() = 1;
  
  // now clean up
  done_vector();

  extrap = 0;
  tim_exit("vector");
  //tim_print(0);
}

////////////////////////////////////////////////////////////////////////////

class ExtrapErrorOp : public BlockedSCElementOp {
  private:
    SCF *scf_;

  public:
    ExtrapErrorOp(SCF *s) : scf_(s) {}
    ~ExtrapErrorOp() {}

    int has_side_effects() { return 1; }

    void process(SCMatrixBlockIter& bi) {
      int ir=current_block();
      
      for (bi.reset(); bi; bi++) {
        int i=bi.i();
        int j=bi.j();
        if (scf_->occupation(ir,i) == scf_->occupation(ir,j))
          bi.set(0.0);
      }
    }
};

RefSCExtrapError
SCF::extrap_error()
{
  RefSymmSCMatrix mofock = effective_fock();
  
  RefSCElementOp op = new ExtrapErrorOp(this);
  mofock.element_op(op);
  
  RefSymmSCMatrix aoerror = mofock.clone();
  aoerror.assign(0.0);
  aoerror.accumulate_transform(scf_vector_,mofock);
  mofock=0;

  RefSCExtrapError error = new SymmSCMatrixSCExtrapError(aoerror);
  aoerror=0;

  return error;
}
