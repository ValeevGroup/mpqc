
#include <util/misc/timer.h>

#include <math/scmat/offset.h>
#include <math/optimize/diis.h>
#include <math/optimize/scextrapmat.h>

#include <chemistry/qc/basis/symmint.h>

#include <chemistry/qc/scf/scf.h>

///////////////////////////////////////////////////////////////////////////

void
SCF::compute_vector(double& eelec)
{
  int i;

  // one day this should be in the input
  RefSelfConsistentExtrapolation extrap = new DIIS;
  
  // set up subclass for vector calculation
  init_vector();
  
  // calculate the nuclear repulsion energy
  double nucrep = molecule()->nuclear_repulsion_energy();

  for (int iter=0; iter < maxiter_; iter++) {
    // form the density from the current vector 
    double delta = new_density();
    
    // check convergence
    if (delta < desired_value_accuracy())
      break;

    // reset the density from time to time
    if (iter && !(iter%dens_reset_freq_))
      reset_density();
      
    // form the AO basis fock matrix
    ao_fock();

    // calculate the electronic energy
    eelec = scf_energy();
    printf("  iter %5d energy = %20.15f delta = %15.10g\n",
           iter+1,eelec+nucrep,delta);

    // now extrapolate the fock matrix
    RefSCExtrapData data = extrap_data();
    RefSCExtrapError error = extrap_error();
    extrap->extrapolate(data,error);
    data=0;
    error=0;

    // diagonalize effective MO fock to get MO vector
    RefSCMatrix nvector = scf_vector_.clone();
    RefDiagSCMatrix evals(basis_dimension(), basis_matrixkit());
  
    RefSymmSCMatrix eff = effective_fock();
    eff.diagonalize(evals,nvector);
    set_occupations(evals);

    eff=0;
    evals=0;
    
    // transform MO vector to AO basis
    scf_vector_ = scf_vector_ * nvector;
    nvector=0;
    
    // and orthogonalize vector
    scf_vector_->schmidt_orthog(overlap().pointer(),basis()->nbasis());
  }
      
  eigenvectors_ = scf_vector_;
  eigenvectors_.computed() = 1;
  
  // now clean up
  done_vector();

  extrap = 0;
}
