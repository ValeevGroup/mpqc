
#ifdef __GNUC__
#pragma implementation
#endif

#include <math/scmat/offset.h>
#include <math/optimize/diis.h>
#include <math/optimize/scextrapmat.h>

#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/basis/tbint.h>

#include <chemistry/qc/scf/scf.h>

///////////////////////////////////////////////////////////////////////////
// SCF

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
    printf("  total scf energy = %20.15f\n",eelec+nucrep);

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

///////////////////////////////////////////////////////////////////////////

void
SCF::compute_vector(double& eelec)
{
  int i;

  RefSelfConsistentExtrapolation extrap = new DIIS;
  
  init_vector();
  
  tbi = integral()->electron_repulsion();

  // initialize some junk
  double nucrep = molecule()->nuclear_repulsion_energy();

  for (int iter=0; iter < maxiter_; iter++) {
    // form the density from the current vector 
    double delta = new_density();
    
    // check convergence
    if (delta < desired_value_accuracy())
      break;

    // reset the density from time to time
    if (iter && !(iter%10))
      reset_density();
      
    // form the AO basis fock matrix
    ao_fock();

    // calculate the electronic energy
    eelec = scf_energy();
    printf("iter %5d energy = %20.15f delta = %15.10g\n",
           iter,eelec+nucrep,delta);

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
  tbi=0;
}

void
SCF::ao_gmat()
{
  double tnint=0;
  
  TwoBodyIntIter tbii(tbi);

  for (tbii.start(); tbii.ready(); tbii.next()) {
    ShellQuartetIter& q = tbii.current_quartet();

    int i=tbii.ishell();
    int j=tbii.jshell();
    int k=tbii.kshell();
    int l=tbii.lshell();

    int e12 = (i==j);
    int e34 = (k==l);
    int e13e24 = (i==k) && (j==l);
    int e_any = e12||e34||e13e24;
    
    for (q.start(); q.ready(); q.next()) {
      int ii=q.i();
      int jj=q.j();
      int kk=q.k();
      int ll=q.l();

      double pki_int = q.val();

      if (e_any) {
        int lij,lkl;
        double pkval;
                      
        if (jj == kk) {
          /*
           * if i=j=k or j=k=l, then this integral contributes
           * to J, K1, and K2 of G(ij), so
           * pkval = (ijkl) - 0.25 * ((ikjl)-(ilkj))
           *       = 0.5 * (ijkl)
           */
          if (ii == jj || kk == ll) {
            lij=i_offset(ii)+jj;
            lkl=i_offset(kk)+ll;
            pkval = (lij==lkl) ? 0.5*pki_int: pki_int;

            make_contribution(ii, jj, kk, ll, pkval, 5);

          } else {
            /*
             * if j=k, then this integral contributes
             * to J and K1 of G(ij)
             *
             * pkval = (ijkl) - 0.25 * (ikjl)
             *       = 0.75 * (ijkl)
             */
            lij=i_offset(ii)+jj;
            lkl=i_offset(kk)+ll;
            pkval = (lij==lkl) ? 0.5*pki_int: pki_int;
            
            make_contribution(ii, jj, kk, ll, pkval, 4);

            /*
             * this integral also contributes to K1 and K2 of
             * G(il)
             *
             * pkval = -0.25 * ((ijkl)+(ikjl))
             *       = -0.5 * (ijkl)
             */
            lij=i_offset(ii)+ll;
            lkl=ij_offset(kk,jj);
            pkval = (lij==lkl) ? 0.5*pki_int: pki_int;
            
            make_contribution(ii, ll, kk, jj, pkval, 3);
          }
        } else if (ii == kk || jj == ll) {
          /*
           * if i=k or j=l, then this integral contributes
           * to J and K2 of G(ij)
           *
           * pkval = (ijkl) - 0.25 * (ilkj)
           *       = 0.75 * (ijkl)
           */
          lij=i_offset(ii)+jj;
          lkl=i_offset(kk)+ll;
          pkval = (lij==lkl) ? 0.5*pki_int: pki_int;

          make_contribution(ii, jj, kk, ll, pkval, 4);

          /*
           * this integral also contributes to K1 and K2 of
           * G(ik)
           *
           * pkval = -0.25 * ((ijkl)+(ilkj))
           *       = -0.5 * (ijkl)
           */
          lij=i_offset(ii)+kk;
          lkl=ij_offset(jj,ll);
          pkval = (lij==lkl) ? 0.5*pki_int: pki_int;

          make_contribution(ii, kk, jj, ll, pkval, 3);

        } else {
          /*
           * This integral contributes to J of G(ij)
           *
           * pkval = (ijkl)
           */
          lij=i_offset(ii)+jj;
          lkl=i_offset(kk)+ll;
          pkval = (lij==lkl) ? 0.5*pki_int: pki_int;

          make_contribution(ii, jj, kk, ll, pkval, 1);

          /*
           * and to K1 of G(ik)
           *
           * pkval = -0.25 * (ijkl)
           */
          lij=i_offset(ii)+kk;
          lkl=ij_offset(jj,ll);
          pkval = (lij==lkl) ? 0.5*pki_int: pki_int;

          make_contribution(ii, kk, jj, ll, pkval, 2);

          if ((ii != jj) && (kk != ll)) {
            /*
             * if i!=j and k!=l, then this integral also
             * contributes to K2 of G(il)
             *
             * pkval = -0.25 * (ijkl)
             *
             * note: if we get here, then ik can't equal jl,
             * so pkval wasn't multiplied by 0.5 above.
             */
            lij=i_offset(ii)+ll;
            lkl=ij_offset(kk,jj);

            make_contribution(ii, ll, kk, jj, pkval, 2);
          }
        }
      } else {
        int lij,lkl;

        if (jj == kk) {
          /*
           * if j=k, then this integral contributes
           * to J and K1 of G(ij)
           *
           * pkval = (ijkl) - 0.25 * (ikjl)
           *       = 0.75 * (ijkl)
           */
          lij=i_offset(ii)+jj;
          lkl=i_offset(kk)+ll;

          make_contribution(ii, jj, kk, ll, pki_int, 4);

          /*
           * this integral also contributes to K1 and K2 of
           * G(il)
           *
           * pkval = -0.25 * ((ijkl)+(ikjl))
           *       = -0.5 * (ijkl)
           */
          lij=i_offset(ii)+ll;
          lkl=ij_offset(kk,jj);

          make_contribution(ii, ll, kk, jj, pki_int, 3);

        } else if (ii == kk || jj == ll) {
          /*
           * if i=k or j=l, then this integral contributes
           * to J and K2 of G(ij)
           *
           * pkval = (ijkl) - 0.25 * (ilkj)
           *       = 0.75 * (ijkl)
           */
          lij=i_offset(ii)+jj;
          lkl=i_offset(kk)+ll;

          make_contribution(ii, jj, kk, ll, pki_int, 4);

          /*
           * this integral also contributes to K1 and K2 of
           * G(ik)
           *
           * pkval = -0.25 * ((ijkl)+(ilkj))
           *       = -0.5 * (ijkl)
           */
          lij=i_offset(ii)+kk;
          lkl=ij_offset(jj,ll);

          make_contribution(ii, kk, jj, ll, pki_int, 3);

        } else {
          /*
           * This integral contributes to J of G(ij)
           *
           * pkval = (ijkl)
           */
          lij=i_offset(ii)+jj;
          lkl=i_offset(kk)+ll;

          make_contribution(ii, jj, kk, ll, pki_int, 1);

          /*
           * and to K1 of G(ik)
           *
           * pkval = -0.25 * (ijkl)
           */
          lij=i_offset(ii)+kk;
          lkl=ij_offset(jj,ll);
          
          make_contribution(ii, kk, jj, ll, pki_int, 2);

          /*
           * and to K2 of G(il)
           *
           * pkval = -0.25 * (ijkl)
           */
          lij=i_offset(ii)+ll;
          lkl=ij_offset(kk,jj);

          make_contribution(ii, ll, kk, jj, pki_int, 2);
        }
      }
    }

    tnint += (double) q.nint();
  }

  printf("%20.0f integrals\n",tnint);
}

//////////////////////////////////////////////////////////////////////////////

static void
set_scale(double& coulombscale, double& exchangescale,
          int i, int j, int k, int l)
{
  double scale = 1.0;

  if ((i!=k)||(j!=l))
    scale *= 2.0;

  if (i!=j)
    scale *= 2.0;

  coulombscale = 0.5*scale;
  exchangescale = -0.25*scale;

  if (k!=l)
    coulombscale *= 2.0;

  if ((k!=l)&&(i==j))
    exchangescale *= 2.0;
}

void
SCF::compute_gradient(const RefSCVector& gradient)
{
  init_gradient();
  gradient.assign(0.0);

  // handy things
  GaussianBasisSet& gbs = *basis().pointer();

  // form one electron contribution
  RefSymmSCMatrix lag = lagrangian();
  RefSymmSCMatrix dens = gradient_density();

  RefOneBodyDerivInt obints = integral()->deriv();
  for (int x=0; x < molecule()->natom(); x++) {
    for (int ish=0; ish < gbs.nshell(); ish++) {
      GaussianShell& gsi = gbs(ish);
      
      int istart = gbs.shell_to_function(ish);
      int iend = istart + gsi.nfunction();

      for (int jsh=0; jsh <= ish; jsh++) {
        GaussianShell& gsj = gbs(jsh);

        int jstart = gbs.shell_to_function(jsh);
        int jend = jstart + gsj.nfunction();

        obints->compute_overlap_shell(x,ish,jsh);

        const double *buf = obints->buffer();
        
        int index=0;
        double dx=0, dy=0, dz=0;
        for (int i=istart; i < iend; i++) {
          for (int j=jstart; j < jend; j++) {
            double ewij = lag.get_element(i,j);
            
            dx += buf[index++] * ewij;
            dy += buf[index++] * ewij;
            dz += buf[index++] * ewij;
          }
        }

        obints->compute_hcore_shell(x,ish,jsh);
        buf = obints->buffer();
        
        index=0;
        for (int i=istart; i < iend; i++) {
          for (int j=jstart; j < jend; j++) {
            double dij = dens.get_element(i,j);
            
            dx += buf[index++] * dij;
            dy += buf[index++] * dij;
            dz += buf[index++] * dij;
          }
        }

        if (ish != jsh) {
          dx *= 2.0;
          dy *= 2.0;
          dz *= 2.0;
        }

        gradient.accumulate_element(x*3+0, dx);
        gradient.accumulate_element(x*3+1, dy);
        gradient.accumulate_element(x*3+2, dz);
      }
    }
  }

  lag=0;
  dens=0;
  obints=0;
  
  gradient.print("one electron contribution");

  RefTwoBodyDerivInt tbints = integral()->electron_repulsion_deriv();
  
  // now the two electron part
  for (int i=0; i < gbs.nshell(); i++) {
    for (int j=0; j <= i; j++) {
      for (int k=0; k <= i; k++) {
        for (int l=0; l <= ((i==k)?j:k); l++) {
          tbints->compute_shell(i,j,k,l);
          const double * buf = tbints->buffer();
          
          double coulombscale, exchangescale;

          set_scale(coulombscale, exchangescale, i, j, k, l);
          make_gradient_contribution();
        }
      }
    }
  }
  
  // clean up
  tbints=0;
  
  done_gradient();
}

void
SCF::compute_hessian(const RefSymmSCMatrix& hessian)
{
  init_hessian();

  done_hessian();
}
