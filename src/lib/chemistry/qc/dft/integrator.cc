//
// integrator.cc
//
// Copyright (C) 1997 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#include <util/misc/math.h>

#include <util/misc/regtime.h>
#include <util/misc/formio.h>
#include <util/state/stateio.h>
#include <util/container/carray.h>
#include <chemistry/qc/dft/integrator.h>

using namespace std;
using namespace sc;

namespace sc {

#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
template void delete_c_array2<Ref<RadialIntegrator> >(Ref<RadialIntegrator>**);
template Ref<RadialIntegrator>**
  new_c_array2<Ref<RadialIntegrator> >(int,int,Ref<RadialIntegrator>);
template void delete_c_array3<Ref<RadialIntegrator> >(Ref<RadialIntegrator>***);
template Ref<RadialIntegrator>***
  new_c_array3<Ref<RadialIntegrator> >(int,int,int,Ref<RadialIntegrator>);

template void delete_c_array2<Ref<AngularIntegrator> >(Ref<AngularIntegrator>**);
template Ref<AngularIntegrator>**
  new_c_array2<Ref<AngularIntegrator> >(int,int,Ref<AngularIntegrator>);
template void delete_c_array3<Ref<AngularIntegrator> >(Ref<AngularIntegrator>***);
template Ref<AngularIntegrator>***
  new_c_array3<Ref<AngularIntegrator> >(int,int,int,Ref<AngularIntegrator>);
#endif

//#define CHECK_ALIGN(v) if(int(&v)&7)ExEnv::outn()<<"Bad Alignment: "<< ## v <<endl;
#define CHECK_ALIGN(v)

///////////////////////////////////////////////////////////////////////////
// utility functions

inline static double
dot(double v[3], double w[3])
{
  return v[0]*w[0] + v[1]*w[1] + v[2]*w[2];
}

inline static double
norm(double v[3])
{
  return sqrt(dot(v,v));
}

static double
get_radius(const Ref<Molecule> &mol, int iatom)
{
  double r = mol->atominfo()->maxprob_radius(mol->Z(iatom));
  if (r == 0) {
      static bool warned = false;
      if (!warned) {
          ExEnv::out0()
              << indent << "WARNING: BeckeIntegrationWeight usually uses the atomic maximum" << std::endl
              << indent << "         probability radius, however this is not available for" << std::endl
              << indent << "         one of the atoms in your system. The Bragg radius will" << std::endl
              << indent << "         be used instead." << std::endl;
          warned = true;
        }
      r = mol->atominfo()->bragg_radius(mol->Z(iatom));
    }
  if (r == 0) {
      ExEnv::out0()
          << indent << "ERROR: BeckeIntegrationWeight could not find a maximum probability" << std::endl
          << indent << "       or a Bragg radius for an atom" << std::endl;
      abort();
    }
  return r;
}

///////////////////////////////////////////////////////////////////////////
// DenIntegratorThread

//ThreadLock *tlock;

class DenIntegratorThread: public Thread {
  protected:
    // data common to all threads
    int nthread_;
    int nshell_;
    int nbasis_;
    int natom_;
    int n_integration_center_;
    DenIntegrator *integrator_;
    int spin_polarized_;
    int need_hessian_;
    int need_gradient_;
    GaussianBasisSet *basis_;
    bool linear_scaling_;
    bool use_dmat_bound_;
    double value_;
    DenFunctional *func_;

    double accuracy_;
    int compute_potential_integrals_;

    // data local to thread
    int ithread_;
    Ref<BatchElectronDensity> den_;
    double *alpha_vmat_; // lower triangle of xi_i(r) v(r) xi_j(r) integrals
    double *beta_vmat_; // lower triangle of xi_i(r) v(r) xi_j(r) integrals
    double *nuclear_gradient_;
    double *w_gradient_;
    double *f_gradient_;

  public:
    DenIntegratorThread(int ithread, int nthread,
                        DenIntegrator *integrator,
                        DenFunctional *func,
                        const Ref<BatchElectronDensity> &den,
                        int linear_scaling,
                        int use_dmat_bound,
                        double accuracy,
                        int compute_potential_integrals,
                        int need_nuclear_gradient);
    virtual ~DenIntegratorThread();
    double do_point(int iatom, const SCVector3 &r,
                    double weight, double multiplier,
                    double *nuclear_gradient,
                    double *f_gradient, double *w_gradient);
    double *nuclear_gradient() { return nuclear_gradient_; }
    double *alpha_vmat() { return alpha_vmat_; }
    double *beta_vmat() { return beta_vmat_; }
    double value() { return value_; }
};

DenIntegratorThread::DenIntegratorThread(int ithread, int nthread,
                                         DenIntegrator *integrator,
                                         DenFunctional *func,
                                         const Ref<BatchElectronDensity> &den,
                                         int linear_scaling,
                                         int use_dmat_bound,
                                         double accuracy,
                                         int compute_potential_integrals,
                                         int need_nuclear_gradient)
{
  value_ = 0.0;

  den_ = den;
  ithread_ = ithread;
  nthread_ = nthread;
  integrator_ = integrator;
  spin_polarized_ = integrator->spin_polarized();
  nshell_ = integrator->basis()->nshell();
  nbasis_ = integrator->basis()->nbasis();
  natom_ = integrator->basis()->molecule()->natom();
  n_integration_center_
      = integrator->basis()->molecule()->n_non_q_atom();
  need_gradient_ = func->need_density_gradient();
  need_hessian_ = func->need_density_hessian();
  basis_ = integrator->basis().pointer();
  linear_scaling_ = linear_scaling;
  use_dmat_bound_ = use_dmat_bound;
  func_ = func;
  accuracy_ = accuracy;
  compute_potential_integrals_ = compute_potential_integrals;
  den_->set_accuracy(accuracy);

  if (need_nuclear_gradient) {
      nuclear_gradient_ = new double[3*natom_];
      memset(nuclear_gradient_, 0, 3*natom_*sizeof(double));
    }
  else nuclear_gradient_ = 0;

  alpha_vmat_ = 0;
  beta_vmat_ = 0;
  if (compute_potential_integrals_) {
      int ntri = (nbasis_*(nbasis_+1))/2;
      alpha_vmat_ = new double[ntri];
      memset(alpha_vmat_, 0, sizeof(double)*ntri);
      if (spin_polarized_) {
          beta_vmat_ = new double[ntri];
          memset(beta_vmat_, 0, sizeof(double)*ntri);
        }
    }

  w_gradient_ = 0;
  f_gradient_ = 0;
  if (nuclear_gradient_) {
      w_gradient_ = new double[n_integration_center_*3];
      f_gradient_ = new double[natom_*3];
    }
}

DenIntegratorThread::~DenIntegratorThread()
{
  delete[] alpha_vmat_;
  delete[] beta_vmat_;
  delete[] nuclear_gradient_;
  delete[] f_gradient_;
  delete[] w_gradient_;
}

///////////////////////////////////////////////////////////////////////////
// DenIntegrator

static ClassDesc DenIntegrator_cd(
  typeid(DenIntegrator),"DenIntegrator",1,"public SavableState",
  0, 0, 0);

DenIntegrator::DenIntegrator(StateIn& s):
  SavableState(s)
{
  init_object();
  s.get(linear_scaling_);
  s.get(use_dmat_bound_);
}

DenIntegrator::DenIntegrator()
{
  init_object();
}

DenIntegrator::DenIntegrator(const Ref<KeyVal>& keyval)
{
  init_object();

  linear_scaling_ = keyval->booleanvalue("linear_scaling",
                                         KeyValValueboolean(linear_scaling_));
  use_dmat_bound_ = keyval->booleanvalue("use_dmat_bound",
                                         KeyValValueboolean(use_dmat_bound_));
}

DenIntegrator::~DenIntegrator()
{
}

void
DenIntegrator::save_data_state(StateOut& s)
{
  s.put(linear_scaling_);
  s.put(use_dmat_bound_);
}

void
DenIntegrator::init_object()
{
  threadgrp_ = ThreadGrp::get_default_threadgrp();
  messagegrp_ = MessageGrp::get_default_messagegrp();
  compute_potential_integrals_ = 0;
  accuracy_ = DBL_EPSILON;
  linear_scaling_ = 1;
  use_dmat_bound_ = 1;
  alpha_vmat_ = 0;
  beta_vmat_ = 0;
}

void
DenIntegrator::set_compute_potential_integrals(int i)
{
  compute_potential_integrals_=i;
}

void
DenIntegrator::init(const Ref<Wavefunction> &wfn)
{
  init(wfn->basis(), wfn->integral());
}

void
DenIntegrator::init(const Ref<GaussianBasisSet> &basis,
                    const Ref<Integral> &integral)
{
  basis_ = basis;
  integral_ = integral;
  den_ = new BatchElectronDensity(basis_, integral_,
                                  accuracy_);
  den_->set_linear_scaling(linear_scaling_);
  den_->set_use_dmat_bound(use_dmat_bound_);
  den_->init();
}

void
DenIntegrator::set_accuracy(double a)
{
  accuracy_ = a;
  if (den_) den_->set_accuracy(a);
}

void
DenIntegrator::done()
{
  basis_ = 0;
  integral_ = 0;
  den_ = 0;
}

void
DenIntegrator::init_integration(const Ref<DenFunctional> &func,
                                const RefSymmSCMatrix& densa,
                                const RefSymmSCMatrix& densb,
                                double *nuclear_gradient)
{
  int i;
  value_ = 0.0;

  den_->set_densities(densa,densb);
  spin_polarized_ = den_->spin_polarized();

  func->set_compute_potential(
      compute_potential_integrals_ || nuclear_gradient != 0);

  func->set_spin_polarized(spin_polarized_);

  natom_ = basis_->molecule()->natom();
  n_integration_center_
      = basis_->molecule()->n_non_q_atom();

  nshell_ = basis_->nshell();
  nbasis_ = basis_->nbasis();

  delete[] alpha_vmat_;
  delete[] beta_vmat_;
  alpha_vmat_ = 0;
  beta_vmat_ = 0;
  if (compute_potential_integrals_) {
      int ntri = (nbasis_*(nbasis_+1))/2;
      alpha_vmat_ = new double[ntri];
      memset(alpha_vmat_, 0, sizeof(double)*ntri);
      if (spin_polarized_) {
          beta_vmat_ = new double[ntri];
          memset(beta_vmat_, 0, sizeof(double)*ntri);
        }
    }
}

void
DenIntegrator::done_integration()
{
  messagegrp_->sum(value_);
  if (compute_potential_integrals_) {
      int ntri = (nbasis_*(nbasis_+1))/2;
      messagegrp_->sum(alpha_vmat_,ntri);
      if (spin_polarized_) {
          messagegrp_->sum(beta_vmat_,ntri);
        }
    }
}

double
DenIntegratorThread::do_point(int iatom, const SCVector3 &r,
                        double weight, double multiplier,
                        double *nuclear_gradient,
                        double *f_gradient, double *w_gradient)
{
  int i,j;
  double w_mult = weight * multiplier;

  CHECK_ALIGN(w_mult);

  PointInputData id(r);

  den_->compute_density(r,
                        &id.a.rho,
                        (need_gradient_?id.a.del_rho:0),
                        (need_hessian_?id.a.hes_rho:0),
                        &id.b.rho,
                        (need_gradient_?id.b.del_rho:0),
                        (need_hessian_?id.b.hes_rho:0));

  id.compute_derived(spin_polarized_, need_gradient_, need_hessian_);

  int ncontrib = den_->ncontrib();
  int *contrib = den_->contrib();
  int ncontrib_bf = den_->ncontrib_bf();
  int *contrib_bf = den_->contrib_bf();
  double *bs_values = den_->bs_values();
  double *bsg_values = den_->bsg_values();
  double *bsh_values = den_->bsh_values();

  PointOutputData od;
  if ( (id.a.rho + id.b.rho) > 1e2*DBL_EPSILON) {
      if (nuclear_gradient == 0) {
          func_->point(id, od);
        }
      else {
          func_->gradient(id, od, f_gradient, iatom, basis_,
                          den_->alpha_density_matrix(),
                          den_->beta_density_matrix(),
                          ncontrib, contrib,
                          ncontrib_bf, contrib_bf,
                          bs_values, bsg_values,
                          bsh_values);
        }
    }
  else {
      return id.a.rho + id.b.rho;
    }
  
  value_ += od.energy * w_mult;

  if (compute_potential_integrals_) {
      // the contribution to the potential integrals
      if (need_gradient_) {
          double gradsa[3], gradsb[3];
          gradsa[0] = w_mult*(2.0*od.df_dgamma_aa*id.a.del_rho[0] +
                                  od.df_dgamma_ab*id.b.del_rho[0]);
          gradsa[1] = w_mult*(2.0*od.df_dgamma_aa*id.a.del_rho[1] +
                                 od.df_dgamma_ab*id.b.del_rho[1]);
          gradsa[2] = w_mult*(2.0*od.df_dgamma_aa*id.a.del_rho[2] +
                                  od.df_dgamma_ab*id.b.del_rho[2]);
          double drhoa = w_mult*od.df_drho_a, drhob=0.0;
          if (spin_polarized_) {
              drhob = w_mult*od.df_drho_b;
              gradsb[0] = w_mult*(2.0*od.df_dgamma_bb*id.b.del_rho[0] +
                                      od.df_dgamma_ab*id.a.del_rho[0]);
              gradsb[1] = w_mult*(2.0*od.df_dgamma_bb*id.b.del_rho[1] +
                                      od.df_dgamma_ab*id.a.del_rho[1]);
              gradsb[2] = w_mult*(2.0*od.df_dgamma_bb*id.b.del_rho[2] +
                                      od.df_dgamma_ab*id.a.del_rho[2]);
            }

          for (int j=0; j<ncontrib_bf; j++) {
              int jt = contrib_bf[j];
              double dfdra_phi_m = drhoa*bs_values[j];
              double dfdga_phi_m = gradsa[0]*bsg_values[j*3+0] +
                                   gradsa[1]*bsg_values[j*3+1] +
                                   gradsa[2]*bsg_values[j*3+2];
              double vamu = dfdra_phi_m + dfdga_phi_m, vbmu=0.0;
              double dfdrb_phi_m, dfdgb_phi_m;
              if (spin_polarized_) {
                  dfdrb_phi_m = drhob*bs_values[j];
                  dfdgb_phi_m = gradsb[0]*bsg_values[j*3+0] +
                                       gradsb[1]*bsg_values[j*3+1] +
                                       gradsb[2]*bsg_values[j*3+2];
                  vbmu = dfdrb_phi_m + dfdgb_phi_m;
                }

              int jtoff = (jt*(jt+1))>>1;

              for (int k=0; k <= j; k++) {
                  int kt = contrib_bf[k];
                  int jtkt = jtoff + kt;

                  double dfdga_phi_n = gradsa[0]*bsg_values[k*3+0] +
                                       gradsa[1]*bsg_values[k*3+1] +
                                       gradsa[2]*bsg_values[k*3+2];
                  alpha_vmat_[jtkt] += vamu * bs_values[k] +
                                     dfdga_phi_n * bs_values[j];
                  if (spin_polarized_) {
                      double dfdgb_phi_n = gradsb[0]*bsg_values[k*3+0] +
                                           gradsb[1]*bsg_values[k*3+1] +
                                           gradsb[2]*bsg_values[k*3+2];
                      beta_vmat_[jtkt] += vbmu * bs_values[k] +
                                          dfdgb_phi_n * bs_values[j];
                    }
                }
            }
        }
      else {
          double drhoa = w_mult*od.df_drho_a;
          double drhob = w_mult*od.df_drho_b;
          for (int j=0; j<ncontrib_bf; j++) {
              int jt = contrib_bf[j];
              double bsj = bs_values[j];
              double dfa_phi_m = drhoa * bsj;
              double dfb_phi_m = drhob * bsj;
              int jtoff = (jt*(jt+1))>>1;
              for (int k=0; k <= j; k++) {
                  int kt = contrib_bf[k];
                  int jtkt = jtoff + kt;
                  double bsk = bs_values[k];
                  alpha_vmat_[jtkt] += dfa_phi_m * bsk;
                  if (spin_polarized_)
                      beta_vmat_[jtkt] += dfb_phi_m * bsk;
                }
            }
        }
    }

  if (nuclear_gradient != 0) {
      // the contribution from f dw/dx
      if (w_gradient) {
          for (int icenter = 0; icenter<n_integration_center_; icenter++) {
              int iatom = basis_->molecule()->non_q_atom(icenter);
              for (int ixyz=0; ixyz<3; ixyz++) {
                  nuclear_gradient[iatom*3+ixyz]
                      += w_gradient[icenter*3+ixyz] * od.energy * multiplier;
                }
            }
        }
      // the contribution from (df/dx) w
      for (i=0; i<natom_*3; i++) {
          nuclear_gradient[i] += f_gradient[i] * w_mult;
        }
    }

  return id.a.rho + id.b.rho;
}

///////////////////////////////////////////////////////////////////////////
// IntegrationWeight

static ClassDesc IntegrationWeight_cd(
  typeid(IntegrationWeight),"IntegrationWeight",1,"public SavableState",
  0, 0, 0);

IntegrationWeight::IntegrationWeight(StateIn& s):
  SavableState(s)
{
}

IntegrationWeight::IntegrationWeight()
{
}

IntegrationWeight::IntegrationWeight(const Ref<KeyVal>& keyval)
{
}

IntegrationWeight::~IntegrationWeight()
{
}

void
IntegrationWeight::save_data_state(StateOut& s)
{
}

void
IntegrationWeight::init(const Ref<Molecule> &mol, double tolerance)
{
  mol_ = mol;
  tol_ = tolerance;
}

void
IntegrationWeight::done()
{
}

void
IntegrationWeight::fd_w(int icenter, SCVector3 &point,
                        double *fd_grad_w)
{
  if (!fd_grad_w) return;
  double delta = 0.001;
  int natom = mol_->natom();
  Ref<Molecule> molsav = mol_;
  Ref<Molecule> dmol = new Molecule(*mol_.pointer());
  for (int i=0; i<natom; i++) {
      for (int j=0; j<3; j++) {
          dmol->r(i,j) += delta;
          if (icenter == i) point[j] += delta;
          init(dmol,tol_);
          double w_plus = w(icenter, point);
          dmol->r(i,j) -= 2*delta;
          if (icenter == i) point[j] -= 2*delta;
          init(dmol,tol_);
          double w_minus = w(icenter, point);
          dmol->r(i,j) += delta;
          if (icenter == i) point[j] += delta;
          fd_grad_w[i*3+j] = (w_plus-w_minus)/(2.0*delta);
//            ExEnv::outn() << scprintf("%d,%d %12.10f %12.10f %12.10f",
//                             i,j,w_plus,w_minus,fd_grad_w[i*3+j])
//                 << endl;
        }
    }
  init(molsav, tol_);
}

void
IntegrationWeight::test(int icenter, SCVector3 &point)
{
  int natom = mol_->natom();
  int natom3 = natom*3;

  // tests over sums of weights 
  int i;
  double sum_weight = 0.0;
  for (i=0; i<natom; i++) {
      double weight = w(i,point);
      sum_weight += weight;
    }
  if (fabs(1.0 - sum_weight) > DBL_EPSILON) {
      ExEnv::out0() << "IntegrationWeight::test: failed on weight" << endl;
      ExEnv::out0() << "sum_w = " << sum_weight << endl;
    }

  // finite displacement tests of weight gradients
  double *fd_grad_w = new double[natom3];
  double *an_grad_w = new double[natom3];
  w(icenter, point, an_grad_w);
  fd_w(icenter, point, fd_grad_w);
  for (i=0; i<natom3; i++) {
      double mag = fabs(fd_grad_w[i]);
      double err = fabs(fd_grad_w[i]-an_grad_w[i]);
      int bad = 0;
      if (mag > 0.00001 && err/mag > 0.01) bad = 1;
      else if (err > 0.00001) bad = 1;
      if (bad) {
          ExEnv::out0() << "iatom = " << i/3
               << " ixyx = " << i%3
               << " icenter = " << icenter << " point = " << point << endl;
          ExEnv::out0() << scprintf("dw/dx bad: fd_val=%16.13f an_val=%16.13f err=%16.13f",
                           fd_grad_w[i], an_grad_w[i],
                           fd_grad_w[i]-an_grad_w[i])
               << endl;
        }
    }
  delete[] fd_grad_w;
  delete[] an_grad_w;  
}

void
IntegrationWeight::test()
{
  SCVector3 point;
  for (int icenter=0; icenter<mol_->natom(); icenter++) {
      for (point[0]=-1; point[0]<=1; point[0]++) {
          for (point[1]=-1; point[1]<=1; point[1]++) {
              for (point[2]=-1; point[2]<=1; point[2]++) {
                  test(icenter, point);
                }
            }
        }
    }
}

///////////////////////////////////////////////////////////////////////////
// BeckeIntegrationWeight

// utility functions

inline static double
calc_s(double m)
{
  double m1 = 1.5*m - 0.5*m*m*m;
  double m2 = 1.5*m1 - 0.5*m1*m1*m1;
  double m3 = 1.5*m2 - 0.5*m2*m2*m2;
  return 0.5*(1.0-m3);
}

inline static double
calc_f3_prime(double m)
{
  double m1 = 1.5*m - 0.5*m*m*m;
  double m2 = 1.5*m1 - 0.5*m1*m1*m1;
  double m3 = 1.5 *(1.0 - m2*m2);
  double n2 = 1.5 *(1.0 - m1*m1);
  double o1 = 1.5 *(1.0 - m*m);
  return m3*n2*o1;
}

static ClassDesc BeckeIntegrationWeight_cd(
  typeid(BeckeIntegrationWeight),"BeckeIntegrationWeight",1,"public IntegrationWeight",
  0, create<BeckeIntegrationWeight>, create<BeckeIntegrationWeight>);

BeckeIntegrationWeight::BeckeIntegrationWeight(StateIn& s):
  SavableState(s),
  IntegrationWeight(s)
{
  n_integration_centers = 0;
  atomic_radius = 0;
  a_mat = 0;
  oorab = 0;
  centers = 0;
}

BeckeIntegrationWeight::BeckeIntegrationWeight()
{
  n_integration_centers = 0;
  centers = 0;
  atomic_radius = 0;
  a_mat = 0;
  oorab = 0;
}

BeckeIntegrationWeight::BeckeIntegrationWeight(const Ref<KeyVal>& keyval):
  IntegrationWeight(keyval)
{
  n_integration_centers = 0;
  centers = 0;
  atomic_radius = 0;
  a_mat = 0;
  oorab = 0;
}

BeckeIntegrationWeight::~BeckeIntegrationWeight()
{
  done();
}

void
BeckeIntegrationWeight::save_data_state(StateOut& s)
{
  IntegrationWeight::save_data_state(s);
}

void
BeckeIntegrationWeight::init(const Ref<Molecule> &mol, double tolerance)
{
  done();
  IntegrationWeight::init(mol, tolerance);

  // We only want to include to include "atoms" that correspond to nuclei,
  // not charges.
  n_integration_centers = mol->n_non_q_atom();

  vector<double> atomic_radius(n_integration_centers);
  centers = new SCVector3[n_integration_centers];

  for (int icenter=0; icenter<n_integration_centers; icenter++) {
      int iatom = mol->non_q_atom(icenter);

      atomic_radius[icenter] = get_radius(mol, iatom);

      centers[icenter].x() = mol->r(iatom,0);
      centers[icenter].y() = mol->r(iatom,1);
      centers[icenter].z() = mol->r(iatom,2);
    }
      
  a_mat = new double*[n_integration_centers];
  a_mat[0] = new double[n_integration_centers*n_integration_centers];
  oorab = new double*[n_integration_centers];
  oorab[0] = new double[n_integration_centers*n_integration_centers];

  for (int icenter=0; icenter<n_integration_centers; icenter++) {
      if (icenter) {
          a_mat[icenter] = &a_mat[icenter-1][n_integration_centers];
          oorab[icenter] = &oorab[icenter-1][n_integration_centers];
        }

      double atomic_radius_a = atomic_radius[icenter];
      for (int jcenter=0; jcenter < n_integration_centers; jcenter++) {
          double chi=atomic_radius_a/atomic_radius[jcenter];
          double uab=(chi-1.)/(chi+1.);
          a_mat[icenter][jcenter] = uab/(uab*uab-1.);
          if (icenter!=jcenter) {
              oorab[icenter][jcenter]
                  = 1./centers[icenter].dist(centers[jcenter]);
            }
          else {
              oorab[icenter][jcenter] = 0.0;
            }
        }
    }

}

void
BeckeIntegrationWeight::done()
{
  delete[] atomic_radius;
  atomic_radius = 0;

  delete[] centers;
  centers = 0;

  if (a_mat) {
      delete[] a_mat[0];
      delete[] a_mat;
      a_mat = 0;
    }

  if (oorab) {
      delete[] oorab[0];
      delete[] oorab;
      oorab = 0;
    }

  n_integration_centers = 0;
}

double
BeckeIntegrationWeight::compute_p(int icenter, SCVector3&point)
{
  double ra = point.dist(centers[icenter]);
  double *ooraba = oorab[icenter];
  double *aa = a_mat[icenter];

  double p = 1.0;
  for (int jcenter=0; jcenter < n_integration_centers; jcenter++) {
      if (icenter != jcenter) {
          double mu = (ra-point.dist(centers[jcenter]))*ooraba[jcenter];

          if (mu <= -1.)
              continue; // s(-1) == 1.0
          else if (mu >= 1.) {
              return 0.0; // s(1) == 0.0
            }
          else
              p *= calc_s(mu + aa[jcenter]*(1.-mu*mu));
        }
    }

  return p;
}

// compute derivative of mu(grad_center,bcenter) wrt grad_center;
// NB: the derivative is independent of the (implicit) wcenter
// provided that wcenter!=grad_center
void
BeckeIntegrationWeight::compute_grad_nu(int grad_center, int bcenter,
                                        SCVector3 &point, SCVector3 &grad)
{
  SCVector3 r_g = point - centers[grad_center];
  SCVector3 r_b = point - centers[bcenter];
  SCVector3 r_gb = centers[grad_center] - centers[bcenter];
  double mag_r_g = r_g.norm();
  double mag_r_b = r_b.norm();
  double oorgb = oorab[grad_center][bcenter];
  double mu = (mag_r_g-mag_r_b)*oorgb;
  double a_gb = a_mat[grad_center][bcenter];
  double coef = 1.0-2.0*a_gb*mu;
  double r_g_coef;

  if (mag_r_g < 10.0 * DBL_EPSILON) r_g_coef = 0.0;
  else r_g_coef = -coef*oorgb/mag_r_g;
  int ixyz;
  for (ixyz=0; ixyz<3; ixyz++) grad[ixyz] = r_g_coef * r_g[ixyz];
  double r_gb_coef = coef*(mag_r_b - mag_r_g)*oorgb*oorgb*oorgb;
  for (ixyz=0; ixyz<3; ixyz++) grad[ixyz] += r_gb_coef * r_gb[ixyz];
}

// compute t(nu_ij)
double
BeckeIntegrationWeight::compute_t(int icenter, int jcenter, SCVector3 &point)
{
  // Cf. Johnson et al., JCP v. 98, p. 5612 (1993) (Appendix B)
  // NB: t is zero if s is zero
  
  SCVector3 r_i = point - centers[icenter];
  SCVector3 r_j = point - centers[jcenter];
  SCVector3 r_ij = centers[icenter] - centers[jcenter];
  double t;
  double mag_r_j = r_j.norm();
  double mag_r_i = r_i.norm();
  double mu = (mag_r_i-mag_r_j)*oorab[icenter][jcenter];
  if (mu >= 1.0-100*DBL_EPSILON) {
      t = 0.0;
      return t;
    }

  double a_ij = a_mat[icenter][jcenter];
  double nu = mu + a_ij*(1.-mu*mu);
  double s;
  if (mu <= -1.0) s = 1.0;
  else s = calc_s(nu);
  if (fabs(s) < 10*DBL_EPSILON) {
      t = 0.0;
      return t;
    }
  double p1 = 1.5*nu - 0.5*nu*nu*nu;
  double p2 = 1.5*p1 - 0.5*p1*p1*p1;

  t = -(27.0/16.0) * (1 - p2*p2) * (1 - p1*p1) * (1 - nu*nu) / s;
  return t;
}

void
BeckeIntegrationWeight::compute_grad_p(int grad_center, int bcenter,
                                          int wcenter, SCVector3&point,
                                          double p, SCVector3&grad)
{
  // the gradient of p is computed using the formulae from
  // Johnson et al., JCP v. 98, p. 5612 (1993) (Appendix B)

  if (grad_center == bcenter) {
      grad = 0.0;
      for (int dcenter=0; dcenter<n_integration_centers; dcenter++) {
          if (dcenter == bcenter) continue;
          SCVector3 grad_nu;
          compute_grad_nu(grad_center, dcenter, point, grad_nu);
          double t = compute_t(grad_center,dcenter,point);
          for (int ixyz=0; ixyz<3; ixyz++) grad[ixyz] += t * grad_nu[ixyz];
        }
    }
  else {
      SCVector3 grad_nu;
      compute_grad_nu(grad_center, bcenter, point, grad_nu);
      double t = compute_t(bcenter,grad_center,point);
      for (int ixyz=0; ixyz<3; ixyz++) grad[ixyz] = -t * grad_nu[ixyz];
    }
  grad *= p;
}

double
BeckeIntegrationWeight::w(int acenter, SCVector3 &point,
                          double *w_gradient)
{
  int icenter;
  double p_sum=0.0, p_a=0.0;
    
  for (icenter=0; icenter<n_integration_centers; icenter++) {
      double p_tmp = compute_p(icenter, point);
      if (icenter==acenter) p_a=p_tmp;
      p_sum += p_tmp;
    }
  double w_a = p_a/p_sum;

  if (w_gradient) {
      // w_gradient is computed using the formulae from
      // Johnson et al., JCP v. 98, p. 5612 (1993) (Appendix B)
      int i,j;
      for (i=0; i<n_integration_centers*3; i++ ) w_gradient[i] = 0.0;
//      fd_w(acenter, point, w_gradient);  // imbn commented out for debug
//        ExEnv::outn() << point << " ";
//        for (int i=0; i<n_integration_centers*3; i++) {
//            ExEnv::outn() << scprintf(" %10.6f", w_gradient[i]);
//          }
//        ExEnv::outn() << endl;
//      return w_a;  // imbn commented out for debug
      for (int ccenter = 0; ccenter < n_integration_centers; ccenter++) {
          // NB: for ccenter==acenter, use translational invariance
          // to get the corresponding component of the gradient
          if (ccenter != acenter) {
              SCVector3 grad_c_w_a;
              SCVector3 grad_c_p_a;
              compute_grad_p(ccenter, acenter, acenter, point, p_a, grad_c_p_a);
              for (i=0; i<3; i++) grad_c_w_a[i] = grad_c_p_a[i]/p_sum;
              for (int bcenter=0; bcenter<n_integration_centers; bcenter++) {
                  SCVector3 grad_c_p_b;
                  double p_b = compute_p(bcenter,point);
                  compute_grad_p(ccenter, bcenter, acenter, point, p_b,
                                 grad_c_p_b);
                  for (i=0; i<3; i++) grad_c_w_a[i] -= w_a*grad_c_p_b[i]/p_sum;
                }
              for (i=0; i<3; i++) w_gradient[ccenter*3+i] = grad_c_w_a[i];
            }
        }
      // fill in w_gradient for ccenter==acenter
      for (j=0; j<3; j++) {
          for (i=0; i<n_integration_centers; i++) {
              if (i != acenter) {
                  w_gradient[acenter*3+j] -= w_gradient[i*3+j];
                }
            }
        }
    }

  return w_a;
}

///////////////////////////////////////////////////
// RadialIntegrator
static ClassDesc RadialIntegrator_cd(
  typeid(RadialIntegrator),"RadialIntegrator",1,"public SavableState",
  0, 0, 0);

RadialIntegrator::RadialIntegrator(StateIn& s):
  SavableState(s)
{
}

RadialIntegrator::RadialIntegrator()
{
}

RadialIntegrator::RadialIntegrator(const Ref<KeyVal>& keyval)
{
}

RadialIntegrator::~RadialIntegrator()
{
}

void
RadialIntegrator::save_data_state(StateOut& s)
{
}

///////////////////////////////////////
//  AngularIntegrator

static ClassDesc AngularIntegrator_cd(
  typeid(AngularIntegrator),"AngularIntegrator",1,"public SavableState",
  0, 0, 0);

AngularIntegrator::AngularIntegrator(StateIn& s):
  SavableState(s)
{
}

AngularIntegrator::AngularIntegrator()
{
}

AngularIntegrator::AngularIntegrator(const Ref<KeyVal>& keyval)
{
}

AngularIntegrator::~AngularIntegrator()
{
}

void
AngularIntegrator::save_data_state(StateOut& s)
{
}

///////////////////////////////////////
//  EulerMaclaurinRadialIntegrator

static ClassDesc EulerMaclaurinRadialIntegrator_cd(
  typeid(EulerMaclaurinRadialIntegrator),"EulerMaclaurinRadialIntegrator",1,"public RadialIntegrator",
  0, create<EulerMaclaurinRadialIntegrator>, create<EulerMaclaurinRadialIntegrator>);

EulerMaclaurinRadialIntegrator::EulerMaclaurinRadialIntegrator(StateIn& s):
  SavableState(s),
  RadialIntegrator(s)
{
  s.get(nr_);
}

EulerMaclaurinRadialIntegrator::EulerMaclaurinRadialIntegrator()
{
  nr_ = 75;
}

EulerMaclaurinRadialIntegrator::EulerMaclaurinRadialIntegrator(int nr_points)
{
  nr_ = nr_points;
}

EulerMaclaurinRadialIntegrator::EulerMaclaurinRadialIntegrator(const Ref<KeyVal>& keyval):
  RadialIntegrator(keyval)
{
  nr_ = keyval->intvalue("nr", KeyValValueint(75));
}

EulerMaclaurinRadialIntegrator::~EulerMaclaurinRadialIntegrator()
{
}

void
EulerMaclaurinRadialIntegrator::save_data_state(StateOut& s)
{
  RadialIntegrator::save_data_state(s);
  s.put(nr_);
}

int
EulerMaclaurinRadialIntegrator::nr() const
{
  return nr_;
}

double
EulerMaclaurinRadialIntegrator::radial_value(int ir, int nr, double radii,
                                             double &multiplier)
{
  double q = (double)ir/(double)nr;
  double value = q/(1.-q);
  double r = radii*value*value;
  double dr_dq =  2.*radii*q*pow(1.-q,-3.);
  double dr_dqr2 = dr_dq*r*r;
  multiplier = dr_dqr2/nr;
  return r;
}

void
EulerMaclaurinRadialIntegrator::print(ostream &o) const
{
  o << indent
    << scprintf("%s: nr = %d", class_name(), nr()) << endl;
}

//////////////////////////////////////////////////////////////////////////
// LebedevLaikovIntegrator

static ClassDesc LebedevLaikovIntegrator_cd(
  typeid(LebedevLaikovIntegrator),"LebedevLaikovIntegrator",1,"public AngularIntegrator",
  0, create<LebedevLaikovIntegrator>, create<LebedevLaikovIntegrator>);

LebedevLaikovIntegrator::LebedevLaikovIntegrator(StateIn& s):
  SavableState(s),
  AngularIntegrator(s)
{
  s.get(npoint_);
  init(npoint_);
}

LebedevLaikovIntegrator::LebedevLaikovIntegrator()
{
  init(302);
}

LebedevLaikovIntegrator::LebedevLaikovIntegrator(int npoint)
{
  init(npoint);
}

LebedevLaikovIntegrator::LebedevLaikovIntegrator(const Ref<KeyVal>& keyval)
{
  KeyValValueint defnpoint(302);
  
  init(keyval->intvalue("n", defnpoint));
}

LebedevLaikovIntegrator::~LebedevLaikovIntegrator()
{
  delete [] x_;
  delete [] y_;
  delete [] z_;
  delete [] w_;
}

void
LebedevLaikovIntegrator::save_data_state(StateOut& s)
{
  AngularIntegrator::save_data_state(s);
  s.put(npoint_);
}

extern "C" {
    int Lebedev_Laikov_sphere (int N, double *X, double *Y, double *Z,
                               double *W);
    int Lebedev_Laikov_npoint (int lvalue);
}

int
LebedevLaikovIntegrator::nw(void) const
{
  return npoint_;
}

void
LebedevLaikovIntegrator::init(int n)
{
  // ExEnv::outn() << " LebedevLaikovIntegrator::init -> before x_, y_, z_, and w_ malloc's " << endl;
  // ExEnv::outn() << " n = " << n << endl;
  
  x_ = new double[n];
  y_ = new double[n];
  z_ = new double[n];
  w_ = new double[n];

  // ExEnv::outn() << " LebedevLaikovIntegrator::init -> nw_points = " << n << endl;
  
  npoint_ = Lebedev_Laikov_sphere(n, x_, y_, z_, w_);
  if (npoint_ != n) {
      ExEnv::outn() << class_name() << ": bad number of points given: " << n << endl;
      abort();
    }
}

int
LebedevLaikovIntegrator::num_angular_points(double r_value, int ir)
{
  if (ir == 0) return 1;
  return npoint_;
}

double
LebedevLaikovIntegrator
::angular_point_cartesian(int iangular, double r,
                          SCVector3 &integration_point) const
{
  integration_point.x() = r*x_[iangular];
  integration_point.y() = r*y_[iangular];
  integration_point.z() = r*z_[iangular];

  return 4.0*M_PI*w_[iangular];
}

void
LebedevLaikovIntegrator::print(ostream &o) const
{
  o << indent
    << scprintf("%s:  n = %d", class_name(), npoint_) << endl;
}

/////////////////////////////////
//  GaussLegendreAngularIntegrator

static ClassDesc GaussLegendreAngularIntegrator_cd(
  typeid(GaussLegendreAngularIntegrator),"GaussLegendreAngularIntegrator",1,"public AngularIntegrator",
  0, create<GaussLegendreAngularIntegrator>, create<GaussLegendreAngularIntegrator>);

GaussLegendreAngularIntegrator::GaussLegendreAngularIntegrator(StateIn& s):
  SavableState(s),
  AngularIntegrator(s)
{
  s.get(ntheta_);
  s.get(nphi_);
  s.get(Ktheta_);
  theta_quad_weights_ = new double[ntheta_];
  theta_quad_points_ = new double[ntheta_];
}

GaussLegendreAngularIntegrator::GaussLegendreAngularIntegrator()
{
  set_ntheta(16);
  set_nphi(32);
  set_Ktheta(5);
  int ntheta = get_ntheta();
  theta_quad_weights_ = new double [ntheta];
  theta_quad_points_ = new double [ntheta];
}

GaussLegendreAngularIntegrator::GaussLegendreAngularIntegrator(const Ref<KeyVal>& keyval)
{
  set_ntheta( keyval->intvalue("ntheta") );
  if (keyval->error() != KeyVal::OK) set_ntheta(16);
  set_nphi( keyval->intvalue("nphi") );
  if (keyval->error() != KeyVal::OK) set_nphi(2*get_ntheta());
  set_Ktheta( keyval->intvalue("Ktheta") );
  if (keyval->error() != KeyVal::OK) set_Ktheta(5);

  int ntheta = get_ntheta();
  theta_quad_weights_ = new double [ntheta];
  theta_quad_points_ = new double [ntheta];
}

GaussLegendreAngularIntegrator::~GaussLegendreAngularIntegrator()
{
  delete [] theta_quad_points_;
  delete [] theta_quad_weights_;
}

void
GaussLegendreAngularIntegrator::save_data_state(StateOut& s)
{
  AngularIntegrator::save_data_state(s);
  s.put(ntheta_);
  s.put(nphi_);
  s.put(Ktheta_);
}

int
GaussLegendreAngularIntegrator::get_ntheta(void) const
{
  return ntheta_;
}

void
GaussLegendreAngularIntegrator::set_ntheta(int i)
{
  ntheta_ = i;
}

int
GaussLegendreAngularIntegrator::get_nphi(void) const
{
  return nphi_;
}

void
GaussLegendreAngularIntegrator::set_nphi(int i)
{
  nphi_ = i;
}

int
GaussLegendreAngularIntegrator::get_Ktheta(void) const
{
  return Ktheta_;
}

void
GaussLegendreAngularIntegrator::set_Ktheta(int i)
{
  Ktheta_ = i;
}

int
GaussLegendreAngularIntegrator::get_ntheta_r(void) const
{
  return ntheta_r_;
}

void
GaussLegendreAngularIntegrator::set_ntheta_r(int i)
{
  ntheta_r_ = i;
}

int
GaussLegendreAngularIntegrator::get_nphi_r(void) const
{
  return nphi_r_;
}

void
GaussLegendreAngularIntegrator::set_nphi_r(int i)
{
  nphi_r_ = i;
}

int
GaussLegendreAngularIntegrator::get_Ktheta_r(void) const
{
  return Ktheta_r_;
}

void
GaussLegendreAngularIntegrator::set_Ktheta_r(int i)
{
  Ktheta_r_ = i;
}

int
GaussLegendreAngularIntegrator::nw(void) const
{
  return nphi_*ntheta_;
}

double
GaussLegendreAngularIntegrator::sin_theta(SCVector3 &point) const
{
  return sin(point.theta());
}

int
GaussLegendreAngularIntegrator::num_angular_points(double r_value,
                                                   int ir)
{
  int Ktheta, ntheta, ntheta_r;
  
  if (ir == 0) {
      set_ntheta_r(1);
      set_nphi_r(1);
    }
  else {
      Ktheta = get_Ktheta();
      ntheta = get_ntheta();
      ntheta_r= (int) (r_value*Ktheta*ntheta);
      set_ntheta_r(ntheta_r);
      if (ntheta_r > ntheta) set_ntheta_r(ntheta);
      if (ntheta_r < 6) set_ntheta_r(6);
      set_nphi_r(2*get_ntheta_r());
    }

  gauleg(0.0, M_PI, get_ntheta_r());

  return get_ntheta_r()*get_nphi_r();
}

void
GaussLegendreAngularIntegrator::gauleg(double x1, double x2, int n)
{
  int m,j,i;
  double z1,z,xm,xl,pp,p3,p2,p1;
  const double EPS = 10.0 * DBL_EPSILON;

  m=(n+1)/2;
  xm=0.5*(x2+x1);
  xl=0.5*(x2-x1);
  for (i=1;i<=m;i++)  {
      z=cos(M_PI*(i-0.25)/(n+0.5));
      do {
          p1=1.0;
          p2=0.0;
          for (j=1;j<=n;j++) {
              p3=p2;
              p2=p1;
              p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
	    }
          pp=n*(z*p1-p2)/(z*z-1.0);
          z1=z;
          z=z1-p1/pp;
	} while (fabs(z-z1) > EPS);
      theta_quad_points_[i-1]=xm-xl*z;
      theta_quad_points_[n-i]=xm+xl*z;
      theta_quad_weights_[i-1]=2.0*xl/((1.0-z*z)*pp*pp);
      theta_quad_weights_[n-i]=theta_quad_weights_[i-1];
    }
}

double
GaussLegendreAngularIntegrator
::angular_point_cartesian(int iangular, double r,
                          SCVector3 &integration_point) const
{
  int itheta, iphi, nphi_r;

  nphi_r = get_nphi_r();
  itheta = iangular/nphi_r;
  iphi = iangular - itheta*nphi_r;
  SCVector3 point;
  point.theta() = theta_quad_points_[itheta];
  point.phi() = (double) iphi/ (double) nphi_r * 2.0 * M_PI;
  point.r() = r;
  point.spherical_to_cartesian(integration_point);
  return ( sin_theta(point)*theta_quad_weights_[itheta]*2.0*M_PI/(double)nphi_r );
}

void
GaussLegendreAngularIntegrator::print(ostream &o) const
{
  o << indent << class_name() << ":" << endl;
  o << incindent;
  o << indent << scprintf("ntheta   = %5d", get_ntheta()) << endl;
  o << indent << scprintf("nphi     = %5d", get_nphi()) << endl;
  o << indent << scprintf("Ktheta   = %5d", get_Ktheta()) << endl;
  o << decindent;
}

//////////////////////////////////////////////
//  RadialAngularIntegratorThread

class RadialAngularIntegratorThread: public DenIntegratorThread {
  protected:
    SCVector3 *centers_;
    int *nr_;
    double *atomic_radius_;
    Molecule *mol_;
    RadialAngularIntegrator *ra_integrator_;
    IntegrationWeight *weight_;
    int point_count_total_;
    double total_density_;
  public:
    RadialAngularIntegratorThread(int ithread, int nthread,
                                  RadialAngularIntegrator *integrator,
                                  DenFunctional *func,
                                  const Ref<BatchElectronDensity> &den,
                                  int linear_scaling, int use_dmat_bound,
                                  double accuracy,
                                  int compute_potential_integrals,
                                  int need_nuclear_gradient);
    ~RadialAngularIntegratorThread();
    void run();
    double total_density() { return total_density_; }
    int point_count() { return point_count_total_; }
};

RadialAngularIntegratorThread
::RadialAngularIntegratorThread(int ithread, int nthread,
                                RadialAngularIntegrator *integrator,
                                DenFunctional *func,
                                const Ref<BatchElectronDensity> &den,
                                int linear_scaling, int use_dmat_bound,
                                double accuracy,
                                int compute_potential_integrals,
                                int need_nuclear_gradient):
  DenIntegratorThread(ithread,nthread,
                      integrator, func,
                      den,
                      linear_scaling, use_dmat_bound,
                      accuracy,
                      compute_potential_integrals,
                      need_nuclear_gradient)
{
  int icenter;
  ra_integrator_ = integrator;
  int deriv_order = (need_nuclear_gradient==0?0:1);

  mol_ = integrator_->basis()->molecule().pointer();

  weight_ = ra_integrator_->weight().pointer();

  nr_ = new int[n_integration_center_];
  
  for (icenter=0; icenter<n_integration_center_; icenter++) {
      int iatom = mol_->non_q_atom(icenter);
      nr_[icenter]
	= ra_integrator_->get_radial_grid(mol_->Z(iatom),deriv_order)->nr();
    }

  centers_ = new SCVector3[n_integration_center_];
  for (icenter=0; icenter<n_integration_center_; icenter++) {
      int iatom = mol_->non_q_atom(icenter);
      centers_[icenter].x() = mol_->r(iatom,0);
      centers_[icenter].y() = mol_->r(iatom,1);
      centers_[icenter].z() = mol_->r(iatom,2);
    }

  atomic_radius_ = new double[n_integration_center_];
  for (icenter=0; icenter<n_integration_center_; icenter++) {
      int iatom = mol_->non_q_atom(icenter);
      atomic_radius_[icenter] = get_radius(mol_, iatom);
    }

  point_count_total_ = 0;
  total_density_ = 0.0;
}

RadialAngularIntegratorThread::~RadialAngularIntegratorThread()
{
  delete[] centers_;
  delete[] atomic_radius_;
  delete[] nr_;
}

void
RadialAngularIntegratorThread::run()
{
  int icenter;
  int nangular;
  int ir, iangular;           // Loop indices for diff. integration dim
  int point_count;            // Counter for # integration points per center
  int nr;

  SCVector3 center;           // Cartesian position of center
  SCVector3 integration_point;

  double w,radial_multiplier,angular_multiplier;
  int deriv_order = (nuclear_gradient_==0?0:1);
        
  int parallel_counter = 0;

  for (icenter=0; icenter < n_integration_center_; icenter++) {
      int iatom = mol_->non_q_atom(icenter);
      point_count=0;
      center = centers_[icenter];
      // get current radial grid: depends on convergence threshold
      RadialIntegrator *radial
          = ra_integrator_->get_radial_grid(mol_->Z(iatom), deriv_order);
      nr = radial->nr();
      for (ir=0; ir < nr; ir++) {
          if (! (parallel_counter++%nthread_ == ithread_)) continue;
          double r = radial->radial_value(ir, nr, atomic_radius_[icenter],
                                          radial_multiplier);
          // get current angular grid: depends on radial point and threshold
          AngularIntegrator *angular
              = ra_integrator_->get_angular_grid(r, atomic_radius_[icenter],
                                                 mol_->Z(iatom),
						 deriv_order);
          nangular = angular->num_angular_points(r/atomic_radius_[icenter],ir);
          for (iangular=0; iangular<nangular; iangular++) {
              angular_multiplier
                  = angular->angular_point_cartesian(iangular,r,
                                                     integration_point);
              integration_point += center;
              w=weight_->w(icenter, integration_point, w_gradient_);
              point_count++;
              double multiplier = angular_multiplier * radial_multiplier;
              total_density_
                  += w * multiplier
                  * do_point(iatom, integration_point,
                             w, multiplier,
                             nuclear_gradient_, f_gradient_, w_gradient_);
            }
        }
      point_count_total_ += point_count;
    }
}

//////////////////////////////////////////////
//  RadialAngularIntegrator

static ClassDesc RadialAngularIntegrator_cd(
  typeid(RadialAngularIntegrator),"RadialAngularIntegrator",1,"public DenIntegrator",
  0, create<RadialAngularIntegrator>, create<RadialAngularIntegrator>);

RadialAngularIntegrator::RadialAngularIntegrator(StateIn& s):
  SavableState(s),
  DenIntegrator(s)
{
  s.get(natomic_rows_);
  s.get(max_gridtype_);
  s.get(prune_grid_);
  s.get(gridtype_);
  s.get(npruned_partitions_);
  s.get(dynamic_grids_);

//  ExEnv::outn() << "natomic_rows_ = " << natomic_rows_ << endl;
//  ExEnv::outn() << "max_gridtype_ = " << max_gridtype_ << endl;
//  ExEnv::outn() << "prune_grid_ = " << prune_grid_ << endl;
//  ExEnv::outn() << "gridtype_ = " << gridtype_ << endl;
//  ExEnv::outn() << "npruned_partitions_ = " << npruned_partitions_ << endl;
//  ExEnv::outn() << "dynamic_grids_ = " << dynamic_grids_ << endl;
  

//  ExEnv::outn() << "In StateIn Constructor!" << endl;
  weight_ = new BeckeIntegrationWeight;

  int i;
  grid_accuracy_ = new double[max_gridtype_];
  grid_accuracy_[0] = 1e-4;
  for (i=1; i<max_gridtype_; i++) grid_accuracy_[i] = grid_accuracy_[i-1]*1e-1;

  Alpha_coeffs_ = new_c_array2(natomic_rows_,npruned_partitions_-1,
                                       double(0));
  s.get_array_double(Alpha_coeffs_[0], natomic_rows_*(npruned_partitions_-1));

  radial_user_ << SavableState::restore_state(s);
  angular_user_ << SavableState::restore_state(s);

  init_default_grids();    
  set_grids();

}

RadialAngularIntegrator::RadialAngularIntegrator()
{
  weight_  = new BeckeIntegrationWeight;
  
  init_parameters();
  init_default_grids();
  set_grids();

}

RadialAngularIntegrator::RadialAngularIntegrator(const Ref<KeyVal>& keyval):
  DenIntegrator(keyval)
{
  
  radial_user_ << keyval->describedclassvalue("radial");
  angular_user_ << keyval->describedclassvalue("angular");

  weight_ << keyval->describedclassvalue("weight");
  if (weight_ == 0) weight_ = new BeckeIntegrationWeight;
//  ExEnv::outn() << "In Ref<KeyVal> Constructor" << endl;
  
  init_parameters(keyval);
  init_default_grids();
  set_grids();

}

RadialAngularIntegrator::~RadialAngularIntegrator()
{
  delete_c_array2(Alpha_coeffs_);
  delete_c_array2(radial_grid_);
  delete_c_array3(angular_grid_);
  delete_c_array2(nr_points_);
  delete[] xcoarse_l_;
  delete[] grid_accuracy_;
}

void
RadialAngularIntegrator::save_data_state(StateOut& s)
{
  DenIntegrator::save_data_state(s);
  s.put(natomic_rows_);
  s.put(max_gridtype_);
  s.put(prune_grid_);
  s.put(gridtype_);
//  s.put(nr_points_[0], natomic_rows_*max_gridtype_);
//  s.put(xcoarse_l_, natomic_rows_);
  s.put(npruned_partitions_);
  s.put(dynamic_grids_);
  s.put_array_double(Alpha_coeffs_[0],natomic_rows_*(npruned_partitions_-1));
  
//  ExEnv::outn() << "natomic_rows_ = " << natomic_rows_ << endl;
//  ExEnv::outn() << "max_gridtype_ = " << max_gridtype_ << endl;
//  ExEnv::outn() << "prune_grid_ = " << prune_grid_ << endl;
//  ExEnv::outn() << "gridtype_ = " << gridtype_ << endl;
//  ExEnv::outn() << "npruned_partitions_ = " << npruned_partitions_ << endl;
//  ExEnv::outn() << "dynamic_grids_ = " << dynamic_grids_ << endl;
  
  SavableState::save_state(radial_user_.pointer(),s);
  SavableState::save_state(angular_user_.pointer(),s);
}

void
RadialAngularIntegrator::init_parameters(void)
{

  prune_grid_ = 1;
  gridtype_ = 3;
  npruned_partitions_ = 5;
  dynamic_grids_ = 1;
  max_gridtype_ = 6;
  natomic_rows_ = 5;
  grid_accuracy_ = new double[max_gridtype_];
  
  int i;
  grid_accuracy_[0] = 1e-4;
  for (i=1; i<max_gridtype_; i++) grid_accuracy_[i] = grid_accuracy_[i-1]*1e-1;
    
  init_pruning_coefficients();

  // ExEnv::outn() << "gridtype_ = " << gridtype_ << endl;
  
}

void
RadialAngularIntegrator::init_parameters(const Ref<KeyVal>& keyval)
{
  max_gridtype_ = 6;
  natomic_rows_ = 5;

  if (keyval->exists("grid")) {
      std::string grid = keyval->stringvalue("grid");
      if (grid == "xcoarse")        gridtype_ = 0;
      else if (grid == "coarse")    gridtype_ = 1;
      else if (grid == "medium")    gridtype_ = 2;
      else if (grid == "fine")      gridtype_ = 3;
      else if (grid == "xfine")     gridtype_ = 4;
      else if (grid == "ultrafine") gridtype_ = 5;
      else {
          ExEnv::out0()
                       << indent
                       << "ERROR: grid = \"" << grid << "\" not recognized."
                       << endl
                       << indent
                       << "Choices are: xcoarse, coarse, medium, fine, xfine."
                       << endl
                       << indent
                       << "The default is grid = fine."
                       << endl;
          abort();
        }

    }
  else {
      gridtype_ = 3;
    }


  
  //ExEnv::outn() << " gridtype = " << gridtype_ << endl;
  //ExEnv::outn() << " max_gridtype = " << max_gridtype_ << endl;
  dynamic_grids_ = keyval->intvalue("dynamic");
  if (keyval->error() != KeyVal::OK) dynamic_grids_ = 1;
  grid_accuracy_ = new double[max_gridtype_];
  //ExEnv::outn() << "init_parameters:: max_gridtype_ = " << max_gridtype_;
  
  int i;
  grid_accuracy_[0] = 1e-4;
  for (i=1; i<max_gridtype_; i++) grid_accuracy_[i] = grid_accuracy_[i-1]*1e-1;
  
  init_pruning_coefficients(keyval);

}


void
RadialAngularIntegrator::set_grids(void)
{
  int i, j, k;

  radial_grid_ = new_c_array2(natomic_rows_,gridtype_+1,
                              Ref<RadialIntegrator>());
  angular_grid_ = new_c_array3(natomic_rows_, npruned_partitions_,
                               gridtype_+1, Ref<AngularIntegrator>());
  
  int prune_formula_1[5] = {26, 18, 12, 0, 12};  // for H to Ne
  int prune_formula_2[5] = {36, 24, 12, 0, 12};  // for Na and up

  double prune_factor[5] = {0.1, 0.4, 0.6, 1.0, 0.6};

  if (npruned_partitions_ == 1) {
      prune_factor[0] = 1.0;
      prune_formula_1[0] = 0;
      prune_formula_2[0] = 0;
    }
  else if (npruned_partitions_ != 5) {
      ExEnv::errn() << "RadialAngularIntegrator::set_grids: "
                   << "npruned_partitations must be 1 or 5" << endl;
      abort();
    }

  for (i=0; i<natomic_rows_; i++) {
      for (j=0; j<=gridtype_; j++)
          radial_grid_[i][j]
              = new EulerMaclaurinRadialIntegrator(nr_points_[i][j]);

      for (j=0; j<npruned_partitions_; j++) {
          for (k=0; k<=gridtype_; k++) {
              int grid_l = xcoarse_l_[i] + angular_grid_offset(k);
              int use_l;

              // the constant shift depends on the row and the partition
              int prune_formula;
              if (i<=1) prune_formula = prune_formula_1[j]; // H to Ne
              else prune_formula = prune_formula_2[j];      // Na and up

              // Compute the l to be used from both the constant shift and
              // the multiplicative factor method.  Use the method giving
              // the highest l.
              int use_l_formula = grid_l - prune_formula;
              int use_l_factor = int(grid_l*prune_factor[j] + 0.5);
              if (use_l_formula > use_l_factor) use_l = use_l_formula;
              else use_l = use_l_factor;

              angular_grid_[i][j][k]
                  = new LebedevLaikovIntegrator(Lebedev_Laikov_npoint(use_l));
//                ExEnv::outn() << " angular_grid_["
//                             << i << "]["
//                             << j << "]["
//                             << k
//                             << "]->nw = " << angular_grid_[i][j][k]->nw()
//                             << " xc_l = " << xcoarse_l_[i]
//                             << " off = " << angular_grid_offset(k)
//                             << " grid_l = " << grid_l
//                             << " use_l = " << use_l
//                             << endl;
            }
        }
    }
}

void
RadialAngularIntegrator::init_pruning_coefficients(void)
{
  // Set up Alpha arrays for pruning
  //ExEnv::outn() << "npruned_partitions = " << npruned_partitions_ << endl;
  //ExEnv::outn() << "natomic_rows = " << natomic_rows_ << endl;
  int num_boundaries = npruned_partitions_-1;
  Alpha_coeffs_ = new_zero_c_array2(natomic_rows_, num_boundaries,
                                            double(0));

  // set pruning cutoff variables - Alpha -> radial shell cutoffs
  init_alpha_coefficients();
   
}

void
RadialAngularIntegrator::init_pruning_coefficients(const Ref<KeyVal>& keyval)
{
  int i, j;

  prune_grid_ = keyval->booleanvalue("prune_grid");
  if (keyval->error() != KeyVal::OK) prune_grid_ = 1;

  // ExEnv::outn() << "prune_grid = " << prune_grid_ << endl;
  
  // Need to generalize this to parse input for any number of grids
  if (prune_grid_) {
      npruned_partitions_ = keyval->count("angular_points");
      if (keyval->error() != KeyVal::OK) npruned_partitions_ = 5;

      // set pruning cutoff variables - Alpha -> radial shell cutoffs
      int num_boundaries = npruned_partitions_-1;
      Alpha_coeffs_ = new_zero_c_array2(natomic_rows_, num_boundaries,
                                                double(0));
      int alpha_rows = keyval->count("alpha_coeffs");
      if (keyval->error() != KeyVal::OK) {
          if (npruned_partitions_ != 5) {
              ExEnv::outn() << " RadialAngularIntegrator:: Need to supply alpha coefficients "
                   << "for the " << num_boundaries << " partition boundaries " << endl;
              abort();
            }
          init_alpha_coefficients();
        }
      else { // alpha coefficients defined in input
          int check;
          for (i=0; i<alpha_rows; i++) {
              check = keyval->count("alpha_coeffs", i);
              if (check != num_boundaries) {
                  ExEnv::outn() << "RadialAngularIntegrator:: Number of alpha coefficients does "
                       << "not match the number of boundaries (" << check << " != "
                       << num_boundaries << ")" << endl;
                  abort();
                }
              for (j=0; j<num_boundaries; j++) 
                  Alpha_coeffs_[i][j] = keyval->doublevalue("alpha_coeffs", i, j);
            }
        }
    }
  else {
      npruned_partitions_ = 1;
      Alpha_coeffs_ = new_zero_c_array2(natomic_rows_,0,
                                                double(0));
    }
}

void
RadialAngularIntegrator::init_alpha_coefficients(void)
{
  // assumes Alpha_coeffs_ is allocated and zeroed.

  Alpha_coeffs_[0][0] = 0.25;   Alpha_coeffs_[0][1] = 0.5;
  Alpha_coeffs_[0][2] = 0.9;    Alpha_coeffs_[0][3] = 4.5;
  Alpha_coeffs_[1][0] = 0.1667; Alpha_coeffs_[1][1] = 0.5;
  Alpha_coeffs_[1][2] = 0.8;    Alpha_coeffs_[1][3] = 4.5;
  Alpha_coeffs_[2][0] = 0.1;    Alpha_coeffs_[2][1] = 0.4;
  Alpha_coeffs_[2][2] = 0.7;    Alpha_coeffs_[2][3] = 2.5;

  //  No pruning for atoms past second row
  int i;
  for (i=3; i<natomic_rows_; i++) Alpha_coeffs_[i][2] = DBL_MAX;

}

void
RadialAngularIntegrator::init_default_grids(void)
{
  xcoarse_l_ = new int[natomic_rows_];

  nr_points_ = new_c_array2(natomic_rows_,max_gridtype_,int(0));

  // Set angular momentum level of reference xcoarse grids for each atomic row
  xcoarse_l_[0] = 17; xcoarse_l_[1] = 17; xcoarse_l_[2] = 21;
  xcoarse_l_[3] = 25; xcoarse_l_[4] = 31;

  // Set number of radial points for each atomic row and grid type
  nr_points_[0][0] = 30; nr_points_[0][1] = 50; nr_points_[0][2] = 75;
  nr_points_[0][3] = 85; nr_points_[0][4] = 100; nr_points_[0][5] = 140;

  nr_points_[1][0] = 30; nr_points_[1][1] = 50; nr_points_[1][2] = 75;
  nr_points_[1][3] = 85; nr_points_[1][4] = 100; nr_points_[1][5] = 140;

  nr_points_[2][0] = 45; nr_points_[2][1] = 75; nr_points_[2][2] = 95;
  nr_points_[2][3] = 110; nr_points_[2][4] = 125; nr_points_[2][5] = 175;

  nr_points_[3][0] = 75; nr_points_[3][1] = 95; nr_points_[3][2] = 110;
  nr_points_[3][3] = 130; nr_points_[3][4] = 160; nr_points_[3][5] = 210;

  nr_points_[4][0] = 105; nr_points_[4][1] = 130; nr_points_[4][2] = 155;
  nr_points_[4][3] = 180; nr_points_[4][4] = 205; nr_points_[4][5] = 235;

  // prune_grid_ = 1; npruned_partitions_ = 5; gridtype_ = 2;
}

int
RadialAngularIntegrator::angular_grid_offset(int gridtype)
{
  switch (gridtype) {
  case 0: return 0;
  case 1: return 6;
  case 2: return 12;
  case 3: return 18;
  case 4: return 30;
  case 5: return 42;
  default: abort();
    }
  return 0;
}

RadialIntegrator *
RadialAngularIntegrator::get_radial_grid(int charge, int deriv_order)
{
  if (radial_user_ == 0) {

      int select_grid;
      
      if (dynamic_grids_ && deriv_order == 0)
	  select_grid = select_dynamic_grid();
      else select_grid = gridtype_;
      //ExEnv::out << "RAI::get_radial_grid -> select_grid = " << select_grid;
      
      if (charge<3) return radial_grid_[0][select_grid].pointer();
      else if (charge<11) return radial_grid_[1][select_grid].pointer();
      else if (charge<19) return radial_grid_[2][select_grid].pointer();
      else if (charge<37) return radial_grid_[3][select_grid].pointer();
      else if (charge<55) return radial_grid_[4][select_grid].pointer();
      else {
          ExEnv::outn() << " No default radial grids for atomic charge " << charge << endl;
          abort();
        }
    }

  return radial_user_.pointer();

}

int
RadialAngularIntegrator::select_dynamic_grid(void)
{
  double accuracy = get_accuracy();
  // accurate_grid = gridtype_ to get original non-dynamic version
  // ExEnv::out << " accuracy = " << accuracy << endl;
  int select_grid;
  int i;

  if (accuracy >= grid_accuracy_[0]) select_grid=0;
  else if (accuracy <= grid_accuracy_[gridtype_]) select_grid=gridtype_;
  else {
      for (i=gridtype_; i>=0; i--)
          if (accuracy >= grid_accuracy_[i]) select_grid=i;
    }
      
  // ExEnv::out << " select_grid = " << select_grid << endl;
  return select_grid;
}  

int
RadialAngularIntegrator::get_atomic_row(int i)
{
  if (i<3) return 0;
  else if (i<11) return 1;
  else if (i<19) return 2;
  else if (i<37) return 3;
  else if (i<55) return 4;
  else if (i<87) return 5;

  ExEnv::outn() << " RadialAngularIntegrator::get_atomic_row: Z too large: "
               << i << endl;
  abort();
  return 0;
}

AngularIntegrator *
RadialAngularIntegrator::get_angular_grid(double radius, double atomic_radius,
                                          int Z, int deriv_order)
{
  int atomic_row, i;

  int select_grid;
  if (dynamic_grids_ && deriv_order == 0) select_grid = select_dynamic_grid();
  else select_grid = gridtype_;

  //ExEnv::out << "RAI::get_angular_grid -> select_grid = " << select_grid;
  //ExEnv::out << " prune_grid_ = " << prune_grid_ << endl;
  atomic_row = get_atomic_row(Z);
  if (angular_user_ == 0) {
      // Seleted Alpha's
      double *Alpha = Alpha_coeffs_[atomic_row];
      // gridtype_ will need to be adjusted for dynamic grids  
      for (i=0; i<npruned_partitions_-1; i++) {
          if (radius/atomic_radius < Alpha[i]) {
              return angular_grid_[atomic_row][i][select_grid].pointer();
            }
        }
      return angular_grid_[atomic_row][npruned_partitions_-1][select_grid]
          .pointer();
    }
  else {
      return angular_user_.pointer();
    }
}

void
RadialAngularIntegrator::integrate(const Ref<DenFunctional> &denfunc,
                              const RefSymmSCMatrix& densa,
                              const RefSymmSCMatrix& densb,
                              double *nuclear_gradient)
{
  int i;

  Timer tim("integrate");

  init_integration(denfunc, densa, densb, nuclear_gradient);

  weight_->init(basis()->molecule(), DBL_EPSILON);

  int me = messagegrp_->me();
  int nthread = threadgrp_->nthread();
  int nthread_overall = nthread;
  messagegrp_->sum(nthread_overall);
  int ithread_overall = 0;
  if (me > 0) {
      messagegrp_->recv(me - 1,ithread_overall);
    }
  if (me < messagegrp_->n() - 1) {
      int ithread_overall_next = ithread_overall + nthread;
      messagegrp_->send(me + 1, ithread_overall_next);
    }

  // create threads
  //cout << "creating test lock" << endl;
  //Ref<ThreadLock> reflock = threadgrp_->new_lock();
  //tlock = reflock.pointer();
  RadialAngularIntegratorThread **threads =
      new RadialAngularIntegratorThread*[nthread];
  for (i=0; i<nthread; i++) {
      Ref<BatchElectronDensity> bed = new BatchElectronDensity(den_, true);
      if (nuclear_gradient != 0) {
          bed->set_need_basis_gradient(true);
          if (denfunc->need_density_gradient()) bed->set_need_basis_hessian(true);
        }
      threads[i] = new RadialAngularIntegratorThread(
          i + ithread_overall, nthread_overall,
          this, denfunc.pointer(),
          bed,
          linear_scaling_, use_dmat_bound_,
          accuracy_, compute_potential_integrals_,
          nuclear_gradient != 0);
      threadgrp_->add_thread(i, threads[i]);
    }

  //for (i=0; i<nthread; i++) {
  //    ExEnv::outn() << "Intial Thread " << i
  //                 << ": point count = " << threads[i]->point_count()
  //                 << " total density = " << threads[i]->total_density()
  //                 << " value = " << threads[i]->value()
  //                 << endl;
  //}

  // run threads
  threadgrp_->start_threads();
  threadgrp_->wait_threads();
  //for (i=0; i<nthread; i++) {
  //    cout << "Running thread " << i << endl;
  //    threads[i]->run();
  //  }

  // sum results
  int point_count_total = 0;
  double total_density = 0.0;
  value_ = 0.0;
  for (i=0; i<nthread; i++) {
      //ExEnv::outn() << "Thread " << i
      //             << ": point count = " << threads[i]->point_count()
      //             << " total density = " << threads[i]->total_density()
      //             << " value = " << threads[i]->value()
      //             << endl;
      point_count_total += threads[i]->point_count();
      total_density += threads[i]->total_density();
      value_ += threads[i]->value();
      if (compute_potential_integrals_) {
          int ntri = (nbasis_*(nbasis_+1))/2;
          double *alpha_vmat_i = threads[i]->alpha_vmat();
          for (int j=0; j<ntri; j++) alpha_vmat_[j] += alpha_vmat_i[j];
          if (spin_polarized_) {
              double *beta_vmat_i = threads[i]->beta_vmat();
              for (int j=0; j<ntri; j++) beta_vmat_[j] += beta_vmat_i[j];
            }
        }
      if (nuclear_gradient != 0) {
          int natom3 = 3 * basis()->molecule()->natom();
          double *th_nuclear_gradient = threads[i]->nuclear_gradient();
          for (int j=0; j<natom3; j++) {
              nuclear_gradient[j] += th_nuclear_gradient[j];
            }
        }
    }

  threadgrp_->delete_threads();
  delete[] threads;

  messagegrp_->sum(point_count_total);
  messagegrp_->sum(total_density);
  done_integration();
  weight_->done();

  ExEnv::out0() << indent
       << "Total integration points = " << point_count_total << endl;
  ExEnv::out0() << indent
               << "Integrated electron density = "
               << scprintf("%16.12f", total_density)
               << endl;

  tim.exit("integrate");
}

void
RadialAngularIntegrator::print(ostream &o) const
{
  o << indent << class_name() << ":" << endl;
  o << incindent;
  if (radial_user_) {
      o << indent << "User defined radial grid:" << endl;
      o << incindent;
      radial_user_->print(o);
      o << decindent;
    }
  if (angular_user_) {
      o << indent << "User defined angular grid:" << endl;
      o << incindent;
      angular_user_->print(o);
      o << decindent;
    }
  if (angular_user_ == 0 || radial_user_ == 0) {
      if (prune_grid_) o << indent << "Pruned ";
      switch (gridtype_) {
      case 0: o << "xcoarse"; break;
      case 1: o << "coarse"; break;
      case 2: o << "medium"; break;
      case 3: o << "fine"; break;
      case 4: o << "xfine"; break;
      case 5: o << "ultrafine"; break;
      default: o << "unknown"; break;
        }
      o << " grid employed" << endl;
    }
  
  o << decindent;
}

/////////////////////////////////////////////////////////////////////////////

}

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

