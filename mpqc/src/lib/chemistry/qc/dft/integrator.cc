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

#ifdef __GNUC__
#pragma implementation
#endif

#include <math.h>

#include <util/misc/timer.h>
#include <util/misc/formio.h>
#include <util/state/stateio.h>
#include <chemistry/qc/dft/integrator.h>

#define COUNT_CONTRIBUTIONS 0 // not mt-safe if 1
//#define CHECK_ALIGN(v) if(int(&v)&7)cout<<"Bad Alignment: "<< ## v <<endl;
#define CHECK_ALIGN(v)

///////////////////////////////////////////////////////////////////////////
// DenIntegrator

#define CLASSNAME DenIntegrator
#define PARENTS public SavableState
#include <util/state/statei.h>
#include <util/class/classia.h>
void *
DenIntegrator::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

DenIntegrator::DenIntegrator(StateIn& s):
  SavableState(s)
{
}

DenIntegrator::DenIntegrator()
{
  contrib_ = 0;
  contrib_bf_ = 0;
  bs_values_ = 0;
  bsg_values_ = 0;
  bsh_values_ = 0;
  alpha_dmat_ = 0;
  beta_dmat_ = 0;
  dmat_bound_ = 0;
  alpha_vmat_ = 0;
  beta_vmat_ = 0;
  compute_potential_integrals_ = 0;
  accuracy_ = DBL_EPSILON;
  linear_scaling_ = 1;
  use_dmat_bound_ = 1;
}

DenIntegrator::DenIntegrator(const RefKeyVal& keyval)
{
  contrib_ = 0;
  contrib_bf_ = 0;
  bs_values_ = 0;
  bsg_values_ = 0;
  bsh_values_ = 0;
  alpha_dmat_ = 0;
  beta_dmat_ = 0;
  dmat_bound_ = 0;
  alpha_vmat_ = 0;
  beta_vmat_ = 0;
  compute_potential_integrals_ = 0;
  accuracy_ = DBL_EPSILON;

  linear_scaling_ = keyval->booleanvalue("linear_scaling",
                                         KeyValValueboolean(1));
  use_dmat_bound_ = keyval->booleanvalue("use_dmat_bound",
                                         KeyValValueboolean(1));
}

DenIntegrator::~DenIntegrator()
{
  delete[] contrib_;
  delete[] contrib_bf_;
  delete[] bs_values_;
  delete[] bsg_values_;
  delete[] bsh_values_;
  delete[] alpha_dmat_;
  delete[] beta_dmat_;
  delete[] dmat_bound_;
  delete[] alpha_vmat_;
  delete[] beta_vmat_;
}

void
DenIntegrator::save_data_state(StateOut& s)
{
  cout << class_name() << ": cannot save state" << endl;
  abort();
}

void
DenIntegrator::set_compute_potential_integrals(int i)
{
  compute_potential_integrals_=i;
}

void
DenIntegrator::init(const RefWavefunction &wfn)
{
  wfn_ = wfn;
  if (linear_scaling_) {
      cout << node0 << indent << "Initializing ShellExtent" << endl;
      extent_ = new ShellExtent;
      extent_->init(wfn_->basis());
      cout << node0 << indent
           << "  nshell = " << wfn_->basis()->nshell() << endl;
      int ncell = extent_->n(0)*extent_->n(1)*extent_->n(2);
      cout << node0 << indent << "  ncell = " << ncell << endl;
      int maxval = 0;
      double ave = 0;
      for (int i=0; i<extent_->n(0); i++) {
          for (int j=0; j<extent_->n(1); j++) {
              for (int k=0; k<extent_->n(2); k++) {
                  const ArrayExtentData &e
                      = extent_->contributing_shells(i,j,k);
                  int val = e.size();
                  if (val > maxval) maxval = val;
                  ave += val;
                }
            }
        }
      ave /= ncell;
      cout << node0 << indent << "  ave nsh/cell = " << ave << endl;
      cout << node0 << indent << "  max nsh/cell = " << maxval << endl;
    }
}

void
DenIntegrator::done()
{
  wfn_ = 0;
  extent_ = 0;
}

void
DenIntegrator::init_integration(const RefDenFunctional &func,
                                const RefSymmSCMatrix& densa,
                                const RefSymmSCMatrix& densb,
                                double *nuclear_gradient)
{
  int i;
  value_ = 0.0;

  wfn_->basis()->set_integral(wfn_->integral());

  func->set_compute_potential(
      compute_potential_integrals_ || nuclear_gradient != 0);

  need_gradient_ = func->need_density_gradient();
  need_hessian_ = 0;
  if (need_hessian_) need_gradient_ = 1;

  spin_polarized_ = wfn_->spin_polarized();
  func->set_spin_polarized(spin_polarized_);

  natom_ = wfn_->molecule()->natom();

  nshell_ = wfn_->basis()->nshell();
  nbasis_ = wfn_->basis()->nbasis();

  delete[] contrib_;
  contrib_ = new int[nshell_];

  delete[] contrib_bf_;
  contrib_bf_ = new int[nbasis_];

  delete[] bs_values_;
  bs_values_ = new double[nbasis_];

  delete[] bsg_values_;
  delete[] bsh_values_;
  bsg_values_=0;
  bsh_values_=0;
  if (need_gradient_ || nuclear_gradient) bsg_values_ = new double[3*nbasis_];
  if (need_hessian_ || (need_gradient_ && nuclear_gradient))
      bsh_values_ = new double[6*nbasis_];
   
  delete[] alpha_dmat_;
  RefSymmSCMatrix adens = densa;
  if (adens.null())
      adens = wfn_->alpha_ao_density();
  alpha_dmat_ = new double[(nbasis_*(nbasis_+1))/2];
  adens->convert(alpha_dmat_);

  delete[] beta_dmat_;
  beta_dmat_ = 0;
  if (spin_polarized_) {
      RefSymmSCMatrix bdens = densb;
      if (bdens.null())
          bdens = wfn_->beta_ao_density();
      beta_dmat_ = new double[(nbasis_*(nbasis_+1))/2];
      bdens->convert(beta_dmat_);
    }

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

  delete[] dmat_bound_;
  dmat_bound_ = new double[(nshell_*(nshell_+1))/2];
  int ij = 0;
  for (i=0; i<nshell_; i++) {
      int ni = wfn_->basis()->shell(i).nfunction();
      for (int j=0; j<=i; j++,ij++) {
          int nj = wfn_->basis()->shell(j).nfunction();
          double bound = 0.0;
          int ibf = wfn_->basis()->shell_to_function(i);
          for (int k=0; k<ni; k++,ibf++) {
              int lmax = nj-1;
              if (i==j) lmax = k;
              int jbf = wfn_->basis()->shell_to_function(j);
              int ijbf = (ibf*(ibf+1))/2 + jbf;
              for (int l=0; l<=lmax; l++,ijbf++) {
                  double a = fabs(alpha_dmat_[ijbf]);
                  if (a > bound) bound = a;
                  if (beta_dmat_) {
                      double b = fabs(beta_dmat_[ijbf]);
                      if (b > bound) bound = b;
                    }
                }
            }
          dmat_bound_[ij] = bound;
        }
    }
}

void
DenIntegrator::done_integration()
{
  RefMessageGrp msg = MessageGrp::get_default_messagegrp();

  msg->sum(value_);
  if (compute_potential_integrals_) {
      int ntri = (nbasis_*(nbasis_+1))/2;
      msg->sum(alpha_vmat_,ntri);
      if (spin_polarized_) {
          msg->sum(beta_vmat_,ntri);
        }
    }
}

inline static double
norm(double v[3])
{
  double x,y,z;
  return sqrt((x=v[0])*x + (y=v[1])*y + (z=v[2])*z);
}

inline static double
dot(double v[3], double w[3])
{
  return v[0]*w[0] + v[1]*w[1] + v[2]*w[2];
}

void
DenIntegrator::get_density(double *dmat, PointInputData::SpinData &d)
{
  int i, j, ish, jsh;

  tim_enter("get_density");

  GaussianBasisSet *basis = wfn_->basis().pointer();

  double tmp = 0.0;
  double densij;
  double bvi, bvix, bviy, bviz;
  double bvixx, bviyx, bviyy, bvizx, bvizy, bvizz;

  const int X = PointInputData::X;
  const int Y = PointInputData::Y;
  const int Z = PointInputData::Z;
  const int XX = PointInputData::XX;
  const int YX = PointInputData::YX;
  const int YY = PointInputData::YY;
  const int ZX = PointInputData::ZX;
  const int ZY = PointInputData::ZY;
  const int ZZ = PointInputData::ZZ;

  double grad[3], hess[6];
  if (need_gradient_) for (i=0; i<3; i++) grad[i] = 0.0;
  if (need_hessian_) for (i=0; i<6; i++) hess[i] = 0.0;

  CHECK_ALIGN(tmp);
  CHECK_ALIGN(grad[0]);

  if (need_gradient_) {
      for (i=0; i < ncontrib_bf_; i++) {
          int it = contrib_bf_[i];
          bvi = bs_values_[i];
          if (need_gradient_) {
              bvix = bsg_values_[i*3+X];
              bviy = bsg_values_[i*3+Y];
              bviz = bsg_values_[i*3+Z];
            }
          if (need_hessian_) {
              bvixx = bsh_values_[i*6+XX];
              bviyx = bsh_values_[i*6+YX];
              bviyy = bsh_values_[i*6+YY];
              bvizx = bsh_values_[i*6+ZX];
              bvizy = bsh_values_[i*6+ZY];
              bvizz = bsh_values_[i*6+ZZ];
            }
          int j3 = 0, j6 = 0;
          int itoff = (it*(it+1))>>1;
          int itjt;
          double t = 0.0;
          for (j=0; j < i; j++) {
              int jt = contrib_bf_[j];
              itjt = itoff+jt;

              densij = dmat[itjt];
              double bvj = bs_values_[j];

              CHECK_ALIGN(densij);
              CHECK_ALIGN(bvj);

              t += densij*bvi*bvj;

              double bvjx, bvjy, bvjz;
              if (need_gradient_) {
                  bvjx = bsg_values_[j3+X];
                  bvjy = bsg_values_[j3+Y];
                  bvjz = bsg_values_[j3+Z];
                  grad[X] += densij*(bvi*bvjx + bvj*bvix);
                  grad[Y] += densij*(bvi*bvjy + bvj*bviy);
                  grad[Z] += densij*(bvi*bvjz + bvj*bviz);
                  j3 += 3;
                }

              if (need_hessian_) {
                  double bvjxx = bsh_values_[j6+XX];
                  double bvjyx = bsh_values_[j6+YX];
                  double bvjyy = bsh_values_[j6+YY];
                  double bvjzx = bsh_values_[j6+ZX];
                  double bvjzy = bsh_values_[j6+ZY];
                  double bvjzz = bsh_values_[j6+ZZ];

                  hess[XX] += densij*(bvi*bvjxx+bvix*bvjx+bvjx*bvix+bvixx*bvj);
                  hess[YX] += densij*(bvi*bvjyx+bviy*bvjx+bvjy*bvix+bviyx*bvj);
                  hess[YY] += densij*(bvi*bvjyy+bviy*bvjy+bvjy*bviy+bviyy*bvj);
                  hess[ZX] += densij*(bvi*bvjzx+bviz*bvjx+bvjz*bvix+bvizx*bvj);
                  hess[ZY] += densij*(bvi*bvjzy+bviz*bvjy+bvjz*bviy+bvizy*bvj);
                  hess[ZZ] += densij*(bvi*bvjzz+bviz*bvjz+bvjz*bviz+bvizz*bvj);

                  j6 += 6;
                }
            }
          densij = dmat[itoff+it]*bvi;
          tmp += t + 0.5*densij*bvi;
          if (need_gradient_) {
              grad[X] += densij*bvix;
              grad[Y] += densij*bviy;
              grad[Z] += densij*bviz;
            }
          if (need_hessian_) {
              hess[XX] += densij*bvixx;
              hess[YX] += densij*bviyx;
              hess[YY] += densij*bviyy;
              hess[ZX] += densij*bvizx;
              hess[ZY] += densij*bvizy;
              hess[ZZ] += densij*bvizz;
            }
        }
    }
  else {
      for (i=0; i < ncontrib_bf_; i++) {
          int it = contrib_bf_[i];
          bvi = bs_values_[i];
          int itoff = (it*(it+1))>>1;
          int itjt;
          double t = 0.0;
          for (j=0; j < i; j++) {
              int jt = contrib_bf_[j];
              itjt = itoff+jt;

              densij = dmat[itjt];
              double bvj = bs_values_[j];

              CHECK_ALIGN(densij);
              CHECK_ALIGN(bvj);

              t += densij*bvi*bvj;
            }
          densij = dmat[itoff+it]*bvi;
          tmp += t + 0.5*densij*bvi;
        }
    }
  d.rho = 2.0 * tmp;

  if (need_gradient_) {
      for (i=0; i<3; i++) d.del_rho[i] = 2.0 * grad[i];
      d.gamma = dot(d.del_rho,d.del_rho);
    }

  if (need_hessian_) {
      for (i=0; i<6; i++) d.hes_rho[i] = 2.0 * hess[i];
      d.lap_rho = d.hes_rho[XX] + d.hes_rho[YY] + d.hes_rho[ZZ];
    }

  tim_exit("get_density");
}

#if COUNT_CONTRIBUTIONS
static int *contrib_array = 0;
#endif

double
DenIntegrator::do_point(int acenter, const SCVector3 &r,
                        const RefDenFunctional &func,
                        double weight, double multiplier,
                        double *nuclear_gradient,
                        double *f_gradient, double *w_gradient)
{
  tim_enter("do_point");

  int i,j,k;
  double w_mult = weight * multiplier;

  CHECK_ALIGN(w_mult);

  GaussianBasisSet *basis = wfn_->basis().pointer();

  // only consider those shells for which phi_i * (Max_j D_ij phi_j) > tol
  if (linear_scaling_ && use_dmat_bound_) {
      const ArrayExtentData &cs = extent_->contributing_shells(r[0],r[1],r[2]);
      ncontrib_ = 0;
      for (i=0; i<cs.size(); i++) {
          int ish = cs[i].shell;
          int contrib = 0;
          for (j=0; j<cs.size(); j++) {
              int jsh = cs[j].shell;
              int ijsh = (ish>jsh)?((ish*(ish+1))/2+jsh):((jsh*(jsh+1))/2+ish);
              if (cs[i].bound*cs[j].bound*dmat_bound_[ijsh] > accuracy_) {
                  contrib = 1;
                  break;
                }
            }
          if (contrib) {
              contrib_[ncontrib_++] = ish;
            }
        }
    }
  else if (linear_scaling_) {
      const ArrayExtentData &cs = extent_->contributing_shells(r[0],r[1],r[2]);
      ncontrib_ = cs.size();
      for (i=0; i<ncontrib_; i++) {
          contrib_[i] = cs[i].shell;
        }
    }
  else {
      ncontrib_ = nshell_;
      for (i=0; i<nshell_; i++) contrib_[i] = i;
    }
  if (ncontrib_ > nshell_) {
      cout << "DenIntegrator::do_point: ncontrib invalid" << endl;
      abort();
    }
#if COUNT_CONTRIBUTIONS
  contrib_array[ncontrib_]++;
#endif
  if (ncontrib_ == 0) { tim_exit("do_point"); return 0.0; }

  ncontrib_bf_ = 0;
  for (i=0; i<ncontrib_; i++) {
      int nbf = basis->shell(contrib_[i]).nfunction();
      int bf = basis->shell_to_function(contrib_[i]);
      for (j=0; j<nbf; j++, bf++) {
          contrib_bf_[ncontrib_bf_++] = bf;
        }
    }

  // compute the basis set values
  double *bsh = bsh_values_, *bsg = bsg_values_, *bsv = bs_values_;
  for (i=0; i<ncontrib_; i++) {
      basis->hessian_shell_values(r,contrib_[i],bsh,bsg,bsv);
      int shsize = basis->shell(contrib_[i]).nfunction();
      if (bsh) bsh += 6 * shsize;
      if (bsg) bsg += 3 * shsize;
      if (bsv) bsv += shsize;
    }

  // loop over basis functions adding contributions to the density
  PointInputData id(r);

  get_density(alpha_dmat_, id.a);
  if (spin_polarized_) {
      get_density(beta_dmat_, id.b);
      if (need_gradient_) {
          id.gamma_ab = id.a.del_rho[0]*id.b.del_rho[0]
                      + id.a.del_rho[1]*id.b.del_rho[1] 
                      + id.a.del_rho[2]*id.b.del_rho[2];
        }
    }
  id.compute_derived(spin_polarized_);

  PointOutputData od;
  if ( (id.a.rho + id.b.rho) > 1e2*DBL_EPSILON) {
      if (nuclear_gradient == 0) {
          func->point(id, od);
        }
      else {
          func->gradient(id, od, f_gradient, acenter, basis,
                         alpha_dmat_, (spin_polarized_?beta_dmat_:alpha_dmat_),
                         ncontrib_, contrib_, ncontrib_bf_, contrib_bf_,
                         bs_values_, bsg_values_, bsh_values_);
        }
    }
  else { tim_exit("do_point"); return id.a.rho + id.b.rho; }
  
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

          for (int j=0; j<ncontrib_bf_; j++) {
              int jt = contrib_bf_[j];
              double dfdra_phi_m = drhoa*bs_values_[j];
              double dfdga_phi_m = gradsa[0]*bsg_values_[j*3+0] +
                                   gradsa[1]*bsg_values_[j*3+1] +
                                   gradsa[2]*bsg_values_[j*3+2];
              double vamu = dfdra_phi_m + dfdga_phi_m, vbmu=0.0;
              double dfdrb_phi_m, dfdgb_phi_m;
              if (spin_polarized_) {
                  dfdrb_phi_m = drhob*bs_values_[j];
                  dfdgb_phi_m = gradsb[0]*bsg_values_[j*3+0] +
                                       gradsb[1]*bsg_values_[j*3+1] +
                                       gradsb[2]*bsg_values_[j*3+2];
                  vbmu = dfdrb_phi_m + dfdgb_phi_m;
                }

              int jtoff = (jt*(jt+1))>>1;

              for (int k=0; k <= j; k++) {
                  int kt = contrib_bf_[k];
                  int jtkt = jtoff + kt;

                  double dfdga_phi_n = gradsa[0]*bsg_values_[k*3+0] +
                                       gradsa[1]*bsg_values_[k*3+1] +
                                       gradsa[2]*bsg_values_[k*3+2];
                  alpha_vmat_[jtkt] += vamu * bs_values_[k] +
                                     dfdga_phi_n * bs_values_[j];
                  if (spin_polarized_) {
                      double dfdgb_phi_n = gradsb[0]*bsg_values_[k*3+0] +
                                           gradsb[1]*bsg_values_[k*3+1] +
                                           gradsb[2]*bsg_values_[k*3+2];
                      beta_vmat_[jtkt] += vbmu * bs_values_[k] +
                                          dfdgb_phi_n * bs_values_[j];
                    }
                }
            }
        }
      else {
          double drhoa = w_mult*od.df_drho_a;
          double drhob = w_mult*od.df_drho_b;
          for (int j=0; j<ncontrib_bf_; j++) {
              int jt = contrib_bf_[j];
              double bsj = bs_values_[j];
              double dfa_phi_m = drhoa * bsj;
              double dfb_phi_m = drhob * bsj;
              int jtoff = (jt*(jt+1))>>1;
              for (int k=0; k <= j; k++) {
                  int kt = contrib_bf_[k];
                  int jtkt = jtoff + kt;
                  double bsk = bs_values_[k];
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
          for (i=0; i<natom_*3; i++) {
              nuclear_gradient[i] += w_gradient[i] * od.energy * multiplier;
            }
        }
      // the contribution from (df/dx) w
      for (i=0; i<natom_*3; i++) {
          nuclear_gradient[i] += f_gradient[i] * w_mult;
        }
    }

  tim_exit("do_point");
  return id.a.rho + id.b.rho;
}

///////////////////////////////////////////////////////////////////////////
// IntegrationWeight

#define CLASSNAME IntegrationWeight
#define PARENTS public SavableState
#include <util/state/statei.h>
#include <util/class/classia.h>
void *
IntegrationWeight::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

IntegrationWeight::IntegrationWeight(StateIn& s):
  SavableState(s)
{
}

IntegrationWeight::IntegrationWeight()
{
}

IntegrationWeight::IntegrationWeight(const RefKeyVal& keyval)
{
}

IntegrationWeight::~IntegrationWeight()
{
}

void
IntegrationWeight::save_data_state(StateOut& s)
{
  cout << class_name() << ": cannot save state" << endl;
  abort();
}

void
IntegrationWeight::init(const RefMolecule &mol, double tolerance)
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
  RefMolecule molsav = mol_;
  RefMolecule dmol = new Molecule(*mol_.pointer());
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
//            cout << scprintf("%d,%d %12.10f %12.10f %12.10f",
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
      cout << "IntegrationWeight::test: failed on weight" << endl;
          cout << "sum_w = " << sum_weight << endl;
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
          cout << "iatom = " << i/3
               << " ixyx = " << i%3
               << " icenter = " << icenter << " point = " << point << endl;
          cout << scprintf("dw/dx bad: fd_val=%16.13f an_val=%16.13f err=%16.13f",
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

#define CLASSNAME BeckeIntegrationWeight
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public IntegrationWeight
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
BeckeIntegrationWeight::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = IntegrationWeight::_castdown(cd);
  return do_castdowns(casts,cd);
}

BeckeIntegrationWeight::BeckeIntegrationWeight(StateIn& s):
  SavableState(s),
  IntegrationWeight(s)
{
  bragg_radius = 0;
  a_mat = 0;
  oorab = 0;

  abort();
}

BeckeIntegrationWeight::BeckeIntegrationWeight()
{
  centers = 0;
  bragg_radius = 0;
  a_mat = 0;
  oorab = 0;
}

BeckeIntegrationWeight::BeckeIntegrationWeight(const RefKeyVal& keyval):
  IntegrationWeight(keyval)
{
  centers = 0;
  bragg_radius = 0;
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
  cout << ": cannot save state" << endl;
  abort();
}

void
BeckeIntegrationWeight::init(const RefMolecule &mol, double tolerance)
{
  done();
  IntegrationWeight::init(mol, tolerance);

  ncenters = mol->natom();

  double *bragg_radius = new double[ncenters];
  int icenter;
  for (icenter=0; icenter<ncenters; icenter++) {
      bragg_radius[icenter] = mol->atominfo()->bragg_radius(mol->Z(icenter));
    }
  
  centers = new SCVector3[ncenters];
  for (icenter=0; icenter<ncenters; icenter++) {
      centers[icenter].x() = mol->r(icenter,0);
      centers[icenter].y() = mol->r(icenter,1);
      centers[icenter].z() = mol->r(icenter,2);
    }

  a_mat = new double*[ncenters];
  a_mat[0] = new double[ncenters*ncenters];
  oorab = new double*[ncenters];
  oorab[0] = new double[ncenters*ncenters];

  for (icenter=0; icenter < ncenters; icenter++) {
      if (icenter) {
          a_mat[icenter] = &a_mat[icenter-1][ncenters];
          oorab[icenter] = &oorab[icenter-1][ncenters];
        }

      double bragg_radius_a = bragg_radius[icenter];
      
      for (int jcenter=0; jcenter < ncenters; jcenter++) {
          double chi=bragg_radius_a/bragg_radius[jcenter];
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
  delete[] bragg_radius;
  bragg_radius = 0;

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
}

double
BeckeIntegrationWeight::compute_p(int icenter, SCVector3&point)
{
  double ra = point.dist(centers[icenter]);
  double *ooraba = oorab[icenter];
  double *aa = a_mat[icenter];

  double p = 1.0;
  for (int jcenter=0; jcenter < ncenters; jcenter++) {
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
      for (int dcenter=0; dcenter<ncenters; dcenter++) {
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
  int icenter, jcenter;
  double p_sum=0.0, p_a=0.0;
    
  for (icenter=0; icenter<ncenters; icenter++) {
      double p_tmp = compute_p(icenter, point);
      if (icenter==acenter) p_a=p_tmp;
      p_sum += p_tmp;
    }
  double w_a = p_a/p_sum;

  if (w_gradient) {
      // w_gradient is computed using the formulae from
      // Johnson et al., JCP v. 98, p. 5612 (1993) (Appendix B)
      int i,j;
      for (i=0; i<ncenters*3; i++ ) w_gradient[i] = 0.0;
//      fd_w(acenter, point, w_gradient);  // imbn commented out for debug
//        cout << point << " ";
//        for (int i=0; i<ncenters*3; i++) {
//            cout << scprintf(" %10.6f", w_gradient[i]);
//          }
//        cout << endl;
//      return w_a;  // imbn commented out for debug
      for (int ccenter = 0; ccenter < ncenters; ccenter++) {
          // NB: for ccenter==acenter, use translational invariance
          // to get the corresponding component of the gradient
          if (ccenter != acenter) {
              SCVector3 grad_c_w_a;
              SCVector3 grad_c_p_a;
              compute_grad_p(ccenter, acenter, acenter, point, p_a, grad_c_p_a);
              for (i=0; i<3; i++) grad_c_w_a[i] = grad_c_p_a[i]/p_sum;
              for (int bcenter=0; bcenter<ncenters; bcenter++) {
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
          for (i=0; i<ncenters; i++) {
              if (i != acenter) {
                  w_gradient[acenter*3+j] -= w_gradient[i*3+j];
                }
            }
        }
    }

  return w_a;
}

///////////////////////////////////////////////////////////////////////////
// Murray93Integrator

// utility functions

static void
gauleg(double x1, double x2, double x[], double w[], int n)
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
      x[i-1]=xm-xl*z;
      x[n-i]=xm+xl*z;
      w[i-1]=2.0*xl/((1.0-z*z)*pp*pp);
      w[n-i]=w[i-1];
    }
}

#define CLASSNAME Murray93Integrator
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#define PARENTS public DenIntegrator
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
Murray93Integrator::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DenIntegrator::_castdown(cd);
  return do_castdowns(casts,cd);
}

Murray93Integrator::Murray93Integrator(StateIn& s):
  SavableState(s),
  DenIntegrator(s)
{
  abort();
}

Murray93Integrator::Murray93Integrator()
{
  nr_ = 64;
  ntheta_ = 16;
  nphi_ = 32;
  Ktheta_ = 5;
  weight_ = new BeckeIntegrationWeight;
}

Murray93Integrator::Murray93Integrator(const RefKeyVal& keyval):
  DenIntegrator(keyval)
{
  nr_ = keyval->intvalue("nr");
  if (nr_ == 0) nr_ = 64;
  ntheta_ = keyval->intvalue("ntheta");
  if (ntheta_ == 0) ntheta_ = 16;
  nphi_ = keyval->intvalue("nphi");
  if (nphi_ == 0) nphi_ = 2*ntheta_;
  Ktheta_ = keyval->intvalue("Ktheta");
  if (keyval->error() != KeyVal::OK)
      Ktheta_ = 5;
  weight_ = new BeckeIntegrationWeight;
}

Murray93Integrator::~Murray93Integrator()
{
}

void
Murray93Integrator::save_data_state(StateOut& s)
{
  cout << ": cannot save state" << endl;
  abort();
}


// taken from Mike Colvin 1997/07/21 and modified
void
Murray93Integrator::integrate(const RefDenFunctional &denfunc,
                              const RefSymmSCMatrix& densa,
                              const RefSymmSCMatrix& densb,
                              double *nuclear_gradient)
{
  tim_enter("integrate");

  init_integration(denfunc, densa, densb, nuclear_gradient);

#if COUNT_CONTRIBUTIONS
  delete[] contrib_array;
  contrib_array = new int[nshell_+1];
  memset(contrib_array, 0, sizeof(int)*(nshell_+1));
#endif

  RefMolecule mol = wavefunction()->molecule();
  weight_->init(mol, DBL_EPSILON);

  int ncenters=mol->natom();   // number of centers
  int icenter;                 // Loop index over centers

  int *nr = new int[ncenters];
  int ntheta = ntheta_;
  int nphi = nphi_;
  for (icenter=0; icenter<ncenters; icenter++) nr[icenter] = nr_;

  double *w_gradient = 0;
  double *f_gradient = 0;
  if (nuclear_gradient) {
      w_gradient = new double[ncenters*3];
      f_gradient = new double[ncenters*3];
    }

  SCVector3 *centers = new SCVector3[ncenters];
  for (icenter=0; icenter<ncenters; icenter++) {
      centers[icenter].x() = mol->r(icenter,0);
      centers[icenter].y() = mol->r(icenter,1);
      centers[icenter].z() = mol->r(icenter,2);
    }

  int ir, itheta, iphi;       // Loop indices for diff. integration dim
  int point_count;            // Counter for # integration points per center
  int point_count_total=0;    // Counter for # integration points

  SCVector3 point;            // sph coord of current integration point
  SCVector3 center;           // Cartesian position of center
  SCVector3 integration_point;

  double w,q,int_volume,sin_theta;
        
  // Determine maximium # grid points
  int nr_max=0;
  for (icenter=0;icenter<ncenters;icenter++) {
      if (nr[icenter]>nr_max) nr_max=nr[icenter];
    }

  // allocate memory for integration sub-expressions
  double *theta_quad_points =  new double[ntheta];
  double *theta_quad_weights = new double[ntheta];
  double *r_values = new double[nr_max];
  double *dr_dq = new double[nr_max];

  RefMessageGrp msg = MessageGrp::get_default_messagegrp();
  int nproc = msg->n();
  int me = msg->me();
  int parallel_counter = 0;

  double *bragg_radius = new double[ncenters];
  for (icenter=0; icenter<ncenters; icenter++) {
      bragg_radius[icenter] = mol->atominfo()->bragg_radius(mol->Z(icenter));
    }

  for (icenter=0; icenter < ncenters; icenter++) {
      if (! (parallel_counter++%nproc == me)) continue;

      point_count=0;
      center = centers[icenter];
      int r_done = 0;

      for (ir=0; ir < nr[icenter]; ir++) {
          // Mike Colvin's interpretation of Murray's radial grid
          q=(double)ir/(double)nr[icenter];
          double value=q/(1-q);
          double r = bragg_radius[icenter]*value*value;
          double dr_dq = 2.0*bragg_radius[icenter]*q*pow(1.0-q,-3.);
          double drdqr2 = dr_dq*r*r;
          point.r() = r;

          // Precalcute theta Gauss-Legendre abcissas and weights
          if (ir==0) {
              ntheta=1;
              nphi=1;
            }
          else {
              ntheta = (int) (r/bragg_radius[icenter] * Ktheta_*ntheta_);
              if (ntheta > ntheta_)
                  ntheta = ntheta_;
              if (ntheta < 6)
                  ntheta=6;
              nphi = 2*ntheta;
            }
          gauleg(0.0, M_PI, theta_quad_points, theta_quad_weights, ntheta);

          for (itheta=0; itheta < ntheta; itheta++) {
              point.theta() = theta_quad_points[itheta];
              sin_theta = sin(point.theta());
              // calculate integration volume 
              int_volume=drdqr2*sin_theta;
              for (iphi=0; iphi < nphi; iphi++) {
                  point.phi() = (double)iphi/(double)nphi * 2.0 * M_PI;
                  point.spherical_to_cartesian(integration_point);
                  integration_point += center;

                  // calculate weighting factor
                  w=weight_->w(icenter, integration_point, w_gradient);
                  //if (w_gradient) weight_->test(icenter, integration_point);
                    
                  // Update center point count
                  point_count++;

                  // Make contribution to Euler-Maclaurin formula
                  double multiplier = int_volume
                                    * theta_quad_weights[itheta]/nr[icenter]
                                    * 2.0 * M_PI / ((double)nphi);
                  if (do_point(icenter, integration_point, denfunc,
                               w, multiplier,
                               nuclear_gradient, f_gradient, w_gradient)
                      * multiplier < 1e2*DBL_EPSILON
                      && multiplier > 1e2*DBL_EPSILON) {
                      r_done=1;
                      break;
                    }
                }

              if (r_done)
                  break;
            }

          if (r_done)
              break;
        }
      point_count_total+=point_count;
    }

  msg->sum(point_count_total);
  done_integration();
  weight_->done();

     cout << node0 << indent
          << "Total integration points = " << point_count_total << endl;
    //cout << scprintf(" Value of integral = %16.14f", value()) << endl;

    delete[] f_gradient;
    delete[] w_gradient;
    delete[] theta_quad_points;
    delete[] bragg_radius;
    delete[] theta_quad_weights;
    delete[] r_values;
    delete[] dr_dq;
    delete[] nr;
    delete[] centers;

  tim_exit("integrate");

#if COUNT_CONTRIBUTIONS
  int tot = 0;
  double sav1 = 0.0;
  double sav2 = 0.0;
  for (int i=0; i<nshell_+1; i++) {
      cout << "contrib_array[" << setw(2) << i << "] = "
           << contrib_array[i] << endl;
      tot += contrib_array[i];
      sav1 += contrib_array[i]*i;
      sav2 += contrib_array[i]*i*i;
    }
  cout << "tot = " << tot << endl;
  cout << "sav1 = " << sav1/(tot*nshell_) << endl;
  cout << "sav2 = " << sav2/(tot*nshell_*nshell_) << endl;
#endif
}

void
Murray93Integrator::print(ostream &o) const
{
  o << node0 << indent << "Murray93Integrator Parameters:" << endl;
  o << incindent;
  o << node0 << indent << scprintf("nr     = %5d", nr_) << endl;
  o << node0 << indent << scprintf("ntheta = %5d", ntheta_) << endl;
  o << node0 << indent << scprintf("nphi   = %5d", nphi_) << endl;
  o << node0 << indent << scprintf("Ktheta = %5d", Ktheta_) << endl;
  o << decindent;
}

///////////////////////////////////////////////////
// Radial Integrator

#define CLASSNAME RadialIntegrator
#define PARENTS public SavableState
#include <util/state/statei.h>
#include <util/class/classia.h>
void *
RadialIntegrator::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

RadialIntegrator::RadialIntegrator(StateIn& s):
  SavableState(s)
{
}

RadialIntegrator::RadialIntegrator()
{
  set_nr(64);
}

RadialIntegrator::RadialIntegrator(const RefKeyVal& keyval)
{
  set_nr( keyval->intvalue("nr") );
  if (keyval->error() != KeyVal::OK)
      set_nr(64);
}

RadialIntegrator::~RadialIntegrator()
{
}

void
RadialIntegrator::save_data_state(StateOut& s)
{
  cout << class_name() << ": cannot save state" << endl;
  abort();
}

void
RadialIntegrator::set_nr(int i)
{
  nr_ = i;
}

int
RadialIntegrator::get_nr(void) const
{
  return nr_;
}

void
RadialIntegrator::print(ostream &o) const
{
  o << node0 << indent << class_name() << ":" << endl;
  o << incindent;
  o << node0 << indent << scprintf("nr       = %5d", get_nr()) << endl;
  o << decindent;
}

///////////////////////////////////////
//  AngularIntegrator

#define CLASSNAME AngularIntegrator
#define PARENTS public SavableState
#include <util/state/statei.h>
#include <util/class/classia.h>
void *
AngularIntegrator::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = SavableState::_castdown(cd);
  return do_castdowns(casts,cd);
}

AngularIntegrator::AngularIntegrator(StateIn& s):
  SavableState(s)
{
}

AngularIntegrator::AngularIntegrator()
{
}

AngularIntegrator::AngularIntegrator(const RefKeyVal& keyval)
{
}

AngularIntegrator::~AngularIntegrator()
{
}

void
AngularIntegrator::save_data_state(StateOut& s)
{
  cout << class_name() << ": cannot save state" << endl;
  abort();
}

///////////////////////////////////////
//  EulerMaclaurinRadialIntegrator

#define CLASSNAME EulerMaclaurinRadialIntegrator
#define PARENTS public RadialIntegrator
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
EulerMaclaurinRadialIntegrator::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = RadialIntegrator::_castdown(cd);
  return do_castdowns(casts,cd);
}

EulerMaclaurinRadialIntegrator::EulerMaclaurinRadialIntegrator(StateIn& s):
  SavableState(s),
  RadialIntegrator(s)
{
}

EulerMaclaurinRadialIntegrator::EulerMaclaurinRadialIntegrator()
{
}

EulerMaclaurinRadialIntegrator::EulerMaclaurinRadialIntegrator(const RefKeyVal& keyval):
RadialIntegrator(keyval)
{
}

EulerMaclaurinRadialIntegrator::~EulerMaclaurinRadialIntegrator()
{
}

void
EulerMaclaurinRadialIntegrator::save_data_state(StateOut& s)
{
  cout << class_name() << ": cannot save state" << endl;
  abort();
}

double
EulerMaclaurinRadialIntegrator::radial_value(int ir, int nr, double radii)
{
  double q = (double) (double)ir/(double)nr;
  double value = q/(1.-q);
  double r = radii*value*value;
  set_dr_dq( 2.*radii*q*pow(1.-q,-3.) );
  set_dr_dqr2( dr_dq_*r*r );
  return r;
}

double
EulerMaclaurinRadialIntegrator::radial_multiplier(int nr)
{
  double value = get_dr_dqr2();
  return value/((double) nr);
}

void
EulerMaclaurinRadialIntegrator::set_dr_dq(double i)
{
  dr_dq_ = i;
}

double
EulerMaclaurinRadialIntegrator::get_dr_dq(void) const
{
  return dr_dq_;
}

void
EulerMaclaurinRadialIntegrator::set_dr_dqr2(double i)
{
  dr_dqr2_ = i;
}

double
EulerMaclaurinRadialIntegrator::get_dr_dqr2(void) const
{
  return dr_dqr2_;
}

//////////////////////////////////////////////////////////////////////////
// LebedevAngularIntegrator
#define CLASSNAME LebedevAngularIntegrator
#define PARENTS public AngularIntegrator
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
LebedevAngularIntegrator::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = AngularIntegrator::_castdown(cd);
  return do_castdowns(casts,cd);
}

LebedevAngularIntegrator::LebedevAngularIntegrator(StateIn& s):
  SavableState(s),
  AngularIntegrator(s)
{
}

LebedevAngularIntegrator::LebedevAngularIntegrator()
{
  
  set_norder(59);
  set_npoints(1202);
  set_N1(13);
  set_N2(4);
  set_N3(16);
  set_point_count(0);
  
  int npoints = get_npoints();
  x_ = new double [npoints];
  y_ = new double [npoints];
  z_ = new double [npoints];
  lebedev_weights_ = new double [npoints];

  build_grid();
}

LebedevAngularIntegrator::LebedevAngularIntegrator(const RefKeyVal& keyval)
{
  set_norder( keyval->intvalue("norder") );
  if (keyval->error() != KeyVal::OK) set_norder(59);
  set_npoints( keyval->intvalue("npoints") );
  if (keyval->error() != KeyVal::OK) set_npoints(1202);
  set_N1( keyval->intvalue("N1") );
  if (keyval->error() != KeyVal::OK) set_N1(13);
  set_N2( keyval->intvalue("N2") );
  if (keyval->error() != KeyVal::OK) set_N2(4);
  set_N3( keyval->intvalue("N3") );
  if (keyval->error() != KeyVal::OK) set_N3(16);
  
  int npoints = get_npoints();
  x_ = new double [npoints];
  y_ = new double [npoints];
  z_ = new double [npoints];
  lebedev_weights_ = new double [npoints];
  set_point_count(0);

  if (npoints == 1202) {
      set_N1(13); set_N2(4); set_N3(16);
    }
  else if (npoints == 974) {
      set_N1(12); set_N2(4); set_N3(12);
    }
  else if (npoints == 770) {
      set_N1(10); set_N2(3); set_N3(9);
    }
  else if (npoints == 590) {
      set_N1(9); set_N2(3); set_N3(6);
    }
  else if (npoints == 302) {
      set_N1(6); set_N2(2); set_N3(2);
    }
  else if (npoints == 266) {
      set_N1(5); set_N2(1); set_N3(2);
    }
  else if (npoints == 194) {
      set_N1(4); set_N2(1); set_N3(1);
    }
  else if (npoints == 110) {
      set_N1(3); set_N2(1); set_N3(0);
    }
  else if (npoints == 50) {
      set_N1(1); set_N2(0); set_N3(0);
    }
  else if (npoints == 38) {
      set_N1(0); set_N2(1); set_N3(0);
    }
  else if (npoints == 6) {
      set_N1(0); set_N2(0); set_N3(0);
    }
  build_grid();
}

LebedevAngularIntegrator::~LebedevAngularIntegrator()
{
  delete [] x_;
  delete [] y_;
  delete [] z_;
  delete [] lebedev_weights_;
}

void
LebedevAngularIntegrator::save_data_state(StateOut& s)
{
  cout << class_name() << ": cannot save state" << endl;
  abort();
}

int
LebedevAngularIntegrator::get_norder(void) const
{
  return norder_;
}

void
LebedevAngularIntegrator::set_norder(int i)
{
  norder_ = i;
}

int
LebedevAngularIntegrator::get_npoints(void) const
{
  return npoints_;
}

void
LebedevAngularIntegrator::set_npoints(int i)
{
  npoints_ = i;
}

int
LebedevAngularIntegrator::get_N1(void) const
{
  return N1_;
}

void
LebedevAngularIntegrator::set_N1(int i)
{
  N1_ = i;
}

int
LebedevAngularIntegrator::get_N2(void) const
{
  return N2_;
}

void
LebedevAngularIntegrator::set_N2(int i)
{
  N2_ = i;
}

int
LebedevAngularIntegrator::get_N3(void) const
{
  return N3_;
}

void
LebedevAngularIntegrator::set_N3(int i)
{
  N3_ = i;
}

int
LebedevAngularIntegrator::get_point_count(void) const
{
  return point_count_;
}

void
LebedevAngularIntegrator::set_point_count(int i)
{
  point_count_ = i;
}

double
LebedevAngularIntegrator
::angular_point_cartesian(int iangular, SCVector3 &point,
                          SCVector3 &integration_point) const
{
  double r = point.r();
  integration_point.x() = r*x_[iangular];
  integration_point.y() = r*y_[iangular];
  integration_point.z() = r*z_[iangular];

  //point.theta() = acos(z_[iangular]);
  //if (sin(point.theta()) < DBL_EPSILON) point.phi() = 0.0;
  //else point.phi() = asin( y_[iangular]/sin(point.theta()) );
  return ( 4.0 * M_PI * lebedev_weights_[iangular] );
}

int
LebedevAngularIntegrator::num_angular_points(double r_value, int ir)
{
  if (ir == 0) return 1;
  else return get_npoints();
}

void
LebedevAngularIntegrator::angular_weights(void)
{
  return;
}

void
LebedevAngularIntegrator::build_grid(void)
{
  int norder, npoints, N1, N2, N3;
  int i, gen_point[6];

  norder = get_norder();
  npoints = get_npoints();
  N1 = get_N1();
  N2 = get_N2();
  N3 = get_N3();
  for (i=0; i<6; i++) gen_point[i] = 0;

  double A1, A2, A3, *l, *m, *B, *q, *r, *C, *u, *v, *w, *D;
  l = new double [N1]; m = new double [N1]; B = new double [N1];
  q = new double [N2]; r = new double [N2]; C = new double [N2];
  u = new double [N3]; v = new double [N3]; w = new double [N3]; D = new double [N3];
  A1 = A2 = A3 = 0.0;
  
  if (norder == 59 && npoints == 1202) {
      for (i=0; i<6; i++) gen_point[i] = 1;
      A1    = 0.110518923327e-3; A2    = 0.920523273809e-3; A3    = 0.913315978645e-3;
      
      l[0]  = 0.371263644966e-1; m[0]  = 0.998620681800e+0; B[0]  = 0.369042189802e-3;
      l[1]  = 0.914006041226e-1; m[1]  = 0.991610739722e+0; B[1]  = 0.560399092868e-3;
      l[2]  = 0.153107785247e+0; m[2]  = 0.976276606395e+0; B[2]  = 0.686529762928e-3;
      l[3]  = 0.218092889166e+0; m[3]  = 0.951247067481e+0; B[3]  = 0.772033855115e-3;
      l[4]  = 0.283987453220e+0; m[4]  = 0.915806886209e+0; B[4]  = 0.830154595889e-3;
      l[5]  = 0.349117760096e+0; m[5]  = 0.869616915182e+0; B[5]  = 0.868669255018e-3;
      l[6]  = 0.412143146144e+0; m[6]  = 0.812573722300e+0; B[6]  = 0.892707628585e-3;
      l[7]  = 0.471899362715e+0; m[7]  = 0.744729469632e+0; B[7]  = 0.906082023857e-3;
      l[8]  = 0.527314545284e+0; m[8]  = 0.666242253736e+0; B[8]  = 0.911977725494e-3;
      l[9]  = 0.620947533244e+0; m[9]  = 0.478380938077e+0; B[9]  = 0.912872013860e-3;
      l[10] = 0.656972271186e+0; m[10] = 0.369830866459e+0; B[10] = 0.913071493569e-3;
      l[11] = 0.684178830907e+0; m[11] = 0.252583955701e+0; B[11] = 0.915287378455e-3;
      l[12] = 0.701260433012e+0; m[12] = 0.128326186660e+0; B[12] = 0.918743627432e-3;
      
      q[0]  = 0.107238221548e+0; r[0]  = 0.994233354821e+0; C[0]  = 0.517697731297e-3;
      q[1]  = 0.258206895950e+0; r[1]  = 0.966089643296e+0; C[1]  = 0.733114368210e-3;
      q[2]  = 0.417275295531e+0; r[2]  = 0.908780131682e+0; C[2]  = 0.846323283638e-3;
      q[3]  = 0.570036691179e+0; r[3]  = 0.821619237061e+0; C[3]  = 0.903112269425e-3;

      u[0]  = 0.982798601826e+0; v[0]  = 0.177177402262e+0; w[0]  = 0.521063947701e-1; D[0]  = 0.648577845316e-3;
      u[1]  = 0.962424923033e+0; v[1]  = 0.247571646343e+0; w[1]  = 0.111564095716e+0; D[1]  = 0.743503091098e-3;
      u[2]  = 0.940200799413e+0; v[2]  = 0.335461628907e+0; w[2]  = 0.590588885324e-1; D[2]  = 0.799852789184e-3;
      u[3]  = 0.932082204014e+0; v[3]  = 0.317361524661e+0; w[3]  = 0.174655167758e+0; D[3]  = 0.810173148747e-3;
      u[4]  = 0.904367419939e+0; v[4]  = 0.409026842709e+0; w[4]  = 0.121723505110e+0; D[4]  = 0.848338957459e-3;
      u[5]  = 0.891240756007e+0; v[5]  = 0.385429115067e+0; w[5]  = 0.239027847938e+0; D[5]  = 0.855629925731e-3;
      u[6]  = 0.867643562846e+0; v[6]  = 0.493222118485e+0; w[6]  = 0.626625062415e-1; D[6]  = 0.880320867974e-3;
      u[7]  = 0.858197998604e+0; v[7]  = 0.478532067592e+0; w[7]  = 0.185750519455e+0; D[7]  = 0.881104818243e-3;
      u[8]  = 0.839675362405e+0; v[8]  = 0.450742259316e+0; w[8]  = 0.302946697353e+0; D[8]  = 0.885028234127e-3;
      u[9]  = 0.816528856402e+0; v[9]  = 0.563212302076e+0; w[9]  = 0.126777480068e+0; D[9]  = 0.902134229904e-3;
      u[10] = 0.801546937078e+0; v[10] = 0.543430356969e+0; w[10] = 0.249411216236e+0; D[10] = 0.901009167711e-3;
      u[11] = 0.777356306907e+0; v[11] = 0.512351848642e+0; w[11] = 0.364983226060e+0; D[11] = 0.902269293843e-3;
      u[12] = 0.766162121390e+0; v[12] = 0.639427963475e+0; w[12] = 0.642454922422e-1; D[12] = 0.915801617469e-3;
      u[13] = 0.755358414353e+0; v[13] = 0.626980550902e+0; w[13] = 0.190601822278e+0; D[13] = 0.913157800319e-3;
      u[14] = 0.734430575756e+0; v[14] = 0.603116169310e+0; w[14] = 0.311227594715e+0; D[14] = 0.910781357948e-3;
      u[15] = 0.704383718402e+0; v[15] = 0.569370249847e+0; w[15] = 0.423864478152e+0; D[15] = 0.910576025897e-3;
    }
  else if (norder == 53 && npoints == 974) {
      for (i=0; i<6; i++) { if(i!=1) gen_point[i] = 1; }
      A1 = 0.1438294190e-3; A2 = 0.0; A3 = 0.1125772288e-2;
      l[0]  = 0.4292963545e-1; m[0]  = 0.9981553450e+0; B[0]  = 0.4948029342e-3;
      l[1]  = 0.1051426854e+0; m[1]  = 0.9888832243e+0; B[1]  = 0.7357990108e-3;
      l[2]  = 0.1750024867e+0; m[2]  = 0.9688902204e+0; B[2]  = 0.8889132771e-3;
      l[3]  = 0.2477653379e+0; m[3]  = 0.9366027304e+0; B[3]  = 0.9888347838e-3;
      l[4]  = 0.3206567123e+0; m[4]  = 0.8912679426e+0; B[4]  = 0.1053299681e-2;
      l[5]  = 0.3916520749e+0; m[5]  = 0.8325967237e+0; B[5]  = 0.1092778807e-2;
      l[6]  = 0.4590825874e+0; m[6]  = 0.7605829053e+0; B[6]  = 0.1114389394e-2;
      l[7]  = 0.5214563888e+0; m[7]  = 0.6754009691e+0; B[7]  = 0.1123724788e-2;
      l[8]  = 0.6253170244e+0; m[8]  = 0.4668589056e+0; B[8]  = 0.1125239325e-2;
      l[9]  = 0.6637926744e+0; m[9]  = 0.3446136542e+0; B[9]  = 0.1126153271e-2;
      l[10] = 0.6910410398e+0; m[10] = 0.2119541518e+0; B[10] = 0.1130286931e-2;
      l[11] = 0.7052907007e+0; m[11] = 0.7162440144e-1; B[11] = 0.1134986534e-2;

      q[0]  = 0.1236686762e+0; r[0]  = 0.9923235654e+0; C[0]  = 0.6823367927e-3;
      q[1]  = 0.2940777114e+0; r[1]  = 0.9557815124e+0; C[1]  = 0.9454158160e-3;
      q[2]  = 0.4697753849e+0; r[2]  = 0.8827859807e+0; C[2]  = 0.1074429975e-2;
      q[3]  = 0.6334563241e+0; r[3]  = 0.7737784472e+0; C[3]  = 0.1129300086e-2;

      u[0]  = 0.5974048614e-1; v[0]  = 0.2029128752e+0; w[0]  = 0.9773727228e+0; D[0]  = 0.8436884500e-3;
      u[1]  = 0.1375760408e+0; v[1]  = 0.4602621942e+0; w[1]  = 0.8770584618e+0; D[1]  = 0.1075255720e-2;
      u[2]  = 0.3391016526e+0; v[2]  = 0.5030673999e+0; w[2]  = 0.7949422999e+0; D[2]  = 0.1108577236e-2;
      u[3]  = 0.1271675191e+0; v[3]  = 0.2817606422e+0; w[3]  = 0.9510201693e+0; D[3]  = 0.9566475323e-3;
      u[4]  = 0.2693120740e+0; v[4]  = 0.4331561291e+0; w[4]  = 0.8601434616e+0; D[4]  = 0.1080663250e-2;
      u[5]  = 0.1419786452e+0; v[5]  = 0.6256167358e+0; w[5]  = 0.7671021862e+0; D[5]  = 0.1126797131e-2;
      u[6]  = 0.6709284600e-1; v[6]  = 0.3798395216e+0; w[6]  = 0.9226161107e+0; D[6]  = 0.1022568715e-2;
      u[7]  = 0.7057738183e-1; v[7]  = 0.5517505421e+0; w[7]  = 0.8310175524e+0; D[7]  = 0.1108960267e-2;
      u[8]  = 0.2783888477e+0; v[8]  = 0.6029619156e+0; w[8]  = 0.7476206108e+0; D[8]  = 0.1122790653e-2;
      u[9]  = 0.1979578938e+0; v[9]  = 0.3589606329e+0; w[9]  = 0.9121183784e+0; D[9]  = 0.1032401847e-2;
      u[10] = 0.2087307061e+0; v[10] = 0.5348666438e+0; w[10] = 0.8187485362e+0; D[10] = 0.1107249382e-2;
      u[11] = 0.4055122137e+0; v[11] = 0.5674997546e+0; w[11] = 0.7165918454e+0; D[11] = 0.1121780048e-2;
    }
  else if (norder == 47 && npoints == 770) {
      for (i=0; i<6; i++) gen_point[i] = 1;
      A1 = 0.2192942090e-3;  A2 = 0.1436433617e-2;  A3 = 0.1421940344e-2;
      l[0]  = 0.5087204410e-1;  m[0]  = 0.9974086776e+0;  B[0]  = 0.6798123510e-3;
      l[1]  = 0.1228198790e+0;  m[1]  = 0.9847997535e+0;  B[1]  = 0.9913184235e-3;
      l[2]  = 0.2026890814e+0;  m[2]  = 0.9580366759e+0;  B[2]  = 0.1180207833e-2;
      l[3]  = 0.2847745156e+0;  m[3]  = 0.9153179504e+0;  B[3]  = 0.1296599602e-2;
      l[4]  = 0.3656719078e+0;  m[4]  = 0.8559019286e+0;  B[4]  = 0.1365871427e-2;
      l[5]  = 0.4428264886e+0;  m[5]  = 0.7796213195e+0;  B[5]  = 0.1402988604e-2;
      l[6]  = 0.5140619627e+0;  m[6]  = 0.6866444472e+0;  B[6]  = 0.1418645563e-2;
      l[7]  = 0.6306401219e+0;  m[7]  = 0.4523119203e+0;  B[7]  = 0.1421376741e-2;
      l[8]  = 0.6716883332e+0;  m[8]  = 0.3125213050e+0;  B[8]  = 0.1423996475e-2;
      l[9]  = 0.6979792685e+0;  m[9]  = 0.1601558034e+0;  B[9]  = 0.1431554042e-2;

      q[0] = 0.1446865674e+0; r[0] = 0.9894775374e+0; C[0] = 0.9254401499e-3;
      q[1] = 0.3390263475e+0; r[1] = 0.9407768787e+0; C[1] = 0.1250239995e-2;
      q[2] = 0.5335804651e+0; r[2] = 0.8457493051e+0; C[2] = 0.1394365843e-2;

      u[0] = 0.6944024393e-1; v[0] = 0.2355187894e+0; w[0] = 0.9693858634e+0; D[0] = 0.1127089094e-2;
      u[1] = 0.2269004109e+0; v[1] = 0.4102182474e+0; w[1] = 0.8833103605e+0; D[1] = 0.1345753761e-2;
      u[2] = 0.8025574608e-1; v[2] = 0.6214302417e+0; w[2] = 0.7793481057e+0; D[2] = 0.1424957283e-2;
      u[3] = 0.1467999527e+0; v[3] = 0.3245284345e+0; w[3] = 0.9344148270e+0; D[3] = 0.1261523341e-2;
      u[4] = 0.1571507769e+0; v[4] = 0.5224482189e+0; w[4] = 0.8380641334e+0; D[4] = 0.1392547106e-2;
      u[5] = 0.2365702993e+0; v[5] = 0.6017546634e+0; w[5] = 0.7628406246e+0; D[5] = 0.1418761677e-2;
      u[6] = 0.7714815866e-1; v[6] = 0.4346575516e+0; w[6] = 0.8972853361e+0; D[6] = 0.1338366684e-2;
      u[7] = 0.3062936666e+0; v[7] = 0.4908826589e+0; w[7] = 0.8156092232e+0; D[7] = 0.1393700862e-2;
      u[8] = 0.3822477379e+0; v[8] = 0.5648768149e+0; w[8] = 0.7313007936e+0; D[8] = 0.1415914757e-2;
    }
  else if (norder == 41 && npoints == 590) {
      for (i=0; i<6; i++) { if(i!=1) gen_point[i] = 1; }
      A1 = 0.3095121295e-3; A2 = 0.0; A3 = 0.1852379698e-2;
      l[0] = 0.6095034115e-1; m[0] = 0.9962781297e+0; B[0] = 0.9764331164e-3;
      l[1] = 0.1459036449e+0; m[1] = 0.9784805837e+0; B[1] = 0.1384737234e-2;
      l[2] = 0.2384736701e+0; m[2] = 0.9414141582e+0; B[2] = 0.1617210647e-2;
      l[3] = 0.3317920736e+0; m[3] = 0.8830787279e+0; B[3] = 0.1749564657e-2;
      l[4] = 0.4215761784e+0; m[4] = 0.8028368773e+0; B[4] = 0.1818471778e-2;
      l[5] = 0.5044419707e+0; m[5] = 0.7007685753e+0; B[5] = 0.1846715956e-2;
      l[6] = 0.6372546939e+0; m[6] = 0.4333738687e+0; B[6] = 0.1852028828e-2;
      l[7] = 0.6807744066e+0; m[7] = 0.2703560883e+0; B[7] = 0.1858812585e-2;
      l[8] = 0.7040954938e+0; m[8] = 0.9219040707e-1; B[8] = 0.1871790639e-2;

      q[0] = 0.1724782009e+0; r[0] = 0.9850133350e+0; C[0] = 0.1300321685e-2;
      q[1] = 0.3964755348e+0; r[1] = 0.9180452877e+0; C[1] = 0.1705153996e-2;
      q[2] = 0.6116843442e+0; r[2] = 0.7911019296e+0; C[2] = 0.1857161196e-2;

      u[0] = 0.8213021581e-1; v[0] = 0.2778673190e+0; w[0] = 0.9571020743e+0; D[0] = 0.1555213603e-2;
      u[1] = 0.8999205842e-1; v[1] = 0.5033564271e+0; w[1] = 0.8593798558e+0; D[1] = 0.1802239128e-2;
      u[2] = 0.1816640840e+0; v[2] = 0.5984126497e+0; w[2] = 0.7803207424e+0; D[2] = 0.1849830560e-2;
      u[3] = 0.1720795225e+0; v[3] = 0.3791035407e+0; w[3] = 0.9092134750e+0; D[3] = 0.1713904507e-2;
      u[4] = 0.2634716655e+0; v[4] = 0.4742392842e+0; w[4] = 0.8400474883e+0; D[4] = 0.1802658934e-2;
      u[5] = 0.3518280927e+0; v[5] = 0.5610263808e+0; w[5] = 0.7493106119e+0; D[5] = 0.1842866472e-2;
    }
  else if (norder == 29 && npoints == 302) {
      for (i=0; i<6; i++) { if(i!=1) gen_point[i] = 1; }
      A1 = 0.854591172878e-3; A3 = 0.359911928502e-2;
      l[0] = 0.701176641609; m[0] = 0.129238672710; B[0] = 0.365004580768e-2;
      l[1] = 0.656632941022; m[1] = 0.371034178385; B[1] = 0.360482260142e-2;
      l[2] = 0.472905413258; m[2] = 0.743452042987; B[2] = 0.357672966173e-2;
      l[3] = 0.351564034558; m[3] = 0.867643624544; B[3] = 0.344978842429e-2;
      l[4] = 0.221964523631; m[4] = 0.949454317226; B[4] = 0.310895312238e-2;
      l[5] = 0.0961830852303; m[5] = 0.990705621379; B[5] = 0.235210141366e-2;

      q[0] = 0.571895589188; r[0] = 0.820326419828; C[0] = 0.360082093222e-2;
      q[1] = 0.264415288706; r[1] = 0.964408914879; C[1] = 0.298234496317e-2;

      u[0] = 0.251003475177; v[0] = 0.800072749407; w[0] = 0.544867737258; D[0] = 0.357154055427e-2;
      u[1] = 0.902442529533; v[1] = 0.412772408317; w[1] = 0.123354853258; D[1] = 0.339231220501e-2;
    }
  else if (norder == 27 && npoints == 266) {
      for (i=0; i<6; i++) gen_point[i] = 1;
      A1 = -1.31376912733e-3; A2 = -2.52272870489e-3; A3 = 4.18685388170e-3;
      l[0] = 0.703937339159; m[0] = 0.0945750764036; B[0] = 5.31516797782e-3;
      l[1] = 0.662033866370; m[1] = 0.351315128565;  B[1] = 4.25613135143e-3;
      l[2] = 0.464744872642; m[2] = 0.753673939251;  B[2] = 4.11248239441e-3;
      l[3] = 0.327742065497; m[3] = 0.886098344997;  B[3] = 3.59558489976e-3;
      l[4] = 0.101252624857; m[4] = 0.9896948074629; B[4] = 4.04714237709e-3;

      q[0] = 0.525731112119; r[0] = 0.850650808352; C[0] = 4.22958270065e-3;

      u[0] = 0.819343388819; v[0] = 0.524493924092; w[0] = 0.231479015871; D[0] = 4.07146759383e-3;
      u[1] = 0.939227929750; v[1] = 0.323348454269; w[1] = 0.115311196765; D[1] = 4.08091422578e-3;
    }
  else if (norder == 23 && npoints == 194) {
      for (i=0; i<6; i++) gen_point[i] = 1;
      //A1 = pow(2.,7.)*73./5242545.; A2 = pow(2.,14.)*1663./( pow(3.,3.)*pow(11.,2.)*1458821.);
      //A3 = pow(3.,10.)*1599797./( pow(2.,9.)*pow(173.,2.)*pow(13.,2.)*6545.);
      A1 = 0.178234044724e-2; A2 = 0.571690594998e-2; A3 = 0.557338317884e-2;
      l[0] = 0.444693317871; m[0] = 0.777493219315; B[0] = 0.551877146727e-2;
      l[1] = 0.289246562758; m[1] = 0.912509096867; B[1] = 0.515823771181e-2;
      l[2] = 0.671297344270; m[2] = 0.314196994183; B[2] = 0.560870408259e-2;
      l[3] = 0.129933544765; m[3] = 0.982972302707; B[3] = 0.410677702817e-2;

      q[0] = 0.345770219761; r[0] = 0.938319218128;
      // C[0] = pow(38.,4.)/(pow(33.,2.)*pow(7.,3.)*1105.);
      C[0] = 0.505184606462e-2;
      
      u[0] = 0.159041710538; v[0] = 0.836036015482; w[0] = 0.525118572443; D[0] = 0.553024891623e-2;
    }
  else if (norder == 17 && npoints == 110) {
      gen_point[0] = 1; gen_point[1] = 0; gen_point[2] = 1;
      gen_point[3] = 1; gen_point[4] = 1; gen_point[5] = 0;
      A1 = 0.382827049494e-2;
      A2 = 0.0;
      A3 = 0.98550016044e-2;
      l[0] = 0.185115635345; m[0] = 0.965124035087; B[0] = 0.844068048232e-2;
      l[1] = 0.383386152638; m[1] = 0.840255982384; B[1] = 0.959547133607e-2;
      l[2] = 0.690421048382; m[2] = 0.238807866929; B[2] = 0.994281489118e-2;

      q[0] = 0.478369028812; r[0] = 0.878158910604; C[0] = 0.969499636166e-2;
      // C[0] = 4.*pow(17.,3.)/2027025.;
    }
  else if (norder == 11 && npoints == 50) {
      for (i=0; i<4; i++) gen_point[i] = 1;
      A1 = 4./315.; A2 = 64./2835.; A3 = 27./1280.;
      m[0] = 3./sqrt(11.); l[0] = 1./sqrt(2.) * sqrt( (1.-m[0]*m[0]) ); B[0] = pow(11., 4.)/725760.;
    }
  else if (norder == 9 && npoints == 38) {
      gen_point[0] = 1; gen_point[1] = 0; gen_point[2] = 1;
      gen_point[3] = 0; gen_point[4] = 1; gen_point[5] = 0;
      A1 = 1./105.; A2 = 0.0; A3 = 9./280.;
      q[0] = 0.459700843381; r[0] = 0.888073833977; C[0] = 1./35.;
    }
  else if (npoints == 6) {
      A1 = 1./6.;
      gen_point[0] = 1;
    }
  else {
      cout << class_name() << ": (norder, npoints) = " << norder <<"," << npoints
           << " is not a valid Lebedev angular grid." << endl;
      abort();
    }

  double norm = A1*6.0 + A2*12.0 + A3*8.0;
  int j;

/*
  // Norm of the Lebedev Weights
  for (j=0; j<N1; j++) norm += 24.0*B[j];
  for (j=0; j<N2; j++) norm += 24.0*C[j];
  for (j=0; j<N3; j++) norm += 48.0*D[j];
  cout << "norm = " << scprintf("%30.18f",norm) << endl;
  norm /= norm;
  //A1 *= norm; A2 *= norm; A3 *= norm;
  //for (j=0; j<N1; j++) B[j] *= norm;
  //for (j=0; j<N2; j++) C[j] *= norm;
  //for (j=0; j<N3; j++) D[j] *= norm;

  // Check that points are on a unit sphere
  for (j=0; j<N1; j++) {
      norm = 2.0*l[j]*l[j] + m[j]*m[j];
      cout << "Bk points norm[" << j << "] = " << scprintf("%20.15f",norm) << endl;
      //l[j] = 1./sqrt(2.) * sqrt(1.-m[j]*m[j]);
    }
  for (j=0; j<N2; j++) {
      norm = q[j]*q[j] + r[j]*r[j];
      cout << "Ck points norm[" << j << "] = " << scprintf("%20.15f",norm) << endl;
      //r[j] = sqrt(1.-q[j]*q[j]);
    }
  for (j=0; j<N3; j++) {
      norm = u[j]*u[j] + v[j]*v[j] + w[j]*w[j];
      cout << "Dk points norm[" << j << "] = " << scprintf("%20.15f",norm) << endl;
      //u[j] *= norm; v[j] *= norm; w[j] *= norm;
    }
*/  
  double *zero_array;
  zero_array = new double[N2];
  for (i=0; i<N2; i++) zero_array[i] = 0.0;

  double zero = 0.0;
  double one = 1.0;
  double sqrt2 = 1./sqrt(2.);
  double sqrt3 = 1./sqrt(3.);
  
  // call generate_points rouinte
  if (gen_point[0]) generate_points(&A1, 1, 3, &one,   &zero,  &zero);
  if (gen_point[1]) generate_points(&A2, 1, 3, &sqrt2, &sqrt2, &zero);
  if (gen_point[2]) generate_points(&A3, 1, 1, &sqrt3, &sqrt3, &sqrt3);
  if (gen_point[3]) generate_points(B,  N1, 3, l,      m,      l);
  if (gen_point[4]) generate_points(C,  N2, 6, q,      r,      zero_array); 
  if (gen_point[5]) generate_points(D,  N3, 6, u,      v,      w);
  
  delete [] zero_array;
  cout << " Total number of points in grid = " << get_point_count() << endl;    
}

void
LebedevAngularIntegrator::generate_points(double weights[], int N, int nsets,
                                          double  u[], double v[], double w[])
{
  int i,j,k, nsigns;
  double *p1, *p2, *p3;
  double *tmp_array;

  tmp_array = new double[3];
  
  for (i=0; i<nsets; i++) {
      switch(i) {
      case 0: p1 = u; p2 = v; p3 = w;
          break;
      case 1: p1 = w; p2 = u; p3 = v;
          break;
      case 2: p1 = v; p2 = w; p3 = u;
          break;
      case 3: p1 = v; p2 = u; p3 = w;
          break;
      case 4: p1 = u; p2 = w; p3 = v;
          break;
      case 5: p1 = w; p2 = v; p3 = u;
          break;
      default: cout << class_name() << ": i = " << i << " is not a valid option." << endl;
          abort();
        }
      for (j=0; j<N; j++) {
          tmp_array[0] = p1[j]; tmp_array[1] = p2[j]; tmp_array[2] = p3[j];
          expand(tmp_array,0, weights[j]);              
        }
    }
  delete [] tmp_array;
}

void
LebedevAngularIntegrator::expand(double array[], int offset, double weight)
{
  if (offset > 2) {
      int point_count = get_point_count();
      x_[point_count] = array[0];
      y_[point_count] = array[1];
      z_[point_count] = array[2];
      lebedev_weights_[point_count] = weight;
      //cout << "(" << scprintf("%lf,%lf,%lf) w = %lf",
      //     x_[point_count], y_[point_count], z_[point_count], lebedev_weights_[point_count])
      //     << endl;
      set_point_count(point_count+1);
      return;
    }
  expand(array, offset+1, weight);

  if (array[offset] == 0.) return;
  else array[offset] = -array[offset];
  double *copy_array;
  copy_array = new double[3];
  memcpy(copy_array, array, sizeof(double)*3);
  expand(copy_array, offset+1, weight);
  
  delete [] copy_array;
}

void
LebedevAngularIntegrator::print(ostream &o) const
{
  o << node0 << indent << scprintf("norder   = %5d", get_norder()) << endl;
  o << node0 << indent << scprintf("nangular = %5d", get_npoints()) << endl;
}

/////////////////////////////////
//  GaussLegendreAngularIntegrator

#define CLASSNAME GaussLegendreAngularIntegrator
#define PARENTS public AngularIntegrator
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
GaussLegendreAngularIntegrator::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = AngularIntegrator::_castdown(cd);
  return do_castdowns(casts,cd);
}

GaussLegendreAngularIntegrator::GaussLegendreAngularIntegrator(StateIn& s):
  SavableState(s),
  AngularIntegrator(s)
{
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

GaussLegendreAngularIntegrator::GaussLegendreAngularIntegrator(const RefKeyVal& keyval)
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
  cout << class_name() << ": cannot save state" << endl;
  abort();
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
  return get_ntheta_r()*get_nphi_r();
}

void
GaussLegendreAngularIntegrator::angular_weights(void)
{
  gauleg(0.0, M_PI, get_ntheta_r());
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
::angular_point_cartesian(int iangular, SCVector3 &point,
                          SCVector3 &integration_point) const
{
  int itheta, iphi, nphi_r;

  nphi_r = get_nphi_r();
  itheta = iangular/nphi_r;
  iphi = iangular - itheta*nphi_r;
  point.theta() = theta_quad_points_[itheta];
  point.phi() = (double) iphi/ (double) nphi_r * 2.0 * M_PI;
  point.spherical_to_cartesian(integration_point);
  return ( sin_theta(point)*theta_quad_weights_[itheta]*2.0*M_PI/(double)nphi_r );
}

void
GaussLegendreAngularIntegrator::print(ostream &o) const
{
  o << node0 << indent << scprintf("ntheta   = %5d", get_ntheta()) << endl;
  o << node0 << indent << scprintf("nphi     = %5d", get_nphi()) << endl;
  o << node0 << indent << scprintf("Ktheta   = %5d", get_Ktheta()) << endl;
}

//////////////////////////////////////////////
//  RadialAngularIntegrator

#define CLASSNAME RadialAngularIntegrator
#define PARENTS public DenIntegrator
#define HAVE_KEYVAL_CTOR
#define HAVE_STATEIN_CTOR
#include <util/state/statei.h>
#include <util/class/classi.h>
void *
RadialAngularIntegrator::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DenIntegrator::_castdown(cd);
  return do_castdowns(casts,cd);
}

RadialAngularIntegrator::RadialAngularIntegrator(StateIn& s):
  SavableState(s),
  DenIntegrator(s)
{
  abort();
}

RadialAngularIntegrator::RadialAngularIntegrator()
{
  RadInt_ = new EulerMaclaurinRadialIntegrator;
  AngInt_ = new GaussLegendreAngularIntegrator;
  weight_ = new BeckeIntegrationWeight;
}

RadialAngularIntegrator::RadialAngularIntegrator(const RefKeyVal& keyval):
  DenIntegrator(keyval)
{
  RadInt_ = keyval->describedclassvalue("RadInt");
  if (RadInt_.null()) RadInt_ = new EulerMaclaurinRadialIntegrator(keyval);
  AngInt_ = keyval->describedclassvalue("AngInt");
  if (AngInt_.null()) AngInt_ = new GaussLegendreAngularIntegrator(keyval);
  weight_ = keyval->describedclassvalue("weight");
  if (weight_.null()) weight_ = new BeckeIntegrationWeight;

}

RadialAngularIntegrator::~RadialAngularIntegrator()
{
}

void
RadialAngularIntegrator::save_data_state(StateOut& s)
{
  cout << ": cannot save state" << endl;
  abort();
}

void
RadialAngularIntegrator::integrate(const RefDenFunctional &denfunc,
                              const RefSymmSCMatrix& densa,
                              const RefSymmSCMatrix& densb,
                              double *nuclear_gradient)
{
  tim_enter("integrate");

  init_integration(denfunc, densa, densb, nuclear_gradient);

#if COUNT_CONTRIBUTIONS
  delete[] contrib_array;
  contrib_array = new int[nshell_+1];
  memset(contrib_array, 0, sizeof(int)*(nshell_+1));
#endif

  RefMolecule mol = wavefunction()->molecule();
  weight_->init(mol, DBL_EPSILON);

  int ncenters=mol->natom();   // number of centers
  int icenter;                 // Loop index over centers

  int *nr = new int[ncenters];
  int nangular;
  
  for (icenter=0; icenter<ncenters; icenter++) nr[icenter] = RadInt_->get_nr();

  double *w_gradient = 0;
  double *f_gradient = 0;
  if (nuclear_gradient) {
      w_gradient = new double[ncenters*3];
      f_gradient = new double[ncenters*3];
    }

  SCVector3 *centers = new SCVector3[ncenters];
  for (icenter=0; icenter<ncenters; icenter++) {
      centers[icenter].x() = mol->r(icenter,0);
      centers[icenter].y() = mol->r(icenter,1);
      centers[icenter].z() = mol->r(icenter,2);
    }

  int ir, iangular;           // Loop indices for diff. integration dim
  int point_count;            // Counter for # integration points per center
  int point_count_total=0;    // Counter for # integration points

  SCVector3 point;            // sph coord of current integration point
  SCVector3 center;           // Cartesian position of center
  SCVector3 integration_point;

  double w,q,int_volume,radial_multiplier,angular_multiplier,dr_dqr2;
        
  // Determine maximium # grid points
  int nr_max=0;
  for (icenter=0;icenter<ncenters;icenter++) {
      if (nr[icenter]>nr_max) nr_max=nr[icenter];
    }

  RefMessageGrp msg = MessageGrp::get_default_messagegrp();
  int nproc = msg->n();
  int me = msg->me();
  int parallel_counter = 0;

  double *bragg_radius = new double[ncenters];
  for (icenter=0; icenter<ncenters; icenter++) {
      bragg_radius[icenter] = mol->atominfo()->bragg_radius(mol->Z(icenter));
    }

  for (icenter=0; icenter < ncenters; icenter++) {
      if (! (parallel_counter++%nproc == me)) continue;

      point_count=0;
      center = centers[icenter];
      int r_done = 0;
      for (ir=0; ir < nr[icenter]; ir++) {
//      for (ir=nr[icenter]-1; ir >= 0; ir--) {
          double r = RadInt_->radial_value(ir, nr[icenter], bragg_radius[icenter]);
          point.r() = r;
          //dr_dqr2 = RadInt_->get_dr_dqr2();
          radial_multiplier = RadInt_->radial_multiplier(nr[icenter]);
          nangular = AngInt_->num_angular_points(r/bragg_radius[icenter],ir);
          AngInt_->angular_weights();
          //double radial_int_volume = RadInt_->get_dr_dqr2();
          for (iangular=0; iangular<nangular; iangular++) {
              angular_multiplier =
                AngInt_->angular_point_cartesian(iangular, point, integration_point);
              integration_point += center;
              w=weight_->w(icenter, integration_point, w_gradient);
              //if (w_gradient) weight_->test(icenter, integration_point);
              point_count++;
              double multiplier = angular_multiplier * radial_multiplier;
              //double angular_int_volume = AngInt_->sin_theta(point);
              //double int_volume = radial_int_volume * angular_int_volume;
              if (do_point(icenter, integration_point, denfunc,
                           w, multiplier,
                           nuclear_gradient, f_gradient, w_gradient)
                  * multiplier < 1e2*DBL_EPSILON
                  && multiplier > 1e2*DBL_EPSILON) {
                  r_done=1;
                  // break;
                }
            }
          //if (r_done) 
              // break;
        }
      point_count_total+=point_count;
    }

  msg->sum(point_count_total);
  done_integration();
  weight_->done();

     cout << node0 << indent
          << "Total integration points = " << point_count_total << endl;
    //cout << scprintf(" Value of integral = %16.14f", value()) << endl;

    delete[] f_gradient;
    delete[] w_gradient;
    delete[] bragg_radius;
    delete[] nr;
    delete[] centers;

  tim_exit("integrate");

#if COUNT_CONTRIBUTIONS
  int tot = 0;
  double sav1 = 0.0;
  double sav2 = 0.0;
  for (int i=0; i<nshell_+1; i++) {
      cout << "contrib_array[" << setw(2) << i << "] = "
           << contrib_array[i] << endl;
      tot += contrib_array[i];
      sav1 += contrib_array[i]*i;
      sav2 += contrib_array[i]*i*i;
    }
  cout << "tot = " << tot << endl;
  cout << "sav1 = " << sav1/(tot*nshell_) << endl;
  cout << "sav2 = " << sav2/(tot*nshell_*nshell_) << endl;
#endif
}

void
RadialAngularIntegrator::print(ostream &o) const
{
  o << node0 << indent << class_name() << " Parameters:" << endl;
  o << incindent;
  RadInt_->print(o);
  AngInt_->print(o);
  o << decindent;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

