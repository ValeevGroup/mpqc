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
              if (cs[i].bound*cs[j].bound*dmat_bound_[ijsh] > DBL_EPSILON) {
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
                      * int_volume < 1e2*DBL_EPSILON
                      && int_volume > 1e2*DBL_EPSILON) {
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
RadialIntegrator::get_nr(void)
{
  return nr_;
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
EulerMaclaurinRadialIntegrator::radial_multiplier(double value, int nr)
{
  return value/((double) nr);
}

void
EulerMaclaurinRadialIntegrator::set_dr_dq(double i)
{
  dr_dq_ = i;
}

double
EulerMaclaurinRadialIntegrator::get_dr_dq(void)
{
  return dr_dq_;
}

void
EulerMaclaurinRadialIntegrator::set_dr_dqr2(double i)
{
  dr_dqr2_ = i;
}

double
EulerMaclaurinRadialIntegrator::get_dr_dqr2(void)
{
  return dr_dqr2_;
}

/////////////////////////////////
//  GaussLegendreAngularIntegrator

#define CLASSNAME GaussLegendreAngularIntegrator
#define PARENTS public AngularIntegrator
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
GaussLegendreAngularIntegrator::get_ntheta(void)
{
  return ntheta_;
}

void
GaussLegendreAngularIntegrator::set_ntheta(int i)
{
  ntheta_ = i;
}

int
GaussLegendreAngularIntegrator::get_nphi(void)
{
  return nphi_;
}

void
GaussLegendreAngularIntegrator::set_nphi(int i)
{
  nphi_ = i;
}

int
GaussLegendreAngularIntegrator::get_Ktheta(void)
{
  return Ktheta_;
}

void
GaussLegendreAngularIntegrator::set_Ktheta(int i)
{
  Ktheta_ = i;
}

int
GaussLegendreAngularIntegrator::get_ntheta_r(void)
{
  return ntheta_r_;
}

void
GaussLegendreAngularIntegrator::set_ntheta_r(int i)
{
  ntheta_r_ = i;
}

int
GaussLegendreAngularIntegrator::get_nphi_r(void)
{
  return nphi_r_;
}

void
GaussLegendreAngularIntegrator::set_nphi_r(int i)
{
  nphi_r_ = i;
}

int
GaussLegendreAngularIntegrator::get_Ktheta_r(void)
{
  return Ktheta_r_;
}

void
GaussLegendreAngularIntegrator::set_Ktheta_r(int i)
{
  Ktheta_r_ = i;
}

double
GaussLegendreAngularIntegrator::sin_theta(SCVector3 &point)
{
  return sin(point.theta());
}

int
GaussLegendreAngularIntegrator::num_angular_points(double r_value, int ir)
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
GaussLegendreAngularIntegrator::angular_point_cartesian(int iangular, SCVector3 &point,
                                                        SCVector3 &integration_point)
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


// Murray93 Integrator used as template
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
          double r = RadInt_->radial_value(ir, nr[icenter], bragg_radius[icenter]);
          point.r() = r;
          dr_dqr2 = RadInt_->get_dr_dqr2();
          radial_multiplier = RadInt_->radial_multiplier(dr_dqr2, nr[icenter]);
          nangular = AngInt_->num_angular_points(r/bragg_radius[icenter],ir);
          AngInt_->angular_weights();
          double radial_int_volume = RadInt_->get_dr_dqr2();
          for (iangular=0; iangular<nangular; iangular++) {
              angular_multiplier =
                AngInt_->angular_point_cartesian(iangular, point, integration_point);
              integration_point += center;
              w=weight_->w(icenter, integration_point, w_gradient);
              //if (w_gradient) weight_->test(icenter, integration_point);
              point_count++;
              double multiplier = angular_multiplier * radial_multiplier;
              double angular_int_volume = AngInt_->sin_theta(point);
              double int_volume = radial_int_volume * angular_int_volume;
              if (do_point(icenter, integration_point, denfunc,
                           w, multiplier,
                           nuclear_gradient, f_gradient, w_gradient)
                  * int_volume < 1e2*DBL_EPSILON
                  && int_volume > 1e2*DBL_EPSILON) {
                  r_done=1;
                  break;
                }
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
  o << node0 << indent << scprintf("nr     = %5d", RadInt_->get_nr()) << endl;
  o << node0 << indent << scprintf("ntheta = %5d", AngInt_->get_ntheta()) << endl;
  o << node0 << indent << scprintf("nphi   = %5d", AngInt_->get_nphi()) << endl;
  o << node0 << indent << scprintf("Ktheta = %5d", AngInt_->get_Ktheta()) << endl;
  o << decindent;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:

