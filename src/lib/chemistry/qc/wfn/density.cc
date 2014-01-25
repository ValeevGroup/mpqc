//
// density.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
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

#include <stdexcept>

#include <util/misc/formio.h>
#include <util/render/polygons.h>
#include <math/scmat/local.h>
#include <math/scmat/vector3.h>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/wfn/density.h>
#include <util/misc/scexception.h>

using namespace std;
using namespace sc;

/////////////////////////////////////////////////////////////////////////////
// ElectronDensity

static ClassDesc ElectronDensity_cd(
  typeid(ElectronDensity),"ElectronDensity",1,"public Volume",
  0, create<ElectronDensity>, 0);

ElectronDensity::ElectronDensity(const Ref<KeyVal> &keyval):
  Volume(keyval)
{
  wfn_ << keyval->describedclassvalue("wfn");
}

ElectronDensity::ElectronDensity(const Ref<Wavefunction>& wfn):
  Volume(),
  wfn_(wfn)
{
}

ElectronDensity::~ElectronDensity()
{
}

void
ElectronDensity::compute()
{
  SCVector3 r;
  get_x(r);
  // do_gradient will automatically cause the value to be computed
  if (gradient_needed()) {
      double v[3];
      set_value(wfn_->density_gradient(r,v));
      set_actual_value_accuracy(desired_value_accuracy());
      SCVector3 d(v);
      set_gradient(d);
      set_actual_gradient_accuracy(desired_gradient_accuracy());
    }
  else if (value_needed()) {
      set_value(wfn_->density(r));
      set_actual_value_accuracy(desired_value_accuracy());
    }
  if (hessian_needed()) {
      ExEnv::err0() << indent
           << "ElectronDensity::compute(): hessian isn't yet implemented\n";
      abort();
    }
}

// make a wild guess about the bounding box
void
ElectronDensity::boundingbox(double valuemin,
                             double valuemax,
                             SCVector3& p1, SCVector3& p2)
{
  Molecule& mol = *wfn_->molecule();

  if (mol.natom() == 0) {
      for (int i=0; i<3; i++) p1[i] = p2[i] = 0.0;
    }

  int i;
  for (i=0; i<3; i++) p1[i] = p2[i] = mol.r(0,i);
  for (i=1; i<mol.natom(); i++) {
      for (int j=0; j<3; j++) {
          if (mol.r(i,j) < p1[j]) p1[j] = mol.r(i,j);
          if (mol.r(i,j) > p2[j]) p2[j] = mol.r(i,j);
        }
    }
  for (i=0; i<3; i++) {
      p1[i] = p1[i] - 3.0;
      p2[i] = p2[i] + 3.0;
    }
}

/////////////////////////////////////////////////////////////////////////////
// BatchElectronDensity

static ClassDesc BatchElectronDensity_cd(
  typeid(BatchElectronDensity),"BatchElectronDensity",1,"public Volume",
  0, create<BatchElectronDensity>, 0);

BatchElectronDensity::BatchElectronDensity(const Ref<GaussianBasisSet> &basis,
                                           const Ref<Integral> &integral,
                                           double accuracy):
  Volume()
{
  basis_ = basis;
  integral_ = integral;
  accuracy_ = accuracy;
  zero_pointers();
  using_shared_data_ = false;
  linear_scaling_ = true;
  use_dmat_bound_ = true;
  need_basis_gradient_ = false;
  need_basis_hessian_ = false;
  alpha_dmat_ = 0;
  beta_dmat_ = 0;
  initialized_ = false;
}

BatchElectronDensity::BatchElectronDensity(const Ref<Wavefunction> &wfn,
                                           double accuracy):
  Volume()
{
  wfn_ = wfn;

  basis_ = wfn->basis();
  integral_ = wfn->integral();

  accuracy_ = accuracy;
  zero_pointers();
  using_shared_data_ = false;
  linear_scaling_ = true;
  use_dmat_bound_ = true;
  need_basis_gradient_ = false;
  need_basis_hessian_ = false;
  alpha_dmat_ = 0;
  beta_dmat_ = 0;
  initialized_ = false;
}

BatchElectronDensity::BatchElectronDensity(const Ref<KeyVal> &keyval):
  Volume(keyval)
{
  wfn_ << keyval->describedclassvalue("wfn");
  basis_ = wfn_->basis();
  integral_ = wfn_->integral();
  accuracy_ = keyval->doublevalue("accuracy");
  zero_pointers();
  using_shared_data_ = false;
  linear_scaling_ = true;
  use_dmat_bound_ = true;
  need_basis_gradient_ = false;
  need_basis_hessian_ = false;
  alpha_dmat_ = 0;
  beta_dmat_ = 0;
  initialized_ = false;
}

BatchElectronDensity::BatchElectronDensity(const Ref<BatchElectronDensity>&d,
                                           bool reference_parent_data):
  Volume()
{
  zero_pointers();
  using_shared_data_ = reference_parent_data;
  accuracy_ = d->accuracy_;
  need_basis_gradient_ = d->need_basis_gradient_;
  need_basis_hessian_ = d->need_basis_hessian_;
  initialized_ = d->initialized_;
  spin_polarized_ = d->spin_polarized_;
  if (using_shared_data_) {
      if (d->alpha_dmat_ == 0) {
          throw std::runtime_error("BatchElectronDensity: attempted to use shared data, but parent data not initialized");
        }

      nshell_ = d->nshell_;
      nbasis_ = d->nbasis_;
      basis_ = d->basis_;
      integral_ = d->integral_;
      extent_ = d->extent_;
      alpha_dmat_ = d->alpha_dmat_;
      beta_dmat_ = d->beta_dmat_;
      dmat_bound_ = d->dmat_bound_;
      linear_scaling_ = d->linear_scaling_;
      use_dmat_bound_ = d->use_dmat_bound_;

      init_scratch_data();
    }
}

BatchElectronDensity::~BatchElectronDensity()
{
  clear();
}

void
BatchElectronDensity::zero_pointers()
{
  valdat_ = 0;
  extent_ = 0;

  alpha_dmat_ = 0;
  beta_dmat_ = 0;
  dmat_bound_ = 0;
  contrib_ = 0;
  contrib_bf_ = 0;
  bs_values_ = 0;
  bsg_values_ = 0;
  bsh_values_ = 0;
}

void
BatchElectronDensity::clear()
{
  if (!using_shared_data_) {
      delete extent_;
      delete[] alpha_dmat_;
      delete[] beta_dmat_;
      delete[] dmat_bound_;
    }

  delete[] contrib_;
  delete[] contrib_bf_;
  delete[] bs_values_;
  delete[] bsg_values_;
  delete[] bsh_values_;
  delete valdat_;
  initialized_ = false;

  zero_pointers();
}

void
BatchElectronDensity::init()
{
  if (using_shared_data_)
      throw std::runtime_error("BatchElectronDensity::init: should not be called if using_shared_data");

  clear();
  init_common_data();
  init_scratch_data();
  initialized_ = true;
}

void
BatchElectronDensity::init_common_data()
{
  nshell_ = basis_->nshell();
  nbasis_ = basis_->nbasis();

  if (linear_scaling_) {
      extent_ = new ShellExtent;
      extent_->init(basis_);
    }

  dmat_bound_ = new double[(nshell_*(nshell_+1))/2];
}

void
BatchElectronDensity::set_densities(const Ref<Wavefunction> &wfn)
{
  set_densities(wfn->alpha_density(), wfn->beta_density());
}

void
BatchElectronDensity::set_densities(const RefSymmSCMatrix &aden,
                                    const RefSymmSCMatrix &bden)
{
  RefSymmSCMatrix ad = aden;
  RefSymmSCMatrix bd = bden;

  if (aden == bden || bden == 0) spin_polarized_ = 0;
  else spin_polarized_ = 1;
  
  alpha_dmat_ = new double[(nbasis_*(nbasis_+1))/2];
  beta_dmat_ = 0;
  if (spin_polarized_) {
      beta_dmat_ = new double[(nbasis_*(nbasis_+1))/2];
    }

  ad->convert(alpha_dmat_);
  if (spin_polarized_) bd->convert(beta_dmat_);

  int ij = 0;
  for (int i=0; i<nshell_; i++) {
      int ni = basis_->shell(i).nfunction();
      for (int j=0; j<=i; j++,ij++) {
          int nj = basis_->shell(j).nfunction();
          double bound = 0.0;
          int ibf = basis_->shell_to_function(i);
          for (int k=0; k<ni; k++,ibf++) {
              int lmax = nj-1;
              if (i==j) lmax = k;
              int jbf = basis_->shell_to_function(j);
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
BatchElectronDensity::init_scratch_data()
{
  contrib_ = new int[nshell_];
  contrib_bf_ = new int[nbasis_];
  bs_values_ = new double[nbasis_];
  bsg_values_ = new double[3*nbasis_];
  bsh_values_ = new double[6*nbasis_];
  valdat_ = new GaussianBasisSet::ValueData(basis_, integral_);
}

void
BatchElectronDensity::compute_basis_values(const SCVector3&r)
{
  // only consider those shells for which phi_i * (Max_j D_ij phi_j) > tol
  if (linear_scaling_ && use_dmat_bound_ && extent_ != 0) {
      const std::vector<ExtentData> &cs = extent_->contributing_shells(r[0],r[1],r[2]);
      ncontrib_ = 0;
      for (int i=0; i<cs.size(); i++) {
          int ish = cs[i].shell;
          int contrib = 0;
          for (int j=0; j<cs.size(); j++) {
              int jsh = cs[j].shell;
              int ijsh = (ish>jsh)?((ish*(ish+1))/2+jsh):((jsh*(jsh+1))/2+ish);
//               std::cout << "cs[i].bound = " << cs[i].bound << std::endl;
//               std::cout << "cs[j].bound = " << cs[j].bound << std::endl;
//               std::cout << "dmat_bound_[ijsh] = " << dmat_bound_[ijsh] << std::endl;
//               std::cout << "accuracy_ = " << accuracy_ << std::endl;
              if (cs[i].bound*cs[j].bound*dmat_bound_[ijsh] > 0.00001*accuracy_) {
                  contrib = 1;
                  break;
                }
            }
          if (contrib) {
              contrib_[ncontrib_++] = ish;
            }
        }
    }
  else if (linear_scaling_ && extent_ != 0) {
      const std::vector<ExtentData> &cs = extent_->contributing_shells(r[0],r[1],r[2]);
      ncontrib_ = cs.size();
      for (int i=0; i<ncontrib_; i++) {
          contrib_[i] = cs[i].shell;
        }
    }
  else {
      ncontrib_ = nshell_;
      for (int i=0; i<nshell_; i++) contrib_[i] = i;
    }

  ncontrib_bf_ = 0;
  for (int i=0; i<ncontrib_; i++) {
      int nbf = basis_->shell(contrib_[i]).nfunction();
      int bf = basis_->shell_to_function(contrib_[i]);
      for (int j=0; j<nbf; j++, bf++) {
          contrib_bf_[ncontrib_bf_++] = bf;
        }
    }

  // compute the basis set values
  double *bsv = bs_values_;
  double *bsg = ((need_basis_gradient_||need_gradient_)?bsg_values_:0);
  double *bsh = ((need_basis_hessian_||need_hessian_)?bsh_values_:0);
  for (int i=0; i<ncontrib_; i++) {
      basis_->hessian_shell_values(r,contrib_[i],valdat_,bsh,bsg,bsv);
      int shsize = basis_->shell(contrib_[i]).nfunction();

      if (bsh) bsh += 6 * shsize;
      if (bsg) bsg += 3 * shsize;
      if (bsv) bsv += shsize;
    }
}

void
BatchElectronDensity::compute_spin_density(const double *RESTRICT dmat,
                                           double *RESTRICT rho,
                                           double *RESTRICT pgrad,
                                           double *RESTRICT phess)
{
  int i, j;

  double tmp = 0.0;
  double densij;
  double bvi, bvix, bviy, bviz;
  double bvixx, bviyx, bviyy, bvizx, bvizy, bvizz;

  double grad[3];
  double hess[6];

  double *RESTRICT bs_vals = bs_values_;
  double *RESTRICT bsg_vals = bsg_values_;
  double *RESTRICT bsh_vals = bsh_values_;

  if (need_gradient_) for (i=0; i<3; i++) grad[i] = 0.0;
  if (need_hessian_) for (i=0; i<6; i++) hess[i] = 0.0;

  if (need_gradient_ && !need_hessian_) {
      for (i=0; i < ncontrib_bf_; i++) {
          int it = contrib_bf_[i];
          bvi = bs_vals[i];
          bvix = bsg_vals[i*3+X];
          bviy = bsg_vals[i*3+Y];
          bviz = bsg_vals[i*3+Z];
          int j3 = 0, j6 = 0;
          int itoff = (it*(it+1))>>1;
          int itjt;
          double t = 0.0;
          for (j=0; j < i; j++) {
              int jt = contrib_bf_[j];
              itjt = itoff+jt;

              densij = dmat[itjt];
              double bvj = bs_vals[j];

              t += densij*bvi*bvj;

              double bvjx, bvjy, bvjz;
              bvjx = bsg_vals[j3+X];
              bvjy = bsg_vals[j3+Y];
              bvjz = bsg_vals[j3+Z];
              grad[X] += densij*(bvi*bvjx + bvj*bvix);
              grad[Y] += densij*(bvi*bvjy + bvj*bviy);
              grad[Z] += densij*(bvi*bvjz + bvj*bviz);
              j3 += 3;
            }
          densij = dmat[itoff+it]*bvi;
          tmp += t + 0.5*densij*bvi;
          grad[X] += densij*bvix;
          grad[Y] += densij*bviy;
          grad[Z] += densij*bviz;
        }
    }
  else if (need_gradient_ || need_hessian_) {
      for (i=0; i < ncontrib_bf_; i++) {
          int it = contrib_bf_[i];
          bvi = bs_vals[i];
          if (need_gradient_) {
              bvix = bsg_vals[i*3+X];
              bviy = bsg_vals[i*3+Y];
              bviz = bsg_vals[i*3+Z];
            }
          if (need_hessian_) {
              bvixx = bsh_vals[i*6+XX];
              bviyx = bsh_vals[i*6+YX];
              bviyy = bsh_vals[i*6+YY];
              bvizx = bsh_vals[i*6+ZX];
              bvizy = bsh_vals[i*6+ZY];
              bvizz = bsh_vals[i*6+ZZ];
            }
          int j3 = 0, j6 = 0;
          int itoff = (it*(it+1))>>1;
          int itjt;
          double t = 0.0;
          for (j=0; j < i; j++) {
              int jt = contrib_bf_[j];
              itjt = itoff+jt;

              densij = dmat[itjt];
              double bvj = bs_vals[j];

              t += densij*bvi*bvj;

              double bvjx, bvjy, bvjz;
              if (need_gradient_) {
                  bvjx = bsg_vals[j3+X];
                  bvjy = bsg_vals[j3+Y];
                  bvjz = bsg_vals[j3+Z];
                  grad[X] += densij*(bvi*bvjx + bvj*bvix);
                  grad[Y] += densij*(bvi*bvjy + bvj*bviy);
                  grad[Z] += densij*(bvi*bvjz + bvj*bviz);
                  j3 += 3;
                }

              if (need_hessian_) {
                  double bvjxx = bsh_vals[j6+XX];
                  double bvjyx = bsh_vals[j6+YX];
                  double bvjyy = bsh_vals[j6+YY];
                  double bvjzx = bsh_vals[j6+ZX];
                  double bvjzy = bsh_vals[j6+ZY];
                  double bvjzz = bsh_vals[j6+ZZ];

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
          bvi = bs_vals[i];
          int itoff = (it*(it+1))>>1;
          int itjt;
          double t = 0.0;
          for (j=0; j < i; j++) {
              int jt = contrib_bf_[j];
              itjt = itoff+jt;

              densij = dmat[itjt];
              double bvj = bs_vals[j];

              t += densij*bvi*bvj;
            }
          densij = dmat[itoff+it]*bvi;
          tmp += t + 0.5*densij*bvi;
        }
    }
  if (rho!=0) *rho = tmp;
  if (need_gradient_) for (i=0; i<3; i++) pgrad[i] = grad[i];
  if (need_hessian_) for (i=0; i<6; i++) phess[i] = hess[i];
}

void
BatchElectronDensity::compute_density(const SCVector3 &r,
                                      double *adens,
                                      double *agrad,
                                      double *ahess,
                                      double *bdens,
                                      double *bgrad,
                                      double *bhess)
{
  if (alpha_dmat_ == 0) {
      if (wfn_ == 0) {
          throw ProgrammingError("BatchElectronDensity::compute_density: "
                                 "set_densities must be used to initialize "
                                 "object if wfn is not given",
                                 __FILE__, __LINE__);
        }
      else {
          if (!initialized_) {
              init();
            }
          set_densities(wfn_);
        }
    }

  need_gradient_ = (agrad!=0) || (bgrad!=0);
  need_hessian_ = (ahess!=0) || (bhess!=0);

  compute_basis_values(r);

  compute_spin_density(alpha_dmat_,
                       adens,
                       agrad,
                       ahess);

  bool mismatch = (adens==0 && bdens!=0)
                  ||(agrad==0 && bgrad!=0)
                  ||(ahess==0 && bhess!=0);

  if (spin_polarized_ || mismatch) {
      compute_spin_density(beta_dmat_,
                           bdens,
                           bgrad,
                           bhess);
    }
  else {
      if (bdens!=0) *bdens = *adens;
      if (bgrad!=0)
          for (int i=0;i<3;i++) bgrad[i] = agrad[i];
      if (bhess!=0)
          for (int i=0;i<6;i++) bhess[i] = ahess[i];
    }

  if (adens!=0) *adens *= 2.0;
  if (agrad!=0)
      for (int i=0;i<3;i++) agrad[i] *= 2.0;
  if (ahess!=0)
      for (int i=0;i<6;i++) ahess[i] *= 2.0;
  if (bdens!=0) *bdens *= 2.0;
  if (bgrad!=0)
      for (int i=0;i<3;i++) bgrad[i] *= 2.0;
  if (bhess!=0)
      for (int i=0;i<6;i++) bhess[i] *= 2.0;

//   if (agrad) {
//       cout << scprintf("compute_density: agrad = %12.8f %12.8f %12.8f",
//                        agrad[0], agrad[1], agrad[2])
//            << endl;
//     }

//   cout << "compute_density: exiting"
//        << std::endl;

}

void
BatchElectronDensity::compute()
{
  SCVector3 r;
  get_x(r);

  double val;
  double grad[3];
  double hess[6];
  double aval;
  double agrad[3];
  double ahess[6];
  double bval;
  double bgrad[3];
  double bhess[6];
  compute_density(r,
                  &aval,
                  (gradient_needed()?agrad:0),
                  (hessian_needed()?ahess:0),
                  &bval,
                  (gradient_needed()?bgrad:0),
                  (hessian_needed()?bhess:0));
  val = aval + bval;
  for (int i=0; i<3; i++) grad[i] = agrad[i] + bgrad[i];
  for (int i=0; i<6; i++) hess[i] = ahess[i] + bhess[i];

  if (value_needed()) {
      set_value(val);
      set_actual_value_accuracy(desired_value_accuracy());
    }
  if (gradient_needed()) {
      set_value(val);
      set_actual_value_accuracy(desired_value_accuracy());
      SCVector3 d(grad);
      set_gradient(d);
      set_actual_gradient_accuracy(desired_gradient_accuracy());
    }
  if (hessian_needed()) {
      ExEnv::err0() << indent
           << "BatchElectronDensity::compute(): hessian isn't yet implemented\n";
      abort();
    }
}

void
BatchElectronDensity::boundingbox(double valuemin,
                             double valuemax,
                             SCVector3& p1, SCVector3& p2)
{
#if 0
  // this is a very conservative bounding box
  // also, this code is not correct since extent is not
  // necessarily initialized
  if (alpha_dmat_ == 0) init();
  for (int i=0; i<3; i++) p1[i] = extent_->lower(i);
  for (int i=0; i<3; i++) p2[i] = extent_->upper(i);
#else
  Molecule& mol = *basis_->molecule();

  if (mol.natom() == 0) {
      for (int i=0; i<3; i++) p1[i] = p2[i] = 0.0;
    }

  int i;
  for (i=0; i<3; i++) p1[i] = p2[i] = mol.r(0,i);
  for (i=1; i<mol.natom(); i++) {
      for (int j=0; j<3; j++) {
          if (mol.r(i,j) < p1[j]) p1[j] = mol.r(i,j);
          if (mol.r(i,j) > p2[j]) p2[j] = mol.r(i,j);
        }
    }
  for (i=0; i<3; i++) {
      p1[i] = p1[i] - 3.0;
      p2[i] = p2[i] + 3.0;
    }
#endif
}

/////////////////////////////////////////////////////////////////////////////
// WriteElectronDensity

static ClassDesc WriteElectronDensity_cd(
    typeid(WriteElectronDensity),"WriteElectronDensity",1,
    "public WriteGrid", 0, create<WriteElectronDensity>, 0);

WriteElectronDensity::WriteElectronDensity(const Ref<KeyVal> &keyval):
  WriteGrid(keyval)
{
  wfn_ << keyval->describedclassvalue("wfn");
  if (wfn_ == 0) {
      InputError ex("valid \"wfn\" missing",
                    __FILE__, __LINE__, "wfn", "(null)", class_desc());
      try {
          ex.elaborate()
              << "WriteElectronDensity KeyVal ctor requires"
              << " that \"wfn\" specifies an object"
              << " of type Wavefunction" << std::endl;
        }
      catch (...) {}
      throw ex;
    }
  
  if (keyval->exists("type")) {
      type_ = keyval->stringvalue("type");
      if (type_ == "alpha") {
          density_function_ = &WriteElectronDensity::df_alpha;
        }
      else if (type_ == "beta") {
          density_function_ = &WriteElectronDensity::df_beta;
        }
      else if (type_ == "sum") {
          density_function_ = &WriteElectronDensity::df_sum;
        }
      else if (type_ == "spin") {
          density_function_ = &WriteElectronDensity::df_spin;
        }
      else {
          InputError ex("valid \"type\" missing",
                        __FILE__, __LINE__, "type", "(null)", class_desc());
          try {
              ex.elaborate()
                  << "WriteElectronDensity KeyVal ctor requires"
                  << " that \"type\" is one of \"alpha\", \"beta\","
                  << " \"sum\" or \"spin\". It has been set to \""
                  << type_ << "\"." << std::endl;
            }
          catch (...) {}
          throw ex;
        }
    }
  else {
      type_ = "sum";
      density_function_ = &WriteElectronDensity::df_sum;
    }

  KeyValValuedouble default_accuracy(DBL_EPSILON);
  accuracy_ = keyval->doublevalue("accuracy", default_accuracy);
}

void
WriteElectronDensity::initialize()
{
  bed_ = new BatchElectronDensity(wfn_,accuracy_);
}

void
WriteElectronDensity::label(char* buffer)
{
  sprintf(buffer, "WriteElectronDensity_%s", type_.c_str());
}

Ref<Molecule>
WriteElectronDensity::get_molecule()
{
  return wfn_->molecule();
}

double
WriteElectronDensity::calculate_value(SCVector3 point)
{
  double alpha_density, beta_density;
  bed_->compute_density(point,
                        &alpha_density, 0, 0,
                        &beta_density, 0, 0);
  return (*this.*density_function_)(alpha_density, beta_density);
}

double
WriteElectronDensity::df_alpha(double alpha, double beta) {
  return alpha;
}

double
WriteElectronDensity::df_beta(double alpha, double beta) {
  return beta;
}

double
WriteElectronDensity::df_sum(double alpha, double beta) {
  return alpha + beta;
}

double
WriteElectronDensity::df_spin(double alpha, double beta) {
  return alpha - beta;
}


/////////////////////////////////////////////////////////////////////////////
// DensityColorizer

static ClassDesc DensityColorizer_cd(
  typeid(DensityColorizer),"DensityColorizer",1,"public MoleculeColorizer",
  0, create<DensityColorizer>, 0);

DensityColorizer::DensityColorizer(const Ref<KeyVal>&keyval):
  MoleculeColorizer(keyval)
{
  wfn_ << keyval->describedclassvalue("wfn");
  reference_ = keyval->doublevalue("reference");
  if (keyval->error() == KeyVal::OK) have_reference_ = 1;
  else have_reference_ = 0;
  scale_ = keyval->doublevalue("scale");
  if (keyval->error() == KeyVal::OK) have_scale_ = 1;
  else have_scale_ = 0;
}

DensityColorizer::~DensityColorizer()
{
}

void
DensityColorizer::colorize(const Ref<RenderedPolygons> &poly)
{
  const double base = 0.3;
  int i;
  int nvertex = poly->nvertex();

  if (nvertex) {
      double *data = new double[nvertex];

      for (i=0; i<nvertex; i++) {
          SCVector3 v(poly->vertex(i));
          data[i] = wfn_->density(v);
        }

      double min = data[0], max = data[0];
      for (i=1; i<nvertex; i++) {
          if (min > data[i]) min = data[i];
          if (max < data[i]) max = data[i];
        }

      double center, scale;

      if (have_reference_) center = reference_;
      else center = (max+min)/2.0; 

      double maxdiff = fabs(max - center);
      double mindiff = fabs(min - center);

      if (have_scale_) {
          scale = scale_;
        }
      else {
          if (maxdiff>mindiff && maxdiff>1.0e-6) scale = (1.0-base)/maxdiff;
          else if (mindiff>1.0e-6) scale = (1.0-base)/mindiff;
          else scale = (1.0-base);
        }

      ExEnv::out0() << indent << "DensityColorizer:"
           << scprintf(" reference=%6.5f", center)
           << scprintf(" scale=%8.4f",scale)
           << scprintf(" (%6.5f<=rho<=%6.5f)", max, min)
           << endl;
      for (i=0; i<nvertex; i++) {
          data[i] = (data[i]-center)*scale;
        }

      for (i=0; i<nvertex; i++) {
          Color c;
          if (data[i] < 0.0) c.set_rgb(-data[i]+base,0.3,0.3);
          else c.set_rgb(0.3,0.3,data[i]+base);
          poly->set_vertex_color(i,c);
        }

      delete[] data;
    }
}

/////////////////////////////////////////////////////////////////////////////
// GradDensityColorizer

static ClassDesc GradDensityColorizer_cd(
  typeid(GradDensityColorizer),"GradDensityColorizer",1,"public MoleculeColorizer",
  0, create<GradDensityColorizer>, 0);

GradDensityColorizer::GradDensityColorizer(const Ref<KeyVal>&keyval):
  MoleculeColorizer(keyval)
{
  wfn_ << keyval->describedclassvalue("wfn");
  reference_ = keyval->doublevalue("reference");
  if (keyval->error() == KeyVal::OK) have_reference_ = 1;
  else have_reference_ = 0;
  scale_ = keyval->doublevalue("scale");
  if (keyval->error() == KeyVal::OK) have_scale_ = 1;
  else have_scale_ = 0;
}

GradDensityColorizer::~GradDensityColorizer()
{
}

void
GradDensityColorizer::colorize(const Ref<RenderedPolygons> &poly)
{
  const double base = 0.3;
  int i;
  int nvertex = poly->nvertex();

  Ref<BatchElectronDensity> den = new BatchElectronDensity(wfn_);

  if (nvertex) {
      double *data = new double[nvertex];

      for (i=0; i<nvertex; i++) {
          SCVector3 v(poly->vertex(i));
          SCVector3 g;
          den->set_x(v);
          den->get_gradient(g);
          data[i] = g.norm();
        }

      double min = data[0], max = data[0];
      for (i=1; i<nvertex; i++) {
          if (min > data[i]) min = data[i];
          if (max < data[i]) max = data[i];
        }

      double center, scale;

      if (have_reference_) center = reference_;
      else center = (max+min)/2.0; 

      double maxdiff = fabs(max - center);
      double mindiff = fabs(min - center);

      if (have_scale_) {
          scale = scale_;
        }
      else {
          if (maxdiff>mindiff && maxdiff>1.0e-6) scale = (1.0-base)/maxdiff;
          else if (mindiff>1.0e-6) scale = (1.0-base)/mindiff;
          else scale = (1.0-base);
        }

      ExEnv::out0() << indent << "GradDensityColorizer:"
           << scprintf(" reference=%6.5f", center)
           << scprintf(" scale=%6.2f",scale)
           << scprintf(" (%6.5f<=rho<=%6.5f)", max, min)
           << endl;
      for (i=0; i<nvertex; i++) {
          data[i] = (data[i]-center)*scale;
        }

      for (i=0; i<nvertex; i++) {
          Color c;
          if (data[i] > 0.0) c.set_rgb(data[i]+base,0.3,0.3);
          else c.set_rgb(0.3,0.3,-data[i]+base);
          poly->set_vertex_color(i,c);
        }

      delete[] data;
    }
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
