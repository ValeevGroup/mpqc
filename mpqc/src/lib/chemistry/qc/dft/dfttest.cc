//
// dfttest.cc
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

#include <util/misc/formio.h>
#include <util/group/pregtime.h>
#include <chemistry/qc/dft/functional.h>
#include <chemistry/qc/dft/integrator.h>

#include <chemistry/qc/dft/linkage.h>
#include <chemistry/qc/scf/linkage.h>

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

double *
density_matrix(const RefWavefunction &wfn)
{
  int nbasis = wfn->basis()->nbasis();
  RefSymmSCMatrix adens = wfn->alpha_ao_density();
  double * alpha_dmat = new double[(nbasis*(nbasis+1))/2];
  adens->convert(alpha_dmat);
  return alpha_dmat;
}

void
get_density(PointInputData::SpinData &d, const SCVector3 &r,
            const RefWavefunction &wfn, double *pdmat = 0)
{
  double *dmat;
  if (pdmat) dmat = pdmat;
  else dmat = density_matrix(wfn);

  double * bsg_values_ = new double[3*wfn->basis()->nbasis()];
  double * bs_values_ = new double[wfn->basis()->nbasis()];
  wfn->basis()->set_integral(wfn->integral());
  wfn->basis()->grad_values(r,bsg_values_,bs_values_);

  int i, j;

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
  for (i=0; i<3; i++) grad[i] = 0.0;
  for (i=0; i<6; i++) hess[i] = 0.0;

  int nbasis_ = wfn->basis()->nbasis();
  int ij=0;
  for (i=0; i < nbasis_; i++) {
    bvi = bs_values_[i];
    bvix = bsg_values_[i*3+X];
    bviy = bsg_values_[i*3+Y];
    bviz = bsg_values_[i*3+Z];
    for (j=0; j < i; j++,ij++) {
      densij = dmat[ij];
      double bvj = bs_values_[j];
      double bvjx = bsg_values_[j*3+X];
      double bvjy = bsg_values_[j*3+Y];
      double bvjz = bsg_values_[j*3+Z];

      tmp += 2.0*densij*bvi*bvj;

      grad[X] += densij*(bvi*bvjx + bvj*bvix);
      grad[Y] += densij*(bvi*bvjy + bvj*bviy);
      grad[Z] += densij*(bvi*bvjz + bvj*bviz);
      }

    densij = dmat[ij]*bvi;
    tmp += densij*bvi;
    grad[X] += densij*bvix;
    grad[Y] += densij*bviy;
    grad[Z] += densij*bviz;
    ij++;
    }


  d.rho = tmp;
  for (i=0; i<3; i++) d.del_rho[i] = 2.0 * grad[i];
  for (i=0; i<6; i++) d.hes_rho[i] = 2.0 * hess[i];

  d.lap_rho = d.hes_rho[XX] + d.hes_rho[YY] + d.hes_rho[ZZ];

  d.gamma = dot(d.del_rho,d.del_rho);

  delete[] bsg_values_;
  delete[] bs_values_;
  if (!pdmat) delete[] dmat;
}

double
fd_test_do_point(const SCVector3 &point,
                 const RefDenFunctional &func, const RefWavefunction &wfn,
                 double *frozen_dmat = 0)
{
  PointInputData id(point);
  get_density(id.a, point, wfn, frozen_dmat);
  id.compute_derived(0);
  PointOutputData od;
  if ( (id.a.rho + id.b.rho) > 1e2*DBL_EPSILON) func->point(id, od);
  else return 0.0;
  return od.energy;
}

void
fd_test_point(int acenter, const SCVector3 &tpoint,
              const RefDenFunctional &functional, const RefWavefunction &wfn)
{
  SCVector3 point(tpoint);
  RefMolecule mol = wfn->molecule();

  double *fd_grad_f = new double[mol->natom()*3];
  memset(fd_grad_f,0, 3*mol->natom() * sizeof(double));

  double *dmat = density_matrix(wfn);

  // frozen_dmat = dmat makes the overlap derivative contribution 0
  double *frozen_dmat = dmat;

  double delta = 0.0001;
  for (int i=0; i<mol->natom(); i++) {
    for (int j=0; j<3; j++) {
      if (acenter == i) point[j] += delta;
      mol->r(i,j) += delta;
      wfn->obsolete();
      double f_plus = fd_test_do_point(point, functional, wfn, frozen_dmat);
      if (acenter == i) point[j] -= 2*delta;
      mol->r(i,j) -= 2*delta;
      wfn->obsolete();
      double f_minus = fd_test_do_point(point, functional, wfn, frozen_dmat);
      if (acenter == i) point[j] += delta;
      mol->r(i,j) += delta;
      wfn->obsolete();
      fd_grad_f[i*3+j] = (f_plus-f_minus)/(2.0*delta);
      }
    }

  double * bsh_values_ = new double[6*wfn->basis()->nbasis()];
  double * bsg_values_ = new double[3*wfn->basis()->nbasis()];
  double * bs_values_ = new double[wfn->basis()->nbasis()];
  wfn->basis()->hessian_values(point,bsh_values_,bsg_values_,bs_values_);

  PointInputData id(point);
  get_density(id.a, point, wfn);
  id.compute_derived(0);

  PointOutputData od;
  functional->set_compute_potential(1);
  if ( (id.a.rho + id.b.rho) > 1e2*DBL_EPSILON) functional->point(id, od);
  else od.zero();

  double *an_grad_f = new double[mol->natom()*3];
  memset(an_grad_f,0, 3*mol->natom() * sizeof(double));

  functional->gradient(id, od, an_grad_f, acenter, wfn->basis(),
                       dmat, dmat, bs_values_, bsg_values_, bsh_values_);

  cout << " acenter = " << acenter << " point = " << point << endl;
  cout << "FD df/dx:" << endl;
  for (int i=0; i<mol->natom(); i++) {
    cout << scprintf(" % 16.12f % 16.12f % 16.12f",
                     fd_grad_f[3*i+0],
                     fd_grad_f[3*i+1],
                     fd_grad_f[3*i+2])
                     << endl;
    }

  cout << "AN df/dx:" << endl;
  for (int i=0; i<mol->natom(); i++) {
    cout << scprintf(" % 16.12f % 16.12f % 16.12f",
                     an_grad_f[3*i+0],
                     an_grad_f[3*i+1],
                     an_grad_f[3*i+2])
                     << endl;
    }

  delete[] bsh_values_;
  delete[] bsg_values_;
  delete[] bs_values_;
  delete[] dmat;
  delete[] fd_grad_f;
  delete[] an_grad_f;
}

void
fd_test(const RefDenFunctional &functional, const RefWavefunction &wfn)
{
  cout << "fd_test with functional:" << endl;
  cout << functional;
  for (int i=0; i<wfn->molecule()->natom(); i++) {
    for (double x = -1.0; x <= 1.0; x++) {
      for (double y = -1.0; y <= 1.0; y++) {
        for (double z = -1.0; z <= 1.0; z++) {
          fd_test_point(i, SCVector3(x,y,z), functional, wfn);
          }
        }
      }
    }
}

void
fd_e_test(const RefWavefunction &wfn)
{
  RefMolecule mol = wfn->molecule();

  cout << "Testing dE/dx with:" << endl;
  cout << incindent;
  wfn->print();
  cout << decindent;

  double *fd_grad_e = new double[mol->natom()*3];
  memset(fd_grad_e,0, 3*mol->natom() * sizeof(double));

  double delta = 0.0001;
  for (int i=0; i<mol->natom(); i++) {
    for (int j=0; j<3; j++) {
      mol->r(i,j) += delta;
      wfn->obsolete();
      double e_plus = wfn->energy();
      mol->r(i,j) -= 2*delta;
      wfn->obsolete();
      double e_minus = wfn->energy();
      mol->r(i,j) += delta;
      wfn->obsolete();
      fd_grad_e[i*3+j] = (e_plus-e_minus)/(2.0*delta);
      }
    }

  double *an_grad_e = new double[mol->natom()*3];
  memset(an_grad_e,0, 3*mol->natom() * sizeof(double));

  RefSCVector grad = wfn->get_cartesian_gradient();
  grad->convert(an_grad_e);

  cout << "FD dE/dx:" << endl;
  for (int i=0; i<mol->natom(); i++) {
    cout << scprintf(" % 16.12f % 16.12f % 16.12f",
                     fd_grad_e[3*i+0],
                     fd_grad_e[3*i+1],
                     fd_grad_e[3*i+2])
                     << endl;
    }

  cout << "AN dE/dx:" << endl;
  for (int i=0; i<mol->natom(); i++) {
    cout << scprintf(" % 16.12f % 16.12f % 16.12f",
                     an_grad_e[3*i+0],
                     an_grad_e[3*i+1],
                     an_grad_e[3*i+2])
                     << endl;
    }

  delete[] fd_grad_e;
  delete[] an_grad_e;
}

int
main(int argc, char**argv)
{
  int i;
  char *input = (argc > 1)? argv[1] : SRCDIR "/dfttest.in";

  // open keyval input
  RefKeyVal keyval(new ParsedKeyVal(input));

  RefWavefunction  dft        = keyval->describedclassvalue("dft");
  if (dft.nonnull()) {
    cout << "=========== FD dE/dx Tests ===========" << endl;
    fd_e_test(dft);
    }

  cout << "=========== FD df/drho Tests ===========" << endl;
  RefDenFunctional funcs[] = {
    new SlaterXFunctional,
    new Becke88XFunctional,
    new PBECFunctional,
    new PW91CFunctional,
    new P86CFunctional,
    new XalphaFunctional,
    new LYPCFunctional,
    new PW86XFunctional,
    new PW91XFunctional,
    new PBEXFunctional,
    new G96XFunctional,
    new VWN1LCFunctional,
    new VWN2LCFunctional,
    new VWN3LCFunctional,
    new VWN4LCFunctional,
    new VWN5LCFunctional,
    new PW92LCFunctional,
    new PZ81LCFunctional,
    0
  };
  for (i=0; funcs[i]; i++) {
    cout << "-----------------"
         << funcs[i]->class_name()
         << "-----------------"
         << endl;
    funcs[i]->test();
    }

  RefMolecule mol = keyval->describedclassvalue("molecule");
  if (mol.nonnull()) {
    cout << "=========== FD Weight Tests ===========" << endl;
    RefIntegrationWeight weights[] = {
      new BeckeIntegrationWeight,
      0
    };
    for (i=0; weights[i]; i++) {
      cout << "-----------------"
           << weights[i]->class_name()
           << "-----------------"
           << endl;
      weights[i]->init(mol,1.0e-8);
      weights[i]->test();
      }
    }

  RefDenFunctional functional = keyval->describedclassvalue("functional");
  RefWavefunction  wfn        = keyval->describedclassvalue("wfn");
  if (functional.nonnull() && wfn.nonnull()) {
    cout << "=========== FD df/dx Tests ===========" << endl;
    fd_test(functional, wfn);
    }

  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
