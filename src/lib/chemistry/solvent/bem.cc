//
// bem.cc
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

#include <stdio.h>
#include <util/misc/math.h>
#include <util/misc/formio.h>
#include <util/misc/regtime.h>
#include <math/scmat/matrix.h>
#include <math/scmat/vector3.h>
#include <math/scmat/local.h>
#include <chemistry/solvent/bem.h>

using namespace std;
using namespace sc;

static ClassDesc BEMSolvent_cd(
  typeid(BEMSolvent),"BEMSolvent",1,"public DescribedClass",
  0, create<BEMSolvent>, 0);

BEMSolvent::BEMSolvent(const Ref<KeyVal>& keyval)
{
  vertex_area_ = 0;

  matrixkit_ = new LocalSCMatrixKit;
  
  debug_ = keyval->intvalue("debug");
  
  solute_ << keyval->describedclassvalue("solute");

  solvent_ << keyval->describedclassvalue("solvent");
  // Use the aug-cc-pVQZ MP2 optimum geometry for H2O as default
  if (solvent_ == 0) {
      solvent_ = new Molecule;
      solvent_->add_atom(8, 0.0000000000,  0.0000000000, -0.1265941233);
      solvent_->add_atom(1, 0.0000000000,  1.4304840085,  0.9856159541);
      solvent_->add_atom(1, 0.0000000000, -1.4304840085,  0.9856159541);
    }

  solvent_density_ = keyval->doublevalue("solvent_density");
  // use as default the number density of water in au^-3, T=25 C, P=101325 Pa
  if (keyval->error() != KeyVal::OK) solvent_density_ = 0.004938887;

  surf_ << keyval->describedclassvalue("surface");

  dielectric_constant_ = keyval->doublevalue("dielectric_constant");
  if (keyval->error() != KeyVal::OK) dielectric_constant_ = 78.0;

  grp_ = MessageGrp::get_default_messagegrp();
}

BEMSolvent::~BEMSolvent()
{
}

double**
BEMSolvent::alloc_array(int n, int m)
{
  double ** result = new double*[n];
  result[0] = new double[n*m];
  for (int i=1; i<n; i++) {
      result[i] = &result[i-1][m];
    }
  return result;
}

void
BEMSolvent::free_array(double** array)
{
  if (!array) return;
  delete[] array[0];
  delete[] array;
}

void
BEMSolvent::charge_positions(double**pos)
{
  int i,j;
  int n = ncharge();
  for (i=0; i<n; i++) {
      const SCVector3& p = surf_->vertex(i)->point();
      for (j=0; j<3; j++) { 
          pos[i][j] = p[j];
        }
    }
}

void
BEMSolvent::normals(double**norms)
{
  int i,j;
  int n = ncharge();
  for (i=0; i<n; i++) {
      const SCVector3& p = surf_->vertex(i)->normal();
      for (j=0; j<3; j++) { 
          norms[i][j] = p[j];
        }
    }
}

void
BEMSolvent::init()
{
  surf_->clear();
  surf_->init();
  system_matrix_i_ = 0;

  f_ = (1.0-dielectric_constant_)/(2.0*M_PI*(1.0+dielectric_constant_));

  if (vertex_area_) delete[] vertex_area_;
  vertex_area_ = new double[ncharge()];
  for (int i=0; i<ncharge(); i++) vertex_area_[i] = 0.0;
  TriangulatedSurfaceIntegrator triint(surf_.pointer());
  for (triint = 0; triint.update(); triint++) {
      int j0 = triint.vertex_number(0);
      int j1 = triint.vertex_number(1);
      int j2 = triint.vertex_number(2);
      double r = triint.r();
      double s = triint.s();
      double dA = triint.w();
      vertex_area_[j0] += dA * (1 - r - s);
      vertex_area_[j1] += dA * r;
      vertex_area_[j2] += dA * s;
    }
}

void
BEMSolvent::done(int clear_surface)
{
  if (clear_surface) surf_->clear();
  system_matrix_i_ = 0;

  if (vertex_area_) delete[] vertex_area_;
  vertex_area_ = 0;
}

void
BEMSolvent::charges_to_surface_charge_density(double *charges)
{
  for (int i=0; i<ncharge(); i++) charges[i] /= vertex_area_[i];
}

void
BEMSolvent::surface_charge_density_to_charges(double *charges)
{
  for (int i=0; i<ncharge(); i++) charges[i] *= vertex_area_[i];
}

double
BEMSolvent::polarization_charge(double *charges)
{
  double charge = 0.0;
  int n = ncharge();
  for (int i=0; i<n; i++) charge += charges[i];
  return charge;
}

// the passed enclosed_charge is determined by the called and
// might different from the enclosed charge computed by Gauss's
// law, which is stored as computed_enclosed_charge_
void
BEMSolvent::normalize_charge(double enclosed_charge, double* charges)
{
  int i;
  double expected_charge = enclosed_charge
                         * (1.0/dielectric_constant_ - 1.0);
  double charge = 0.0;
  double charge_pos = 0.0;
  double charge_neg = 0.0;
  int n = ncharge();
  for (i=0; i<n; i++) {
      charge += charges[i];
      if (charges[i] > 0.0) charge_pos += charges[i];
      else charge_neg += charges[i];
    }

  double scale_pos = 1.0;
  double scale_neg = 1.0;
  if (charge_pos > 1.0e-4 && charge_neg < -1.0e-4) {
      scale_pos += (expected_charge-charge)/(2.0*charge_pos);
      scale_neg += (expected_charge-charge)/(2.0*charge_neg);
    }
  else if (charge_pos > 1.0e-4) {
      scale_pos += (expected_charge-charge)/charge_pos;
    }
  else if (charge_neg < -1.0e-4) {
      scale_neg += (expected_charge-charge)/charge_neg;
    }

  double new_charge = 0.0;
  for (i=0; i<n; i++) {
      if (charges[i] > 0.0) charges[i] *= scale_pos;
      else charges[i] *= scale_neg;
      new_charge += charges[i];
    }

  if (fabs(new_charge - expected_charge) > 1.0e-3) {
      ExEnv::outn() << "BEMSolvent:normalize_charge: failed:" << endl
           << "new_charge = " << new_charge << endl
           << "expected_charge = " << expected_charge << endl;
      abort();
    }

  if (debug_) {
      ExEnv::out0() << indent
           << "BEMSolvent:normalize_charge:"
           << endl << indent
           << scprintf("  integrated surface charge = %20.15f", charge)
           << endl << indent
           << scprintf("  expected surface charge = %20.15f", expected_charge)
           << endl;
    }
}

void
BEMSolvent::init_system_matrix()
{
  int i, j;
  int n = ncharge();

  RefSCDimension d = new SCDimension(n);
  RefSCMatrix system_matrix(d,d,matrixkit());
  system_matrix.assign(0.0);

  Timer tim("precomp");
  // precompute some arrays
  TriangulatedSurfaceIntegrator triint(surf_.pointer());
  int n_integration_points = triint.n();
  SCVector3 *surfpv = new SCVector3[n_integration_points];
  double *rfdA = new double[n_integration_points];
  double *sfdA = new double[n_integration_points];
  double *rsfdA = new double[n_integration_points];
  int *j0 = new int[n_integration_points];
  int *j1 = new int[n_integration_points];
  int *j2 = new int[n_integration_points];
  for (triint=0, i=0; i<n_integration_points&&triint.update(); i++,triint++) {
      surfpv[i] = triint.current()->point();
      j0[i] = triint.vertex_number(0);
      j1[i] = triint.vertex_number(1);
      j2[i] = triint.vertex_number(2);
      double r = triint.r();
      double s = triint.s();
      double rs = 1 - r - s;
      double dA = triint.w();
      double fdA = - f_ * dA;
      rfdA[i] = r * fdA;
      sfdA[i] = s * fdA;
      rsfdA[i] = rs * fdA;
    }
  tim.exit("precomp");

  tim.enter("sysmat");
  double *sysmati = new double[n];
  RefSCVector vsysmati(system_matrix->rowdim(),system_matrix->kit());
  // loop thru all the vertices
  for (i = 0; i<n; i++) {
      memset(sysmati,0,sizeof(double)*n);
      Ref<Vertex> v = surf_->vertex(i);
      const SCVector3& pv = v->point();
      const SCVector3& nv = v->normal();
      // integrate over the surface
      for (j = 0; j < n_integration_points; j++) {
          SCVector3 diff(pv - surfpv[j]);
          double normal_component = diff.dot(nv);
          double diff2 = diff.dot(diff);
          if (diff2 <= 1.0e-8) {
              // The self term must not be included here.  This
              // case shouldn't occur for the usual integrators
              // so abort.
              ExEnv::errn() << "BEMSolvent: integrator gave the self term" << endl;
              abort();
            }
          double denom = diff2*sqrt(diff2);
          double common_factor = normal_component/denom;
          sysmati[j0[j]] += common_factor * rsfdA[j];
          sysmati[j1[j]] += common_factor * rfdA[j];
          sysmati[j2[j]] += common_factor * sfdA[j];
        }
      vsysmati->assign(sysmati);
      system_matrix->assign_row(vsysmati,i);
    }
  tim.exit("sysmat");

  delete[] surfpv;
  delete[] rfdA;
  delete[] sfdA;
  delete[] rsfdA;
  delete[] j0;
  delete[] j1;
  delete[] j2;
  delete[] sysmati;

  tim.enter("AV");
  double A = 0.0;
  double V = 0.0;
  for (triint = 0; triint.update(); triint++) {
      V += triint.weight()*triint.dA()[2]*triint.current()->point()[2];
      A += triint.w();
    }
  area_ = A;
  volume_ = V;
  tim.exit("AV");

  ExEnv::out0() << indent
       << scprintf("Solvent Accessible Surface:") << endl
       << indent
       << scprintf("  Area = %15.10f ", A)
       << scprintf("Volume = %15.10f ", V)
       << scprintf("Nvertex = %3d", n) << endl;

  // Add I to the system matrix.
  system_matrix->shift_diagonal(1.0);

  //system_matrix->print("System Matrix");

  tim.enter("inv");
  system_matrix->invert_this();
  system_matrix_i_ = system_matrix;
  tim.exit("inv");

  //system_matrix_i_->print("System Matrix Inverse");
}

void
BEMSolvent::compute_charges(double* efield_dot_normals, double* charges)
{
  Timer tim;

  if (system_matrix_i_ == 0) {
      tim.enter("sysmat");
      init_system_matrix();
      tim.exit("sysmat");
    }

  tim.enter("qenq");
  double efield_dot_normal = 0.0;
  int n = ncharge();
  for (int i=0; i<n; i++)
      efield_dot_normal += efield_dot_normals[i] * vertex_area_[i];
  tim.exit("qenq");

  computed_enclosed_charge_ = efield_dot_normal/(4.0*M_PI);

  if (debug_) {
      double computed_expected_charge = computed_enclosed_charge_
                                      * (1.0/dielectric_constant_ - 1.0);

      ExEnv::out0() << indent
         << scprintf("BEMSolvent:compute_charges: encl q = %20.15f",
                     computed_enclosed_charge_)
         << endl << indent
         << scprintf("BEMSolvent:compute_charges: exp surface q = %20.15f",
                     computed_expected_charge) << endl;
    }

  tim.enter("scomp");
  RefSCVector edotn(system_matrix_i_.coldim(),matrixkit());
  edotn.assign(efield_dot_normals);
  //edotn.print("E dot normals");
  edotn.scale(f_);
  RefSCVector chrg = system_matrix_i_ * edotn;
  //chrg.print("Charges");
  chrg.convert(charges);
  tim.exit("scomp");

  tim.enter("stoq");
  surface_charge_density_to_charges(charges);
  tim.exit("stoq");
}

double
BEMSolvent::nuclear_charge_interaction_energy(double *nuclear_charge,
                                              double** charge_positions,
                                              double* charge)
{
  double energy = 0.0;
  int natom = solute_->natom();
  for (int i=0; i<natom; i++) {
      for (int j=0; j<ncharge(); j++) {
          double r2 = 0.0;
          for (int k=0; k<3; k++) {
              double r = charge_positions[j][k] - solute_->r(i,k);
              r2 += r*r;
            }
          energy += nuclear_charge[i] * charge[j] / sqrt(r2);
        }
    }
  return energy;
}

double
BEMSolvent::nuclear_interaction_energy(double** charge_positions,
                                       double* charge)
{
  double energy = 0.0;
  int natom = solute_->natom();
  for (int i=0; i<natom; i++) {
      for (int j=0; j<ncharge(); j++) {
          double r2 = 0.0;
          for (int k=0; k<3; k++) {
              double r = charge_positions[j][k] - solute_->r(i,k);
              r2 += r*r;
            }
          energy += double(solute_->Z(i)) * charge[j] / sqrt(r2);
        }
    }
  return energy;
}

double
BEMSolvent::self_interaction_energy(double** charge_positions,
                                    double* charge)
{
  int i,j;

  charges_to_surface_charge_density(charge);

  TriangulatedSurfaceIntegrator triint(surf_.pointer());
  int n_integration_points = triint.n();
  SCVector3 *points = new SCVector3[n_integration_points];
  double *charges = new double[n_integration_points];

  double energy = 0.0;
  for (triint=0, i=0; i<n_integration_points&&triint.update(); i++,triint++) {
      points[i] = triint.current()->point();
      int v0 = triint.vertex_number(0);
      int v1 = triint.vertex_number(1);
      int v2 = triint.vertex_number(2);
      double r = triint.r();
      double s = triint.s();
      double rs = 1.0 - r - s;
      double dA = triint.w();
      charges[i] = (charge[v0]*rs + charge[v1]*r + charge[v2]*s)*dA;
      energy += 0.0; // is this good enough for the self term?
    }
  for (i=0; i<n_integration_points; i++) {
      double chargesi = charges[i];
      SCVector3 pointsi(points[i]);
      for (j = 0; j<i; j++) {
          energy += chargesi*charges[j]/pointsi.dist(points[j]);
        }
    }

  delete[] points;
  delete[] charges;

  surface_charge_density_to_charges(charge);

  return energy;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
