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
#include <math.h>
#include <util/misc/formio.h>
#include <math/scmat/matrix.h>
#include <math/scmat/vector3.h>
#include <math/scmat/local.h>
#include <chemistry/solvent/bem.h>

#define CLASSNAME BEMSolvent
#define PARENTS public DescribedClass
#define HAVE_KEYVAL_CTOR
//#include <util/state/statei.h>
#include <util/class/classi.h>
void *
BEMSolvent::_castdown(const ClassDesc*cd)
{
  void* casts[1];
  casts[0] = DescribedClass::_castdown(cd);
  return do_castdowns(casts,cd);
}

BEMSolvent::BEMSolvent(const RefKeyVal& keyval)
{
  vertex_area_ = 0;

  matrixkit_ = new LocalSCMatrixKit;
  
  debug_ = keyval->intvalue("debug");
  
  solute_ = keyval->describedclassvalue("solute");

  surf_ = keyval->describedclassvalue("surface");

  dielectric_constant_ = keyval->doublevalue("dielectric_constant");
  if (keyval->error() != KeyVal::OK) dielectric_constant_ = 78.0;
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
  TriangulatedSurfaceIntegrator triint(surf_);
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
BEMSolvent::done()
{
  surf_->clear();
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
  charges_to_surface_charge_density(charges);
  double charge = 0.0;
  
  // integrate over the surface
  TriangulatedSurfaceIntegrator triint(surf_);
  for (triint = 0; triint.update(); triint++) {
      int j0 = triint.vertex_number(0);
      int j1 = triint.vertex_number(1);
      int j2 = triint.vertex_number(2);
      double r = triint.r();
      double s = triint.s();
      double dA = triint.w();
      charge += dA * (charges[j0] * (1 - r - s)
                      + charges[j1] * r
                      + charges[j2] * s);
    }

  surface_charge_density_to_charges(charges);
  return charge;
}

// the passed enclosed_charge is determined by the called and
// might different from the enclosed charge computed by Gauss's
// law, which is stored as computed_enclosed_charge_
void
BEMSolvent::normalize_charge(double enclosed_charge, double* charges)
{
  charges_to_surface_charge_density(charges);

  double expected_charge = enclosed_charge
                         * (1.0/dielectric_constant_ - 1.0);

  double charge = 0.0;
  double area = 0.0;
  // integrate over the surface
  TriangulatedSurfaceIntegrator triint(surf_);
  for (triint = 0; triint.update(); triint++) {
      int j0 = triint.vertex_number(0);
      int j1 = triint.vertex_number(1);
      int j2 = triint.vertex_number(2);
      double r = triint.r();
      double s = triint.s();
      double dA = triint.w();
      charge += dA * (charges[j0] * (1 - r - s)
                      + charges[j1] * r
                      + charges[j2] * s);
      area += dA;
    }

  double charge_correction = (expected_charge - charge)/area;

  surface_charge_density_to_charges(charges);

  // add the correction to the charges
  for (triint = 0; triint.update(); triint++) {
      int j0 = triint.vertex_number(0);
      int j1 = triint.vertex_number(1);
      int j2 = triint.vertex_number(2);
      double r = triint.r();
      double s = triint.s();
      double dA = triint.w();
      charges[j0] += (1 - r - s) * dA * charge_correction;
      charges[j1] += r * dA * charge_correction;
      charges[j2] += s * dA * charge_correction;
    }

  if (debug_) {
      cout << node0 << indent
           << "BEMSolvent:normalize_charge:"
           << endl << indent
           << scprintf("  integrated surface charge = %20.15f", charge)
           << endl << indent
           << scprintf("  expected surface charge = %20.15f", expected_charge)
           << endl << indent
           << scprintf("  surface area = %20.15f", area) << endl;
    }
}

void
BEMSolvent::init_system_matrix()
{
  int i;
  int ntri = surf_->ntriangle();
  int n = ncharge();

  RefSCDimension d = new SCDimension(n);
  RefSCMatrix system_matrix(d,d,matrixkit());
  system_matrix.assign(0.0);

  // integrate over the surface
  double A = 0.0;
  double V = 0.0;
  TriangulatedSurfaceIntegrator triint(surf_);
  for (triint = 0; triint.update(); triint++) {
      V += triint.weight()*triint.dA()[2]*triint.current()->point()[2];
      const SCVector3& surfpv = triint.current()->point();
      int j0 = triint.vertex_number(0);
      int j1 = triint.vertex_number(1);
      int j2 = triint.vertex_number(2);
      double r = triint.r();
      double s = triint.s();
      double rs = 1 - r - s;
      double dA = triint.w();
      A += dA;
      double fdA = - f_ * dA;
      double rfdA = r * fdA;
      double sfdA = s * fdA;
      double rsfdA = rs * fdA;
      // loop thru all the vertices
      for (i = 0; i<n; i++) {
          RefVertex v = surf_->vertex(i);
          const SCVector3& pv = v->point();
          const SCVector3& nv = v->normal();
          SCVector3 diff(pv - surfpv);
          double normal_component = diff.dot(nv);
          double diff2 = diff.dot(diff);
          if (diff2 <= 1.0e-8) {
              // The self term must not be included here.  This
              // case shouldn't occur for the usual integrators
              // so abort.
              cerr << "BEMSolvent: integrator gave the self term" << endl;
              abort();
            }
          double denom = diff2*sqrt(diff2);
          double common_factor = normal_component/denom;
          system_matrix.set_element(i, j0,
                                    common_factor * rsfdA
                                    + system_matrix.get_element(i, j0));
          system_matrix.set_element(i, j1, common_factor * rfdA
                                    + system_matrix.get_element(i, j2));
          system_matrix.set_element(i, j2, common_factor * sfdA
                                    + system_matrix.get_element(i, j2));
        }
    }

  area_ = A;
  volume_ = V;

  cout << node0 << indent
       << scprintf("Solvent Accessible Surface: Area = %15.10f ", A)
       << scprintf("Volume = %15.10f", V) << endl;

  // Add I to the system matrix.
  system_matrix->shift_diagonal(1.0);

  //system_matrix->print("System Matrix");

  system_matrix_i_ = system_matrix.i();

  //system_matrix_i_->print("System Matrix Inverse");
}

void
BEMSolvent::compute_charges(double* efield_dot_normals, double* charges)
{
  if (system_matrix_i_.null()) {
      init_system_matrix();
    }

  // **this surface integral is only needed for computed_enclosed_charge_
  // **perhaps it should be avoided
  // integrate over the surface
  double A = 0.0;
  double efield_dot_normal = 0.0;
  TriangulatedSurfaceIntegrator triint(surf_);
  for (triint = 0; triint.update(); triint++) {
      const SCVector3& surfpv = triint.current()->point();
      int j0 = triint.vertex_number(0);
      int j1 = triint.vertex_number(1);
      int j2 = triint.vertex_number(2);
      double r = triint.r();
      double s = triint.s();
      double rs = 1 - r - s;
      double dA = triint.w();
      A += dA;

      efield_dot_normal += ( efield_dot_normals[j0] * r
                             +efield_dot_normals[j1] * s
                             +efield_dot_normals[j2] * rs) * dA;

      
    }

  computed_enclosed_charge_ = efield_dot_normal/(4.0*M_PI);

  if (debug_) {
      double computed_expected_charge = computed_enclosed_charge_
                                      * (1.0/dielectric_constant_ - 1.0);

      cout << node0 << indent
         << scprintf("BEMSolvent:compute_charges: Surface Area = %20.15f", A)
         << endl << indent
         << scprintf("BEMSolvent:compute_charges: encl nuc + elec q = %20.15f",
                     computed_enclosed_charge_)
         << endl << indent
         << scprintf("BEMSolvent:compute_charges: exp surface q = %20.15f",
                     computed_expected_charge) << endl;
    }

  RefSCVector edotn(system_matrix_i_.coldim(),matrixkit());
  edotn.assign(efield_dot_normals);
  //edotn.print("E dot normals");
  edotn.scale(f_);
  RefSCVector chrg = system_matrix_i_ * edotn;
  //chrg.print("Charges");
  chrg.convert(charges);

  surface_charge_density_to_charges(charges);
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

  // Precompute all of the integration points.
  TriangulatedSurfaceIntegrator itemplate(surf_);
  int n_integration_points = itemplate.n();
  TriangulatedSurfaceIntegrator *integration_points
      = new TriangulatedSurfaceIntegrator[n_integration_points];

  for (itemplate=0, i=0; i<n_integration_points; i++,itemplate++) {
      integration_points[i] = itemplate;
      integration_points[i].update();
    }

  double energy = 0.0;
  for (i=0; i<n_integration_points; i++) {
      const SCVector3& ipv = integration_points[i].current()->point();
      int iv0 = integration_points[i].vertex_number(0);
      int iv1 = integration_points[i].vertex_number(1);
      int iv2 = integration_points[i].vertex_number(2);
      double ir = integration_points[i].r();
      double is = integration_points[i].s();
      double irs = 1 - ir - is;
      double idA = integration_points[i].w();
      double icharge = (charge[iv0]*irs + charge[iv1]*ir + charge[iv2]*is)*idA;
      // interaction of the surface element with itself
      energy += 0.0; // is this good enough?
      // interaction of the surface element with all other surface elements
      for (j = 0; j<i; j++) {
          const SCVector3& jpv = integration_points[j].current()->point();
          int jv0 = integration_points[j].vertex_number(0);
          int jv1 = integration_points[j].vertex_number(1);
          int jv2 = integration_points[j].vertex_number(2);
          double jr = integration_points[j].r();
          double js = integration_points[j].s();
          double jrs = 1 - jr - js;
          double jdA = integration_points[j].w();
          double jcharge = (charge[jv0]*jrs
                            + charge[jv1]*jr
                            + charge[jv2]*js)*jdA;
          energy += icharge * jcharge/ipv.dist(jpv);
        }
    }

  delete[] integration_points;

  surface_charge_density_to_charges(charge);

  return energy;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
