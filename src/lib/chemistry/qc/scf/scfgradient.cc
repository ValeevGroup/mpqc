//
// scfgradient.cc --- implementation of SCF::compute_gradient
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Edward Seidl <seidl@janed.com>
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

#include <iostream>

#include <stdexcept>

#include <util/misc/regtime.h>
#include <util/misc/formio.h>

#include <math/scmat/offset.h>
#include <math/scmat/local.h>

#include <chemistry/qc/basis/obint.h>

#include <chemistry/qc/scf/scf.h>
#include <chemistry/qc/scf/scflocal.h>

using namespace sc;

//////////////////////////////////////////////////////////////////////////////

static void
ob_gradient(const Ref<OneBodyDerivInt>& derint, double * gradient,
            const RefSymmSCMatrix& density, const Ref<GaussianBasisSet>& gbs_,
            const Ref<MessageGrp>& grp)
{
  int gsh=0;
  
  GaussianBasisSet& gbs = *gbs_.pointer();
  Molecule& mol = *gbs_->molecule().pointer();
  
  Ref<SCMatrixSubblockIter> diter =
    density->local_blocks(SCMatrixSubblockIter::Read);

  for (diter->begin(); diter->ready(); diter->next()) {
    SCMatrixBlock *dblk = diter->block();
    double *ddata;
    int istart, iend;
    int jstart, jend;
    int sub=0;

    ddata = get_tri_block(dblk, istart, iend, jstart, jend, sub);
    if (!ddata) {
      ExEnv::errn() << indent <<
        "ob_gradient: can't figure out what density block is\n";
      abort();
    }
    
    if (istart >= iend || jstart >= jend)
      continue;
    
    int ishstart = gbs.function_to_shell(istart);
    int ishend = (iend) ? gbs.function_to_shell(iend-1) : 0;

    int jshstart = gbs.function_to_shell(jstart);
    int jshend = (jend) ? gbs.function_to_shell(jend-1) : 0;
    
    for (int ish=ishstart; ish <= ishend; ish++) {
      GaussianShell& gsi = gbs(ish);
      
      int ist = gbs.shell_to_function(ish);
      int ien = ist + gsi.nfunction();

      for (int jsh=jshstart; jsh <= jshend; jsh++, gsh++) {
        if (jsh > ish)
          break;
        
        GaussianShell& gsj = gbs(jsh);

        int jst = gbs.shell_to_function(jsh);
        int jen = jst + gsj.nfunction();

        for (int x=0; x < mol.natom(); x++) {
          derint->compute_shell(ish,jsh,x);
          const double *buf = derint->buffer();

          int index=0;
          double dx=0, dy=0, dz=0;
          for (int i=ist; i < ien; i++) {
            for (int j=jst; j < jen; j++) {
              if (i < istart || i >= iend || j < jstart || j >= jend
                || j > i) {
                index += 3;
              } else {
                int doff = (sub) ? ij_offset(i,j) :
                                   ij_offset(i-istart,j-jstart);
                double denij = ddata[doff];
                if (j!=i) denij *= 2.0;
                dx += buf[index++] * denij;
                dy += buf[index++] * denij;
                dz += buf[index++] * denij;
              }
            }
          }

          gradient[x*3+0] += dx;
          gradient[x*3+1] += dy;
          gradient[x*3+2] += dz;
        }
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////////////

void
SCF::compute_gradient(const RefSCVector& gradient)
{
  Timer tim("compute gradient");
  int i;
  
  init_gradient();

  int n3=gradient.n();

  if (atom_basis()) {
    throw std::runtime_error("SCF::compute_gradient: atom_basis not supported");
  }

  // do the nuclear contribution
  tim.enter("nuc rep");
  
  double *g = new double[n3];
  nuclear_repulsion_energy_gradient(g);

  if (debug_) {
    gradient.assign(g);
    print_natom_3(gradient,"Nuclear Contribution to the Gradient:");
  }

  double *o = new double[n3];
  memset(o,0,sizeof(double)*gradient.n());

  // form overlap contribution
  tim.change("overlap gradient");
  RefSymmSCMatrix dens = lagrangian();
  Ref<OneBodyDerivInt> derint = integral()->overlap_deriv();
  ob_gradient(derint, o, dens, basis(), scf_grp_);

  scf_grp_->sum(o,n3);

  if (debug_) {
    gradient.assign(o);
    print_natom_3(gradient,"Overlap Contribution to the Gradient:");
  }

  for (i=0; i < n3; i++) g[i] += o[i];
  
  // other one electron contributions
  tim.change("one electron gradient");
  memset(o,0,sizeof(double)*gradient.n());
  dens = gradient_density();
  derint = integral()->hcore_deriv();
  ob_gradient(derint, o, dens, basis(), scf_grp_);

  scf_grp_->sum(o,n3);

  if (debug_) {
    gradient.assign(o);
    print_natom_3(gradient,"One-Electron Contribution to the Gradient:");
  }

  for (i=0; i < n3; i++) g[i] += o[i];
  
  dens=0;
  derint=0;
  
  // now calculate two electron contribution
  tim.change("two electron gradient");
  memset(o,0,sizeof(double)*gradient.n());
  two_body_deriv(o);
  tim.exit("two electron gradient");

  if (debug_) {
    gradient.assign(o);
    print_natom_3(gradient,"Two-Electron Contribution to the Gradient:");
  }

  for (i=0; i < n3; i++) g[i] += o[i];
  
  gradient.assign(g);
  delete[] g;
  delete[] o;

  if (debug_) {
    print_natom_3(gradient,"Total Gradient:");
  }
  
  done_gradient();
  tim.exit("compute gradient");
  //tim.print();
}

//////////////////////////////////////////////////////////////////////////////

void
SCF::compute_hessian(const RefSymmSCMatrix& hessian)
{
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
