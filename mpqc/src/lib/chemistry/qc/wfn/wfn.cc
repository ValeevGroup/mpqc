//
// wfn.cc
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

#ifdef __GNUC__
#pragma implementation
#endif

#include <stdlib.h>
#include <math.h>

#include <iostream>

#include <util/keyval/keyval.h>
#include <util/misc/timer.h>
#include <util/misc/formio.h>
#include <util/state/stateio.h>
#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/basis/symmint.h>
#include <chemistry/qc/intv3/intv3.h>

#include <chemistry/qc/wfn/wfn.h>

using namespace std;

#define CHECK_SYMMETRIZED_INTEGRALS 0

static ClassDesc Wavefunction_cd(
  typeid(Wavefunction),"Wavefunction",4,"public MolecularEnergy",
  0, 0, 0);

Wavefunction::Wavefunction(const Ref<KeyVal>&keyval):
  // this will let molecule be retrieved from basis
  // MolecularEnergy(new AggregateKeyVal(keyval,
  //                                     new PrefixKeyVal("basis", keyval))),
  MolecularEnergy(keyval),
  overlap_(this),
  overlap_eigvec_(this),
  overlap_isqrt_eigval_(this),
  overlap_sqrt_eigval_(this),
  hcore_(this),
  natural_orbitals_(this),
  natural_density_(this),
  bs_values(0),
  bsg_values(0)
{
  overlap_.compute() = 0;
  overlap_eigvec_.compute() = 0;
  overlap_isqrt_eigval_.compute() = 0;
  overlap_sqrt_eigval_.compute() = 0;
  hcore_.compute() = 0;
  natural_orbitals_.compute() = 0;
  natural_density_.compute() = 0;

  overlap_.computed() = 0;
  overlap_eigvec_.computed() = 0;
  overlap_isqrt_eigval_.computed() = 0;
  overlap_sqrt_eigval_.computed() = 0;
  hcore_.computed() = 0;
  natural_orbitals_.computed() = 0;
  natural_density_.computed() = 0;

  print_nao_ = keyval->booleanvalue("print_nao");
  print_npa_ = keyval->booleanvalue("print_npa");
  KeyValValuedouble lindep_tol_def(1.e-8);
  lindep_tol_ = keyval->doublevalue("lindep_tol", lindep_tol_def);
  symm_orthog_ = keyval->booleanvalue("symm_orthog",
                                      KeyValValueboolean(1));

  debug_ = keyval->intvalue("debug");

  gbs_ = require_dynamic_cast<GaussianBasisSet*>(
    keyval->describedclassvalue("basis").pointer(),
    "Wavefunction::Wavefunction\n"
    );

  integral_ << keyval->describedclassvalue("integrals");
  if (integral_.null())
    integral_ = new IntegralV3(gbs_);
  
  integral_->set_basis(gbs_);
  Ref<PetiteList> pl = integral_->petite_list();

  sodim_ = pl->SO_basisdim();
  aodim_ = pl->AO_basisdim();
  basiskit_ = gbs_->so_matrixkit();
}

Wavefunction::Wavefunction(StateIn&s):
  SavableState(s),
  MolecularEnergy(s),
  overlap_(this),
  overlap_isqrt_eigval_(this),
  overlap_sqrt_eigval_(this),
  overlap_eigvec_(this),
  hcore_(this),
  natural_orbitals_(this),
  natural_density_(this),
  bs_values(0),
  bsg_values(0)
{
  debug_ = 0;

  overlap_.compute() = 0;
  overlap_eigvec_.compute() = 0;
  overlap_isqrt_eigval_.compute() = 0;
  overlap_sqrt_eigval_.compute() = 0;
  hcore_.compute() = 0;
  natural_orbitals_.compute() = 0;
  natural_density_.compute() = 0;

  overlap_.computed() = 0;
  overlap_eigvec_.computed() = 0;
  overlap_isqrt_eigval_.computed() = 0;
  overlap_sqrt_eigval_.computed() = 0;
  hcore_.computed() = 0;
  natural_orbitals_.computed() = 0;
  natural_density_.computed() = 0;

  if (s.version(::class_desc<Wavefunction>()) >= 2) {
    s.get(print_nao_);
    s.get(print_npa_);
  }
  else {
    print_nao_ = 0;
    print_npa_ = 0;
  }

  if (s.version(::class_desc<Wavefunction>()) >= 3) {
    s.get(symm_orthog_);
  }
  else {
    symm_orthog_ = 1;
  }

  if (s.version(::class_desc<Wavefunction>()) >= 4) {
    s.get(lindep_tol_);
  }
  else {
    lindep_tol_ = 1.e-6;
  }

  gbs_ << SavableState::restore_state(s);
  integral_ << SavableState::restore_state(s);

  integral_->set_basis(gbs_);
  Ref<PetiteList> pl = integral_->petite_list();

  sodim_ = pl->SO_basisdim();
  aodim_ = pl->AO_basisdim();
  basiskit_ = gbs_->so_matrixkit();
}

void
Wavefunction::symmetry_changed()
{
  MolecularEnergy::symmetry_changed();

  Ref<PetiteList> pl = integral_->petite_list();
  sodim_ = pl->SO_basisdim();
  aodim_ = pl->SO_basisdim();
  osodim_ = 0;
  overlap_.result_noupdate() = 0;
  overlap_isqrt_eigval_.result_noupdate() = 0;
  overlap_sqrt_eigval_.result_noupdate() = 0;
  overlap_eigvec_.result_noupdate() = 0;
  basiskit_ = gbs_->so_matrixkit();
}

Wavefunction::~Wavefunction()
{
  if (bs_values) {
    delete[] bs_values;
    bs_values=0;
  }
  if (bsg_values) {
    delete[] bsg_values;
    bsg_values=0;
  }
}

void
Wavefunction::save_data_state(StateOut&s)
{
  MolecularEnergy::save_data_state(s);

  // overlap and hcore integrals are cheap so don't store them.
  // same goes for natural orbitals

  s.put(print_nao_);
  s.put(print_npa_);
  s.put(symm_orthog_);
  s.put(lindep_tol_);

  SavableState::save_state(gbs_.pointer(),s);
  SavableState::save_state(integral_.pointer(),s);
}

double
Wavefunction::charge()
{
  return molecule()->nuclear_charge() - nelectron();
}

RefSymmSCMatrix
Wavefunction::ao_density()
{
  return integral()->petite_list()->to_AO_basis(density());
}

RefSCMatrix
Wavefunction::natural_orbitals()
{
  if (!natural_orbitals_.computed()) {
      RefSymmSCMatrix dens = density();

      // transform the density into an orthogonal basis
      RefSCMatrix ortho = so_to_orthog_so();

      RefSymmSCMatrix densortho(oso_dimension(), basis_matrixkit());
      densortho.assign(0.0);
      densortho.accumulate_transform(so_to_orthog_so(),dens);

      RefSCMatrix natorb(oso_dimension(), oso_dimension(),
                         basis_matrixkit());
      RefDiagSCMatrix natden(oso_dimension(), basis_matrixkit());

      densortho.diagonalize(natden, natorb);

      // natorb is the ortho SO to NO basis transform
      // make natural_orbitals_ the SO to the NO basis transform
      natural_orbitals_ = so_to_orthog_so().t() * natorb;
      natural_density_ = natden;

      natural_orbitals_.computed() = 1;
      natural_density_.computed() = 1;
    }

  return natural_orbitals_.result_noupdate();
}

RefDiagSCMatrix
Wavefunction::natural_density()
{
  if (!natural_density_.computed()) {
      natural_orbitals();
    }

  return natural_density_.result_noupdate();
}

RefSymmSCMatrix
Wavefunction::overlap()
{
  if (!overlap_.computed()) {
    integral()->set_basis(gbs_);
    Ref<PetiteList> pl = integral()->petite_list();

#if ! CHECK_SYMMETRIZED_INTEGRALS
    // first form skeleton s matrix
    RefSymmSCMatrix s(basis()->basisdim(), basis()->matrixkit());
    Ref<SCElementOp> ov =
      new OneBodyIntOp(new SymmOneBodyIntIter(integral()->overlap(), pl));

    s.assign(0.0);
    s.element_op(ov);
    ov=0;

    // then symmetrize it
    RefSymmSCMatrix sb(so_dimension(), basis_matrixkit());
    pl->symmetrize(s,sb);

    overlap_ = sb;
#else
    ExEnv::out() << "Checking symmetrized overlap" << endl;

    RefSymmSCMatrix s(basis()->basisdim(), basis()->matrixkit());
    Ref<SCElementOp> ov =
      new OneBodyIntOp(new OneBodyIntIter(integral()->overlap()));

    s.assign(0.0);
    s.element_op(ov);
    ov=0;

    overlap_ = pl->to_SO_basis(s);

    //// use petite list to form saopl

    // form skeleton Hcore in AO basis
    RefSymmSCMatrix saopl(basis()->basisdim(), basis()->matrixkit());
    saopl.assign(0.0);

    ov = new OneBodyIntOp(new SymmOneBodyIntIter(integral_->overlap(), pl));
    saopl.element_op(ov);
    ov=0;

    // also symmetrize using pl->symmetrize():
    RefSymmSCMatrix spl(so_dimension(), basis_matrixkit());
    pl->symmetrize(saopl,spl);

    //// compare the answers

    int n = overlap_.result_noupdate().dim().n();
    int me = MessageGrp::get_default_messagegrp()->me();
    for (int i=0; i<n; i++) {
      for (int j=0; j<=i; j++) {
        double val1 = overlap_.result_noupdate().get_element(i,j);
        double val2 = spl.get_element(i,j);
        if (me == 0) {
          if (fabs(val1-val2) > 1.0e-6) {
            ExEnv::out() << "bad overlap vals for " << i << " " << j
                         << ": " << val1 << " " << val2 << endl;
          }
        }
      }
    }
#endif

    overlap_.computed() = 1;
  }

  return overlap_.result_noupdate();
}

RefSymmSCMatrix
Wavefunction::core_hamiltonian()
{
  if (!hcore_.computed()) {
    integral()->set_basis(gbs_);
    Ref<PetiteList> pl = integral()->petite_list();

#if ! CHECK_SYMMETRIZED_INTEGRALS
    // form skeleton Hcore in AO basis
    RefSymmSCMatrix hao(basis()->basisdim(), basis()->matrixkit());
    hao.assign(0.0);

    Ref<SCElementOp> hc =
      new OneBodyIntOp(new SymmOneBodyIntIter(integral_->kinetic(), pl));
    hao.element_op(hc);
    hc=0;

    Ref<OneBodyInt> nuc = integral_->nuclear();
    nuc->reinitialize();
    hc = new OneBodyIntOp(new SymmOneBodyIntIter(nuc, pl));
    hao.element_op(hc);
    hc=0;

    // now symmetrize Hso
    RefSymmSCMatrix h(so_dimension(), basis_matrixkit());
    pl->symmetrize(hao,h);

    hcore_ = h;
#else
    ExEnv::out() << "Checking symmetrized hcore" << endl;

    RefSymmSCMatrix hao(basis()->basisdim(), basis()->matrixkit());
    hao.assign(0.0);

    Ref<SCElementOp> hc =
      new OneBodyIntOp(new OneBodyIntIter(integral_->kinetic()));
    hao.element_op(hc);
    hc=0;

    Ref<OneBodyInt> nuc = integral_->nuclear();
    nuc->reinitialize();
    hc = new OneBodyIntOp(new OneBodyIntIter(nuc));
    hao.element_op(hc);
    hc=0;

    hcore_ = pl->to_SO_basis(hao);

    //// use petite list to form haopl

    // form skeleton Hcore in AO basis
    RefSymmSCMatrix haopl(basis()->basisdim(), basis()->matrixkit());
    haopl.assign(0.0);

    hc = new OneBodyIntOp(new SymmOneBodyIntIter(integral_->kinetic(), pl));
    haopl.element_op(hc);
    hc=0;

    nuc = integral_->nuclear();
    nuc->reinitialize();
    hc = new OneBodyIntOp(new SymmOneBodyIntIter(nuc, pl));
    haopl.element_op(hc);
    hc=0;

    // also symmetrize using pl->symmetrize():
    RefSymmSCMatrix h(so_dimension(), basis_matrixkit());
    pl->symmetrize(haopl,h);

    //// compare the answers

    int n = hcore_.result_noupdate().dim().n();
    int me = MessageGrp::get_default_messagegrp()->me();
    for (int i=0; i<n; i++) {
      for (int j=0; j<=i; j++) {
        double val1 = hcore_.result_noupdate().get_element(i,j);
        double val2 = h.get_element(i,j);
        if (me == 0) {
          if (fabs(val1-val2) > 1.0e-6) {
            ExEnv::out() << "bad hcore vals for " << i << " " << j
                         << ": " << val1 << " " << val2 << endl;
          }
        }
      }
    }
#endif

    hcore_.computed() = 1;
  }

  return hcore_.result_noupdate();
}

// computes intermediates needed to form orthogonalization matrices
// and their inverses.
void
Wavefunction::compute_overlap_eig()
{
  // first calculate S
  RefSymmSCMatrix M = overlap().copy();

  // Diagonalize M to get m and U
  RefSCMatrix U(M.dim(), M.dim(), M.kit());
  RefDiagSCMatrix m(M.dim(), M.kit());
  M.diagonalize(m,U);
  M = 0;

  Ref<SCElementMaxAbs> maxabsop = new SCElementMaxAbs;
  m.element_op(maxabsop.pointer());
  double maxabs = maxabsop->result();
  double s_tol = lindep_tol_ * maxabs;

  double minabs = maxabs;
  BlockedDiagSCMatrix *bm = dynamic_cast<BlockedDiagSCMatrix*>(m.pointer());
  if (bm == 0) {
      ExEnv::out() << node0 << "Wfn: orthog: expected blocked overlap" << endl;
    }
  int i, j;
  double *pm_sqrt = new double[bm->dim()->n()];
  double *pm_isqrt = new double[bm->dim()->n()];
  int *pm_index = new int[bm->dim()->n()];
  int *nfunc = new int[bm->nblocks()];
  int nfunctot = 0;
  int nlindep = 0;
  for (i=0; i<bm->nblocks(); i++) {
      nfunc[i] = 0;
      if (bm->block(i).null()) continue;
      int n = bm->block(i)->dim()->n();
      double *pm = new double[n];
      bm->block(i)->convert(pm);
      for (j=0; j<n; j++) {
          if (pm[j] > s_tol) {
              if (pm[j] < minabs) { minabs = pm[j]; }
              pm_sqrt[nfunctot] = sqrt(pm[j]);
              pm_isqrt[nfunctot] = 1.0/pm_sqrt[nfunctot];
              pm_index[nfunctot] = j;
              nfunc[i]++;
              nfunctot++;
            }
          else if (symm_orthog_) {
              pm_sqrt[nfunctot] = 0.0;
              pm_isqrt[nfunctot] = 0.0;
              pm_index[nfunctot] = j;
              nfunc[i]++;
              nfunctot++;
              nlindep++;
            }
          else {
              nlindep++;
            }
        }
      delete[] pm;
    }

  // make sure all nodes end up with exactly the same data
  MessageGrp::get_default_messagegrp()->bcast(nfunctot);
  MessageGrp::get_default_messagegrp()->bcast(nfunc, bm->nblocks());
  MessageGrp::get_default_messagegrp()->bcast(pm_sqrt,nfunctot);
  MessageGrp::get_default_messagegrp()->bcast(pm_isqrt,nfunctot);
  MessageGrp::get_default_messagegrp()->bcast(pm_index,nfunctot);

  if (symm_orthog_) {
      osodim_ = new SCDimension(bm->dim()->blocks(),
                                "ortho SO (symmetric)");
    }
  else {
      osodim_ = new SCDimension(nfunctot, bm->nblocks(),
                                nfunc, "ortho SO (canonical)");
      for (i=0; i<bm->nblocks(); i++) {
        osodim_->blocks()->set_subdim(i, new SCDimension(nfunc[i]));
      }
    }

  overlap_eigvec_ = basis_matrixkit()->matrix(sodim_, osodim_);
  if (symm_orthog_) {
      overlap_eigvec_.result_noupdate().assign(U);
    }
  else {
      BlockedSCMatrix *bev
          = dynamic_cast<BlockedSCMatrix*>(overlap_eigvec_.result_noupdate()
                                      .pointer());
      BlockedSCMatrix *bU
          = dynamic_cast<BlockedSCMatrix*>(U.pointer());
      int ifunc = 0;
      for (i=0; i<bev->nblocks(); i++) {
          if (bev->block(i).null()) continue;
          for (j=0; j<nfunc[i]; j++) {
              bev->block(i)->assign_column(
                  bU->block(i)->get_column(pm_index[ifunc]),j
                  );
              ifunc++;
            }
        }
    }
  overlap_eigvec_.computed() = 1;

  overlap_sqrt_eigval_ = basis_matrixkit()->diagmatrix(osodim_);
  overlap_sqrt_eigval_.result_noupdate()->assign(pm_sqrt);
  overlap_sqrt_eigval_.computed() = 1;
  overlap_isqrt_eigval_ = basis_matrixkit()->diagmatrix(osodim_);
  overlap_isqrt_eigval_.result_noupdate()->assign(pm_isqrt);
  overlap_isqrt_eigval_.computed() = 1;

  delete[] nfunc;
  delete[] pm_sqrt;
  delete[] pm_isqrt;
  delete[] pm_index;
  
  if (nlindep != 0) {
      ExEnv::out() << node0 << indent
           << "Wavefunction: orthog: WARNING: "
           << nlindep << " linearly dependent basis function"
           << (nlindep>1?"s":"")
           << endl
           << indent
           << scprintf("  tolerance = %12.8f * max", lindep_tol_)
           << endl;
    }

  ExEnv::out() << node0 << indent
       << scprintf("independent overlap eigenvalue max/min = %12.8f/%16.12f",
                   maxabs, minabs)
       << endl;
  max_overlap_eigval_ = maxabs;
  min_overlap_eigval_ = minabs;

  if (debug_ > 1) {
    overlap().print("S");
    overlap_eigvec_.result_noupdate().print("S eigvec");
    overlap_isqrt_eigval_.result_noupdate().print("s^(-1/2) eigvec");
    overlap_sqrt_eigval_.result_noupdate().print("s^(1/2) eigvec");
    so_to_orthog_so_inverse().print("oSO to SO");
    so_to_orthog_so().print("SO to oSO");
    (so_to_orthog_so()
     *so_to_orthog_so_inverse()).print("X*X^(-1)",ExEnv::out(),14);
  }
}

// returns the orthogonalization matrix
RefSCMatrix
Wavefunction::so_to_orthog_so()
{
  if (!overlap_eigvec_.computed()) compute_overlap_eig();

  RefSCMatrix trans;

  if (symm_orthog_) {
    trans = overlap_eigvec_.result_noupdate()
      * overlap_isqrt_eigval_.result_noupdate()
      * overlap_eigvec_.result_noupdate().t();
  }
  else {
    trans = overlap_isqrt_eigval_.result_noupdate()
      * overlap_eigvec_.result_noupdate().t();
  }

  return trans;
}

RefSCMatrix
Wavefunction::so_to_orthog_so_inverse()
{
  if (!overlap_eigvec_.computed()) compute_overlap_eig();

  RefSCMatrix trans;

  if (symm_orthog_) {
    trans = overlap_eigvec_.result_noupdate()
      * overlap_sqrt_eigval_.result_noupdate()
      * overlap_eigvec_.result_noupdate().t();
  }
  else {
    trans = overlap_eigvec_.result_noupdate()
      * overlap_sqrt_eigval_.result_noupdate();
  }

  return trans;
}

Ref<GaussianBasisSet>
Wavefunction::basis() const
{
  return gbs_;
}

Ref<Integral>
Wavefunction::integral()
{
  return integral_;
}

RefSCDimension
Wavefunction::so_dimension()
{
  return sodim_;
}

RefSCDimension
Wavefunction::ao_dimension()
{
  return aodim_;
}

RefSCDimension
Wavefunction::oso_dimension()
{
  if (!overlap_eigvec_.computed()) compute_overlap_eig();
  return osodim_;
}

Ref<SCMatrixKit>
Wavefunction::basis_matrixkit()
{
  return basiskit_;
}

void
Wavefunction::print(ostream&o) const
{
  MolecularEnergy::print(o);
  basis()->print_brief(o);
  // the other stuff is a wee bit too big to print
  if (print_nao_ || print_npa_) {
    tim_enter("NAO");
    RefSCMatrix naos = ((Wavefunction*)this)->nao();
    tim_exit("NAO");
    if (print_nao_) naos.print("NAO");
  }
}
    
RefSymmSCMatrix
Wavefunction::alpha_density()
{
  if (!spin_polarized()) {
    RefSymmSCMatrix result = density().copy();
    result.scale(0.5);
    return result;
  }
  ExEnv::err() << class_name() << "::alpha_density not implemented" << endl;
  abort();
  return 0;
}

RefSymmSCMatrix
Wavefunction::beta_density()
{
  if (!spin_polarized()) {
    RefSymmSCMatrix result = density().copy();
    result.scale(0.5);
    return result;
  }
  ExEnv::err() << class_name() << "::beta_density not implemented" << endl;
  abort();
  return 0;
}

RefSymmSCMatrix
Wavefunction::alpha_ao_density()
{
  return integral()->petite_list()->to_AO_basis(alpha_density());
}

RefSymmSCMatrix
Wavefunction::beta_ao_density()
{
  return integral()->petite_list()->to_AO_basis(beta_density());
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
