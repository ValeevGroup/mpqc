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
  typeid(Wavefunction),"Wavefunction",5,"public MolecularEnergy",
  0, 0, 0);

Wavefunction::Wavefunction(const Ref<KeyVal>&keyval):
  // this will let molecule be retrieved from basis
  // MolecularEnergy(new AggregateKeyVal(keyval,
  //                                     new PrefixKeyVal("basis", keyval))),
  MolecularEnergy(keyval),
  overlap_(this),
  hcore_(this),
  natural_orbitals_(this),
  natural_density_(this),
  bs_values(0),
  bsg_values(0)
{
  overlap_.compute() = 0;
  hcore_.compute() = 0;
  natural_orbitals_.compute() = 0;
  natural_density_.compute() = 0;

  overlap_.computed() = 0;
  hcore_.computed() = 0;
  natural_orbitals_.computed() = 0;
  natural_density_.computed() = 0;

  print_nao_ = keyval->booleanvalue("print_nao");
  print_npa_ = keyval->booleanvalue("print_npa");
  KeyValValuedouble lindep_tol_def(1.e-8);
  lindep_tol_ = keyval->doublevalue("lindep_tol", lindep_tol_def);
  if (keyval->exists("symm_orthog")) {
      ExEnv::out() << node0 << indent
                   << "WARNING: using obsolete \"symm_orthog\" keyword"
                   << endl;
      if (keyval->booleanvalue("symm_orthog")) {
        orthog_method_ = Symmetric;
      }
      else {
        orthog_method_ = Canonical;
      }
  }
  else {
    char *orthog_name = keyval->pcharvalue("orthog_method");
    if (!orthog_name) {
      orthog_method_ = Symmetric;
    }
    else if (::strcmp(orthog_name, "canonical") == 0) {
      orthog_method_ = Canonical;
    }
    else if (::strcmp(orthog_name, "symmetric") == 0) {
      orthog_method_ = Symmetric;
    }
    else if (::strcmp(orthog_name, "gramschmidt") == 0) {
      orthog_method_ = GramSchmidt;
    }
    else {
      ExEnv::err() << "ERROR: bad orthog_method: \""
                   << orthog_name << "\"" << endl;
      abort();
    }
    delete[] orthog_name;
  }

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
  hcore_(this),
  natural_orbitals_(this),
  natural_density_(this),
  bs_values(0),
  bsg_values(0)
{
  debug_ = 0;

  overlap_.compute() = 0;
  hcore_.compute() = 0;
  natural_orbitals_.compute() = 0;
  natural_density_.compute() = 0;

  overlap_.computed() = 0;
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

  if (s.version(::class_desc<Wavefunction>()) >= 5) {
    int orthog_enum;
    s.get(orthog_enum);
    orthog_method_ = (OrthogMethod) orthog_enum;
  }
  else if (s.version(::class_desc<Wavefunction>()) >= 3) {
    int symm_orthog;
    s.get(symm_orthog);
    if (symm_orthog) orthog_method_ = Symmetric;
    else orthog_method_ = Canonical;
  }
  else {
    orthog_method_ = Symmetric;
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
  aodim_ = pl->AO_basisdim();
  osodim_ = 0;
  overlap_.result_noupdate() = 0;
  orthog_trans_ = 0;
  orthog_trans_inverse_ = 0;
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
  s.put((int)orthog_method_);
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

    if (debug_ > 1) s.print("AO skeleton overlap");

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
Wavefunction::compute_overlap_eig(RefSCMatrix& overlap_eigvec,
                    RefDiagSCMatrix& overlap_isqrt_eigval,
                    RefDiagSCMatrix& overlap_sqrt_eigval)
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
          else if (orthog_method_ == Symmetric) {
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

  if (nlindep > 0 && orthog_method_ == Symmetric) {
    ExEnv::out() << node0 << indent
                 << "WARNING: " << nlindep
                 << " basis function"
                 << (sodim_.n()-osodim_.n()>1?"s":"")
                 << " ignored in symmetric orthogonalization."
                 << endl;
  }

  // make sure all nodes end up with exactly the same data
  MessageGrp::get_default_messagegrp()->bcast(nfunctot);
  MessageGrp::get_default_messagegrp()->bcast(nfunc, bm->nblocks());
  MessageGrp::get_default_messagegrp()->bcast(pm_sqrt,nfunctot);
  MessageGrp::get_default_messagegrp()->bcast(pm_isqrt,nfunctot);
  MessageGrp::get_default_messagegrp()->bcast(pm_index,nfunctot);

  if (orthog_method_ == Symmetric) {
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

  overlap_eigvec = basis_matrixkit()->matrix(sodim_, osodim_);
  if (orthog_method_ == Symmetric) {
      overlap_eigvec.assign(U);
    }
  else {
      BlockedSCMatrix *bev
          = dynamic_cast<BlockedSCMatrix*>(overlap_eigvec.pointer());
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

  overlap_sqrt_eigval = basis_matrixkit()->diagmatrix(osodim_);
  overlap_sqrt_eigval->assign(pm_sqrt);
  overlap_isqrt_eigval = basis_matrixkit()->diagmatrix(osodim_);
  overlap_isqrt_eigval->assign(pm_isqrt);

  delete[] nfunc;
  delete[] pm_sqrt;
  delete[] pm_isqrt;
  delete[] pm_index;
  
  max_orthog_res_ = maxabs;
  min_orthog_res_ = minabs;

  if (debug_ > 1) {
    overlap().print("S");
    overlap_eigvec.print("S eigvec");
    overlap_isqrt_eigval.print("s^(-1/2) eigvec");
    overlap_sqrt_eigval.print("s^(1/2) eigvec");
  }
}

void
Wavefunction::compute_symmetric_orthog()
{
  RefSCMatrix overlap_eigvec;
  RefDiagSCMatrix overlap_isqrt_eigval;
  RefDiagSCMatrix overlap_sqrt_eigval;
  compute_overlap_eig(overlap_eigvec,
                      overlap_isqrt_eigval,
                      overlap_sqrt_eigval);

  orthog_trans_ = overlap_eigvec
    * overlap_isqrt_eigval
    * overlap_eigvec.t();
  orthog_trans_inverse_ = overlap_eigvec
    * overlap_sqrt_eigval
    * overlap_eigvec.t();
}

void
Wavefunction::compute_canonical_orthog()
{
  RefSCMatrix overlap_eigvec;
  RefDiagSCMatrix overlap_isqrt_eigval;
  RefDiagSCMatrix overlap_sqrt_eigval;
  compute_overlap_eig(overlap_eigvec,
                      overlap_isqrt_eigval,
                      overlap_sqrt_eigval);

  orthog_trans_ = overlap_isqrt_eigval * overlap_eigvec.t();
  orthog_trans_inverse_ = overlap_eigvec * overlap_sqrt_eigval;
}

void
Wavefunction::compute_gs_orthog()
{
  // Orthogonalize each subblock of the overlap.
  max_orthog_res_ = 1.0;
  min_orthog_res_ = 1.0;
  BlockedSymmSCMatrix *S
    = dynamic_cast<BlockedSymmSCMatrix *>(overlap().pointer());
  int nblock = S->nblocks();
  Ref<BlockedSCMatrixKit> kit
    = dynamic_cast<BlockedSCMatrixKit*>(S->kit().pointer());
  Ref<SCMatrixKit> subkit = kit->subkit();
  RefSCMatrix *blockorthogs = new RefSCMatrix[nblock];
  int *nblockorthogs = new int[nblock];
  int northog = 0;
  for (int i=0; i<nblock; i++) {
    RefSymmSCMatrix Sblock = S->block(i);
    if (Sblock.null()) {
      blockorthogs[i] = 0;
      nblockorthogs[i] = 0;
      continue;
      }
    RefSCDimension dim = Sblock->dim();
    RefSCMatrix blockorthog(dim,dim,subkit);
    blockorthog->unit();
    double res;
    int nblockorthog = blockorthog->schmidt_orthog_tol(Sblock, lindep_tol_,
                                                       &res);
    if (res < min_orthog_res_) min_orthog_res_ = res;
    blockorthogs[i] = blockorthog;
    nblockorthogs[i] = nblockorthog;
    northog += nblockorthog;
  }

  // Construct the orthog SO basis SCDimension object.
  Ref<SCBlockInfo> blockinfo
    = new SCBlockInfo(northog, nblock, nblockorthogs);
  for (int i=0; i<nblock; i++) {
    blockinfo->set_subdim(i, new SCDimension(nblockorthogs[i]));
  }
  osodim_ = new SCDimension(blockinfo, "ortho SO (Gram-Schmidt)");

  // Replace each blockorthog by a matrix with only linear independent columns
  for (int i=0; i<nblock; i++) {
    if (nblockorthogs[i] == 0) continue;
    RefSCMatrix old_blockorthog = blockorthogs[i];
    blockorthogs[i] = subkit->matrix(sodim_->blocks()->subdim(i),
                                     osodim_->blocks()->subdim(i));
    blockorthogs[i].assign_subblock(old_blockorthog,
                                    0, sodim_->blocks()->subdim(i).n()-1,
                                    0, osodim_->blocks()->subdim(i).n()-1);
  }

  // Compute the inverse of each orthogonalization block.
  RefSCMatrix *inverse_blockorthogs = new RefSCMatrix[nblock];
  for (int i=0; i<nblock; i++) {
    if (nblockorthogs[i] == 0) {
      inverse_blockorthogs[i] = 0;
      }
    else {
      inverse_blockorthogs[i] = blockorthogs[i].gi();
      }
  }

  // Construct the complete transformation matrices
  orthog_trans_ = basis_matrixkit()->matrix(sodim_, osodim_);
  orthog_trans_inverse_ = basis_matrixkit()->matrix(osodim_, sodim_);
  orthog_trans_.assign(0.0);
  orthog_trans_inverse_.assign(0.0);
  BlockedSCMatrix *X
    = dynamic_cast<BlockedSCMatrix*>(orthog_trans_.pointer());
  BlockedSCMatrix *Xi
    = dynamic_cast<BlockedSCMatrix*>(orthog_trans_inverse_.pointer());
  for (int i=0; i<nblock; i++) {
    if (nblockorthogs[i] == 0) continue;
    int nrow = blockorthogs[i].rowdim().n();
    int ncol = blockorthogs[i].coldim().n();
    X->block(i).assign_subblock(blockorthogs[i],
                                0, nrow-1, 0, ncol-1,
                                0, 0);
    Xi->block(i).assign_subblock(inverse_blockorthogs[i],
                                 0, ncol-1, 0, nrow-1,
                                 0, 0);
  }
  orthog_trans_ = orthog_trans_.t();
  orthog_trans_inverse_ = orthog_trans_inverse_.t();

  delete[] blockorthogs;
  delete[] inverse_blockorthogs;
  delete[] nblockorthogs;
}

void
Wavefunction::compute_orthog_trans()
{
  switch(orthog_method_) {
  case GramSchmidt:
    ExEnv::out() << node0 << indent
                 << "Using Gram-Schmidt orthogonalization."
                 << endl;
    compute_gs_orthog();
    break;
  case Symmetric:
    compute_symmetric_orthog();
    ExEnv::out() << node0 << indent
                 << "Using symmetric orthogonalization."
                 << endl;
    break;
  case Canonical:
    compute_canonical_orthog();
    ExEnv::out() << node0 << indent
                 << "Using canonical orthogonalization."
                 << endl;
    break;
  default:
    ExEnv::out() << "Wavefunction::compute_orthog_trans(): bad orthog method"
                 << endl;
    abort();
  }

  ExEnv::out() << node0 << indent
               << "n(SO):        ";
  for (int i=0; i<sodim_->blocks()->nblock(); i++) {
    ExEnv::out() << node0 << scprintf(" %5d", sodim_->blocks()->size(i));
  }
  ExEnv::out() << node0 << endl;

  if (sodim_.n() != osodim_.n()) {
    ExEnv::out() << node0 << indent
                 << "n(orthog SO): ";
    for (int i=0; i<osodim_->blocks()->nblock(); i++) {
      ExEnv::out() << node0 << scprintf(" %5d", osodim_->blocks()->size(i));
      }
    ExEnv::out() << node0 << endl;

    ExEnv::out() << node0 << indent
                 << "WARNING: " << sodim_.n() - osodim_.n()
                 << " basis function"
                 << (sodim_.n()-osodim_.n()>1?"s":"")
                 << " discarded."
                 << endl;
    }
  ExEnv::out() << node0 << indent
               << "Maximum orthogonalization residual = "
               << max_orthog_res_ << endl
               << node0 << indent
               << "Minimum orthogonalization residual = "
               << min_orthog_res_ << endl;

  if (debug_ > 0) {
    osodim_.print();
    if (debug_ > 1) {
      orthog_trans_.print("SO to OSO");
      orthog_trans_inverse_.print("SO to OSO inverse");
      (orthog_trans_*overlap()
       *orthog_trans_.t()).print("X*S*X'",ExEnv::out(),14);
      (orthog_trans_inverse_.t()*overlap().gi()
       *orthog_trans_inverse_).print("X'^(-1)*S^(-1)*X^(-1)",
                                     ExEnv::out(),14);
      (orthog_trans_
       *orthog_trans_inverse_).print("X*X^(-1)",ExEnv::out(),14);
    }
  }
}

// returns the orthogonalization matrix
RefSCMatrix
Wavefunction::so_to_orthog_so()
{
  if (orthog_trans_.null()) {
    compute_orthog_trans();
  }
  return orthog_trans_;
}

RefSCMatrix
Wavefunction::so_to_orthog_so_inverse()
{
  if (orthog_trans_inverse_.null()) {
    compute_orthog_trans();
  }
  return orthog_trans_inverse_;
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
  if (osodim_.null()) compute_orthog_trans();
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

void
Wavefunction::obsolete()
{
  osodim_ = 0;
  orthog_trans_ = 0;
  orthog_trans_inverse_ = 0;

  MolecularEnergy::obsolete();
}

void
Wavefunction::copy_orthog_info(const Ref<Wavefunction>&wfn)
{
  if (orthog_trans_.nonnull() || orthog_trans_inverse_.nonnull()) {
    ExEnv::err() << "WARNING: Wavefunction: orthogonalization info changing"
                 << endl;
  }
  orthog_trans_ = wfn->so_to_orthog_so().copy();
  orthog_trans_inverse_ = wfn->so_to_orthog_so_inverse().copy();
  lindep_tol_ = wfn->lindep_tol_;
  orthog_method_ = wfn->orthog_method_;
  osodim_ = wfn->osodim_;
  min_orthog_res_ = wfn->min_orthog_res_;
  max_orthog_res_ = wfn->max_orthog_res_;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
