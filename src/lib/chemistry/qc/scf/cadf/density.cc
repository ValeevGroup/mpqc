//
// density.cc
//
// Copyright (C) 2014 David Hollman
//
// Author: David Hollman
// Maintainer: DSH
// Created: May 20, 2014
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


#include <math/scmat/local.h>
#include <math/scmat/repl.h>
#include <math/scmat/blocked.h>
#include <math/scmat/offset.h>
#include <util/misc/scexception.h>

#include "cadfclhf.h"

using namespace sc;

//////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////

double
CADFCLHF::new_density()
{
  // copy current density into density diff and scale by -1.  later we'll
  // add the new density to this to get the density difference.
  cl_dens_diff_.assign(cl_dens_);
  cl_dens_diff_.scale(-1.0);

  //so_density(cl_dens_, 2.0);
  //============================================================//
  int i,j,k,kk;
  int me=scf_grp_->me();
  int nproc=scf_grp_->n();
  RefSCMatrix oso_vector = oso_scf_vector_;
  if (oso_vector.null()) {
    oso_vector = oso_eigenvectors();
  }
  RefDiagSCMatrix evals = current_evals_;
  RefSCMatrix vector = so_to_orthog_so().t() * oso_vector;
  oso_vector = 0;

  if (so_dimension()->equiv(oso_dimension())) {
    double core_evals_min = core_evals_.get_element(0);
    int ndocc = ndocc_[0];

    BlockedDiagSCMatrix *bevals = require_dynamic_cast<BlockedDiagSCMatrix*>(
        evals, "SCF::so_density: blocked evals");
    BlockedSCMatrix *bvec = require_dynamic_cast<BlockedSCMatrix*>(
        vector, "SCF::so_density: blocked vector");
    BlockedSymmSCMatrix *bd = require_dynamic_cast<BlockedSymmSCMatrix*>(
        cl_dens_, "SCF::so_density: blocked density");
    BlockedSymmSCMatrix *bs = require_dynamic_cast<BlockedSymmSCMatrix*>(
        overlap(), "SCF::so_density: overlap");
    RefSCMatrix vir = bvec->block(0);
    RefSymmSCMatrix dir = bd->block(0);
    RefDiagSCMatrix evir = bevals->block(0);

    //if(prev_evecs_.null()) {
    //  prev_evecs_ = core_evecs_;
    //  for(i = 0; i < ndocc; ++i) {
    //    prev_occ_.push_back(i);
    //  }
    //  prev_evals_ = core_evals_;
    //}
    BlockedDiagSCMatrix *bcevals = require_dynamic_cast<BlockedDiagSCMatrix*>(
        core_evals_, "SCF::so_density: blocked core evals");
    BlockedSCMatrix *bcvec = require_dynamic_cast<BlockedSCMatrix*>(
        core_evecs_, "SCF::so_density: blocked core evecs");
    //BlockedDiagSCMatrix *bpevals = require_dynamic_cast<BlockedDiagSCMatrix*>(
    //    prev_evals_, "SCF::so_density: blocked core evals");
    //BlockedSCMatrix *bpvec = require_dynamic_cast<BlockedSCMatrix*>(
    //    prev_evecs_, "SCF::so_density: blocked core evecs");

    RefSCMatrix cevecs = bcvec->block(0);
    RefDiagSCMatrix cevals = bcevals->block(0);
    //RefSCMatrix cevecs = bpvec->block(0);
    //RefDiagSCMatrix cevals = bpevals->block(0);

    RefSymmSCMatrix sir = bs->block(0);
    int n_SO = so_dimension()->blocks()->size(0);

    //current_evals_.print("Current eigenvalues");
    double *dens;
    if (dynamic_cast<LocalSymmSCMatrix*>(dir.pointer()))
      dens = dynamic_cast<LocalSymmSCMatrix*>(dir.pointer())->get_data();
    else if (dynamic_cast<ReplSymmSCMatrix*>(dir.pointer()))
      dens = dynamic_cast<ReplSymmSCMatrix*>(dir.pointer())->get_data();
    else
      abort();

    if(not match_orbitals_) {
      // Figure out which orbitals to use
      int col0 = -1, coln = -1;
      int icol = 0;
      while(icol < n_SO) {
        if(evir.get_element(icol) > core_evals_min) {
          col0 = icol;
          coln = icol + ndocc - 1;
          break;
        }
        ++icol;
      }
      if(col0 == -1 || coln > n_SO) {
        throw FeatureNotImplemented("Too few valid eigenvalues to construct a reasonable density matrix", __FILE__, __LINE__, class_desc());
      }
      if(col0 > 0) {
        ExEnv::out0() << indent << "Excluding " << col0 << " bad orbitals." << std::endl;
      }

      RefSCMatrix occbits = vir->get_subblock(0, n_SO-1, col0, coln);

      double **c;
      if (dynamic_cast<LocalSCMatrix*>(occbits.pointer()))
        c = dynamic_cast<LocalSCMatrix*>(occbits.pointer())->get_rows();
      else if (dynamic_cast<ReplSCMatrix*>(occbits.pointer()))
        c = dynamic_cast<ReplSCMatrix*>(occbits.pointer())->get_rows();
      else
        abort();

      int ij=0;
      for (i=0; i < n_SO; i++) {
        for (j=0; j <= i; j++, ij++) {
          if (ij%nproc != me)
            continue;

          double dv = 0;

          int kk=0;
          for (k=col0; k <= coln; k++, kk++)
            dv += c[i][kk]*c[j][kk];

          dens[ij] = dv;
        }
      }

      if (nproc > 1) scf_grp_->sum(dens, i_offset(n_SO));
    }

    //----------------------------------------------------------------------------//
    else {
      // match up orbitals with previous iteration
      // TODO distribute this work
      // for now, just map all of them
      // TODO track evecs across iterations by comparing to the previous iteration
      //if(orb_to_hcore_orb_.size() == 0) {

        Eigen::MatrixXd Dtmp(dir.n(), dir.n());
        Dtmp.setZero();

        double *Cmo_data;
        if (dynamic_cast<LocalSCMatrix*>(vir.pointer()))
          Cmo_data = dynamic_cast<LocalSCMatrix*>(vir.pointer())->get_data();
        else if (dynamic_cast<ReplSCMatrix*>(vir.pointer()))
          Cmo_data = dynamic_cast<ReplSCMatrix*>(vir.pointer())->get_data();
        else
          abort();
        double *Chcore_data;
        if (dynamic_cast<LocalSCMatrix*>(cevecs.pointer()))
          Chcore_data = dynamic_cast<LocalSCMatrix*>(cevecs.pointer())->get_data();
        else if (dynamic_cast<ReplSCMatrix*>(cevecs.pointer()))
          Chcore_data = dynamic_cast<ReplSCMatrix*>(cevecs.pointer())->get_data();
        else
          abort();

        Eigen::Map<RowMatrix> Cmo(Cmo_data, vir.nrow(), vir.ncol());
        Eigen::Map<RowMatrix> Chcore(Chcore_data, vir.nrow(), vir.ncol());
        if(vir.ncol() != cevecs.ncol()) {
          throw FeatureNotImplemented("Different numbers of MOs in core hamiltonian and current coefficients", __FILE__, __LINE__, class_desc());
        }
        if(vir.nrow() != cevecs.nrow()) {
          throw FeatureNotImplemented("Different numbers of SOs in core hamiltonian and current coefficients", __FILE__, __LINE__, class_desc());
        }

        orb_to_hcore_orb_.resize(vir.ncol());
        //Eigen::VectorXd hcore_cuts(vir.ncol());
        //Eigen::VectorXd negdiffnorms(vir.ncol());
        Eigen::VectorXd diffnorms(vir.ncol());
        Eigen::VectorXi indices(vir.ncol());
        indices = Eigen::VectorXi::LinSpaced(Eigen::Sequential, vir.ncol(), 0, vir.ncol() - 1);
        Eigen::VectorXi idx_sorted(vir.ncol());

        std::vector<int> occupied;
        int max_occupied = 0;

        const double hcore_homo_energy = cevals.get_element(ndocc);
        const double ecut = hcore_homo_energy + match_orbitals_max_homo_offset_;
        const int nmo = vir.ncol();
        //auto&& sign = [](double val) -> double { return val < 0 ? -1 : (val > 0 ? 1 : 0); };
        //auto&& cull = [&sign](double val) -> double {
        //  if(fabs(val) < 1e-3) {
        //    return 0;
        //  }
        //  else if(fabs(val) < 1.0) {
        //    return sign(val) * 1.0/(-log(fabs(val)));
        //  }
        //  else {
        //    return sign(val) * 10;
        //  }
        //  //if(fabs(val) < 1e-3) return 0;
        //  //else if(fabs(val) < 1e-2) return sign(val)*1;
        //  //else if(fabs(val) < 1e-1) return sign(val)*2;
        //  //else return sign(val)*3;
        //};
        //Eigen::MatrixXd Chcore_tmp = Chcore.unaryExpr(cull);
        //Eigen::MatrixXd Cmo_tmp = Cmo.unaryExpr(cull);

        Eigen::MatrixXd Sso(n_SO, n_SO);
        for(i = 0; i < n_SO; ++i) {
          for(j = 0; j <= i; ++j) {
            Sso(i,j) = Sso(j,i) = sir.get_element(i,j);
          }
        }
        Eigen::MatrixXd Smo_core(nmo, nmo);
        Smo_core = Cmo.transpose() * Sso * Chcore;

        //Eigen::VectorXd Chcore_tmp_norms = Chcore_tmp.colwise().norm();
        //Eigen::Eigen::VectorXd Cmo_tmp_norms = Cmo_tmp.colwise().norm();

        for(int imo = 0; imo < vir.ncol(); ++imo) {
          // Allow for sign changes
          //minpos = (Chcore.colwise() - Cmo.col(imo)).colwise().norm().minCoeff(minposcol);
          //diffnorms = (
          //    (Chcore.array().abs() > 1e-1).select(1, 0).colwise()
          //    - (Cmo.col(imo).array().abs() > 1e-1).select(Cmo.col(imo), 0)
          //).colwise().norm();
          //negdiffnorms = (
          //    (Chcore.array().abs() > 1e-1).select(Chcore, 0).colwise()
          //    + (Cmo.col(imo).array().abs() > 1e-1).select(Cmo.col(imo), 0)
          //).colwise().norm();
          //diffnorms = (Chcore_tmp.colwise() - Cmo_tmp.col(imo)).colwise().norm();
          //negdiffnorms = (Chcore_tmp.colwise() + Cmo_tmp.col(imo)).colwise().norm();
          //diffnorms = (Chcore.colwise() - Cmo.col(imo)).colwise().norm();
          //negdiffnorms = (Chcore.colwise() + Cmo.col(imo)).colwise().norm();
          //diffnorms = (diffnorms.array() < negdiffnorms.array()).select(diffnorms, negdiffnorms);
          //for(int iii = 0; iii < n_SO; ++iii) {
          //  diffnorms[iii] = std::min(diffnorms[iii], negdiffnorms[iii]);
          //}
          diffnorms = Smo_core.row(imo).cwiseAbs();

          bool can_be_occupied = false;

          //for(auto&& occ_idx : prev_occ_) {
          //  if(diffnorms[occ_idx] > match_orbitals_max_diff_) {
          //    can_be_occupied = true;
          //    break;
          //  }
          //}

          idx_sorted = indices;
          std::sort(idx_sorted.data(), idx_sorted.data() + vir.ncol(), [&diffnorms](const int& a, const int& b){
            //return diffnorms[a] < diffnorms[b];
            return diffnorms[a] > diffnorms[b];
          });
          //orb_to_hcore_orb_[imo] = idx_sorted[0];
          double hcore_cut;

          int ispot = 0;

          while(diffnorms[idx_sorted[ispot]] > match_orbitals_max_diff_ and ispot < nmo) {
            const double icut = cevals.get_element(idx_sorted[ispot]);
            if(icut < hcore_cut) hcore_cut = icut;

            if(evir.get_element(imo) > icut) {
              if(icut < ecut) {
                can_be_occupied = true;
              }
            }

            ++ispot;

          }
          //hcore_cuts[imo] = hcore_cut;

          //auto&& printorb = [&](int iorb, bool core) {
          //  printf("%7s=%d:%.7f\n", "Sym", iorb, diffnorms[iorb]);
          //  printf("%7s=%16.10f\n", "Ene", core ? cevals.get_element(iorb) : evir.get_element(iorb));
          //  printf("%7s=%7s\n", "Spin", "Alpha");
          //  printf("%7s=%7s\n", "Occup", "2.0000000000");
          //  for(int iao = 0; iao < n_SO; ++iao) {
          //    printf("%7d %16.10f\n", iao+1, core ? Chcore(iao, iorb) : Cmo(iao, iorb));
          //  }
          //};

          if(can_be_occupied and occupied.size() < ndocc and evir.get_element(imo) > hcore_cut) { //
          //if(can_be_occupied and imo < ndocc) {
            occupied.push_back(imo);
            max_occupied = imo;
          }
          else if(occupied.size() < ndocc) {
          //else if(imo < ndocc) {
            // TODO Reset the DIIS space when orbitals are skipped.
            ExEnv::out0() << indent << "Skipping bad orbital " << imo << "; max overlap was with hcore orbital "
                << idx_sorted[0] << ": overlap of " << diffnorms[idx_sorted[0]] << std::endl;
            //ExEnv::out0() << "MO "<< imo << ":" << std::endl << Cmo.col(imo).transpose() << std::endl;
            //ExEnv::out0() << "HCore MO "<< imo << ":" << std::endl << Chcore.col(imo).transpose() << std::endl;
            //ExEnv::out0() << "HCore MO "<< idx_sorted[0] << ":" << std::endl << Chcore.col(idx_sorted[0]).transpose() << std::endl;
            //ExEnv::out0() << indent << "   " << ispot-1 << " other orbitals considered, diff norm was " << diffnorms[idx_sorted[0]] << std::endl;
            //printorb(imo, false);
            //printorb(imo, true);
            //printorb(idx_sorted[0], true);
            //printorb(idx_sorted[1], true);
            //printorb(idx_sorted[2], true);
            //printorb(idx_sorted[3], true);
            //printorb(idx_sorted[4], true);
            //printorb(idx_sorted[5], true);
            //printorb(idx_sorted[6], true);
            //printorb(idx_sorted[7], true);
            //printorb(idx_sorted[8], true);
            //printorb(idx_sorted[9], true);
            //abort();
            //if(ispot > 1)
            //  ExEnv::out0() << indent << "  max diff norm considered was " << diffnorms[idx_sorted[ispot-1]] << std::endl;
            //else {
            //  ExEnv::out0() << indent << "  next best options were:"
            //      << std::endl << indent << "    "<< idx_sorted[1] << " with diff norm of " << diffnorms[idx_sorted[1]]
            //      << std::endl << indent << "    "<< idx_sorted[2] << " with diff norm of " << diffnorms[idx_sorted[2]]
            //      << std::endl << indent << "    "<< idx_sorted[3] << " with diff norm of " << diffnorms[idx_sorted[3]]
            //      << std::endl << indent << "    "<< idx_sorted[4] << " with diff norm of " << diffnorms[idx_sorted[4]]
            //      << std::endl << indent << "    "<< idx_sorted[5] << " with diff norm of " << diffnorms[idx_sorted[5]]
            //      << std::endl;
            //  ExEnv::out0() << indent << "  Diff with corresponding orbital (" << imo << ") was " << diffnorms[imo] << std::endl;
            //}
          }
        }

        //prev_occ_ = occupied;

        if(occupied.size() < ndocc) {
          throw SCException("Too few valid eigenvalues to construct a reasonable density matrix", __FILE__, __LINE__, class_desc());
        }

        if(match_orbitals_use_svd_ and not max_occupied == ndocc - 1) {
          Eigen::JacobiSVD<decltype(Smo_core)> s_svd(Smo_core.topRows(max_occupied), Eigen::ComputeThinU);
          auto& Umo_occ = s_svd.matrixU();
          Dtmp = Cmo * Umo_occ * Umo_occ.transpose() * Cmo.transpose();
        }
        else {
          for(auto idx : occupied) {
            Dtmp += Cmo.col(idx) * Cmo.col(idx).transpose();
          }
        }

      int ij=0;
      for (i=0; i < n_SO; i++) {
        for (j=0; j <= i; j++, ij++) {
          dens[ij] = Dtmp(i, j);
        }
      }

    }
    //----------------------------------------------------------------------------//
  }
  else{
    throw FeatureNotImplemented("Higher symmetry", __FILE__, __LINE__, class_desc());
  }

  //============================================================//
  cl_dens_.scale(2.0);
  if(!iter_log_.null()){
    iter_log_->log_density(cl_dens_);
  }

  cl_dens_diff_.accumulate(cl_dens_);

  Ref<SCElementScalarProduct> sp(new SCElementScalarProduct);
  cl_dens_diff_.element_op(sp.pointer(), cl_dens_diff_);

  double delta = sp->result();
  delta = sqrt(delta/i_offset(cl_dens_diff_.n()));

  //DEBUG_DELETE_THIS
  //static int istop = 0;
  //if(istop++) {
  //  return 0;
  //}
  //else {
  //  return delta;
  //}
  return delta;
}
