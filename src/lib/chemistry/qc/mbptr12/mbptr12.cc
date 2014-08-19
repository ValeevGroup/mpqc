//
// mbptr12.cc
//
// Copyright (C) 2001 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
// Maintainer: EV
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
#include <sstream>
#include <numeric>
#include <cmath>
#include <cassert>
#include <util/misc/string.h>
#include <util/misc/scexception.h>
#include <util/misc/formio.h>
#include <util/misc/exenv.h>
#include <util/state/stateio.h>
#include <math/scmat/blocked.h>
#include <chemistry/qc/mbptr12/mbptr12.h>
#include <chemistry/qc/mbptr12/r12_amps.h>
#include <util/misc/print.h>
#include <chemistry/qc/scf/clscf.h>
#include <chemistry/qc/scf/hsosscf.h>
#include <chemistry/qc/wfn/orbitalspace_utils.h>

using namespace std;
using namespace sc;


/*--------------------------------
  MBPT2_R12
 --------------------------------*/

static ClassDesc MBPT2_R12_cd(
  typeid(MBPT2_R12),"MBPT2_R12",9,"public MBPT2",
  0, create<MBPT2_R12>, create<MBPT2_R12>);

MBPT2_R12::MBPT2_R12(StateIn& s):
  MBPT2(s)
{
#if 0
  if (s.version(::class_desc<MBPT2_R12>()) < 9)
    throw InputError("Cannot use MBPT2_R12 class prior to version 9",__FILE__,__LINE__);
#endif

  r12eval_ << SavableState::restore_state(s);
  r12world_ << SavableState::restore_state(s);
  r12a_energy_ << SavableState::restore_state(s);
  r12ap_energy_ << SavableState::restore_state(s);
  r12app_energy_ << SavableState::restore_state(s);
  r12b_energy_ << SavableState::restore_state(s);
  r12c_energy_ << SavableState::restore_state(s);

  if((r12world()->r12tech()->ansatz()->orbital_product_GG()==R12Technology::OrbProdGG_pq) ||
     (r12world()->r12tech()->ansatz()->orbital_product_gg()==R12Technology::OrbProdgg_pq)) {
    throw InputError("MBPT2_R12::MBPT2_R12 -- pq Ansatz not allowed",__FILE__,__LINE__);
  }

  s.get(mp2_corr_energy_);
  s.get(cabs_singles_);
  s.get(cabs_singles_energy_);
  s.get(occ_orbs_);
  s.get(uocc_orbs_);
  s.get(uocc_orbs_pno_truncate_threshold_);

  twopdm_grid_ << SavableState::restore_state(s);
  s.get(plot_pair_function_[0]);
  s.get(plot_pair_function_[1]);
}

MBPT2_R12::MBPT2_R12(const Ref<KeyVal>& keyval):
  MBPT2(keyval)
{
  // MP2-R12 calculations can use virtual orbitals expanded in a separate basis set
  Ref<GaussianBasisSet> bs_vir = require_dynamic_cast<GaussianBasisSet*>(
      keyval->describedclassvalue("vir_basis").pointer(),
      "MBPT2_R12::MBPT2_R12\n"
      );
  Ref<OrbitalSpace> vbs;
  if (bs_vir) {
    int nlindep_vir = -1; // will compute the number of linear dependencies automatically
    vbs = orthogonalize("e(sym)", "VBS", bs_vir,
                        integral(),
                        this->orthog_method(), this->lindep_tol(),
                        nlindep_vir);
  }

  KeyValValuestring def_occ_orbitals("canonical");
  occ_orbs_ = keyval->stringvalue("occ_orbitals", def_occ_orbitals);
  if (occ_orbs_ != "pipek-mezey"
      && occ_orbs_ != "canonical") {
      throw InputError("invalid keyword value",
                       __FILE__, __LINE__, "occ_orbitals", occ_orbs_.c_str(),
                       this->class_desc());
    }
  if (occ_orbs_ != "canonical" && this->ref()->molecule()->point_group()->order() != 1)
    throw InputError("localized orbitals can only be used in C1 symmetry",
                     __FILE__, __LINE__, "occ_orbitals", occ_orbs_.c_str(), this->class_desc());
  if (occ_orbs_ != "canonical" && this->ref()->spin_polarized() == true)
    throw InputError("localized occupied orbitals can only be used for closed-shell molecules",
                     __FILE__, __LINE__, "occ_orbitals", occ_orbs_.c_str(), this->class_desc());

  KeyValValuestring def_uocc_orbitals("canonical");
  uocc_orbs_ = keyval->stringvalue("uocc_orbitals", def_uocc_orbitals);
  if (uocc_orbs_ != "pno"
      && uocc_orbs_ != "pno-r12"
      && uocc_orbs_ != "canonical") {
      throw InputError("invalid keyword value",
                       __FILE__, __LINE__, "uocc_orbitals", uocc_orbs_.c_str(),
                       this->class_desc());
    }
  if (uocc_orbs_ != "canonical" && this->ref()->molecule()->point_group()->order() != 1)
    throw InputError("non-canonical unoccupied orbitals can only be used in C1 symmetry",
                     __FILE__, __LINE__, "uocc_orbitals", uocc_orbs_.c_str(), this->class_desc());
  if (uocc_orbs_ != "canonical" && this->ref()->spin_polarized() == true)
    throw InputError("non-canonical unoccupied orbitals can only be used for closed-shell molecules",
                     __FILE__, __LINE__, "uocc_orbitals", uocc_orbs_.c_str(), this->class_desc());
  if (uocc_orbs_ == "pno" || uocc_orbs_ == "pno-r12") {
    uocc_orbs_pno_truncate_threshold_ = keyval->doublevalue("pno_truncate_threshold", KeyValValuedouble(1e-10));
  }

  // if world not given, make this the center of a new World
  Ref<WavefunctionWorld> world; world << keyval->describedclassvalue("world", KeyValValueRefDescribedClass(0));
  if (world.null())
    world = new WavefunctionWorld(keyval);
  if (world.null())
    throw InputError("MBPT2_R12 requires a WavefunctionWorld; input did not specify it, neither could it be constructed",
                     __FILE__, __LINE__, "world");
  if (world->wfn() == 0) world->set_wfn(this);

  const bool spin_restricted = false;   // do not use spin-restricted orbitals -> for ROHF use semicanonical orbitals
  Ref<RefWavefunction> refinfo = new SD_RefWavefunction(world, ref(), spin_restricted,
                                                        nfzcore(), nfzvirt(),
                                                        vbs, occ_orbs_);
  r12world_ = new R12WavefunctionWorld(keyval, refinfo);
  Ref<R12Technology> r12tech = r12world_->r12tech();

  Ref<R12Technology::NullCorrelationFactor> null_cf; null_cf << r12tech->corrfactor();
  const bool default_cabs_singles = null_cf ? false : true;
  cabs_singles_ = keyval->booleanvalue("cabs_singles",KeyValValueboolean(default_cabs_singles));

  twopdm_grid_ = require_dynamic_cast<TwoBodyGrid*>(
                   keyval->describedclassvalue("twopdm_grid").pointer(),
                   "MBPT2_R12::MBPT2_R12\n"
                 );
  if (twopdm_grid_) {
    plot_pair_function_[0] = static_cast<unsigned int>(keyval->intvalue("plot_pair_function", 0,
									KeyValValueint(0)
									));
    plot_pair_function_[1] = static_cast<unsigned int>(keyval->intvalue("plot_pair_function", 1,
									KeyValValueint(0)
									));
  }

  r12eval_ = 0;
  r12a_energy_ = 0;
  r12ap_energy_ = 0;
  r12app_energy_ = 0;
  r12b_energy_ = 0;
  r12c_energy_ = 0;
  mp2_corr_energy_ = 0.0;
  cabs_singles_energy_ = 0.0;

  this->set_desired_value_accuracy(desired_value_accuracy());

  if (bs_vir) { // create AO space for VBS basis, if needed
    Ref<OrbitalSpaceRegistry> idxreg = this->r12world()->world()->tfactory()->orbital_registry();
    Ref<AOSpaceRegistry> aoidxreg = this->r12world()->world()->tfactory()->ao_registry();
    Ref<Integral> localints = this->r12world()->refwfn()->integral()->clone();
    Ref<OrbitalSpace> nu = new AtomicOrbitalSpace("nu", "VBS(AO)", bs_vir, localints);
    idxreg->add(make_keyspace_pair(nu));
    aoidxreg->add(nu->basis(),nu);
  }
}

MBPT2_R12::~MBPT2_R12()
{
}

void
MBPT2_R12::save_data_state(StateOut& s)
{
  MBPT2::save_data_state(s);
  SavableState::save_state(r12eval_.pointer(),s);
  SavableState::save_state(r12world_.pointer(),s);
  SavableState::save_state(r12a_energy_.pointer(),s);
  SavableState::save_state(r12ap_energy_.pointer(),s);
  SavableState::save_state(r12app_energy_.pointer(),s);
  SavableState::save_state(r12b_energy_.pointer(),s);
  SavableState::save_state(r12c_energy_.pointer(),s);

  s.put(mp2_corr_energy_);
  s.put(cabs_singles_);
  s.put(cabs_singles_energy_);
  s.put(occ_orbs_);
  s.put(uocc_orbs_);
  s.put(uocc_orbs_pno_truncate_threshold_);

  SavableState::save_state(twopdm_grid_.pointer(),s);
  s.put((int)plot_pair_function_[0]);
  s.put((int)plot_pair_function_[1]);
}

void
MBPT2_R12::print(ostream&o) const
{
  o << indent << "MBPT2_R12:" << endl;
  o << incindent;

  o << indent << "Occupied orbitals     : " << occ_orbs_ << endl;
  o << indent << "Unoccupied orbitals   : " << uocc_orbs_;
  if (uocc_orbs_ == "pno" || uocc_orbs_ == "pno-r12") {
    o << " (truncate_threshold = " << scprintf("%10.5e",uocc_orbs_pno_truncate_threshold_) << ")";
  }
  o << endl;

  o << indent << "Include CABS singles? : " << (cabs_singles_ ? "true" : "false") << endl;
  if (cabs_singles_) {
    o << indent << "  E(CABS singles) = " << scprintf("%25.15lf", cabs_singles_energy_)
                                          << endl;
  }

  r12world()->print(o);

  MBPT2::print(o);
  o << decindent;
}

void
MBPT2_R12::set_desired_value_accuracy(double acc)
{
  Function::set_desired_value_accuracy(acc);
  if (ref()->desired_value_accuracy_set_to_default()) {
    // reference should be computed to higher accuracy
    const double ref_acc = acc * ref_to_mp2r12_acc();
    ref()->set_desired_value_accuracy(ref_acc);
  }
}

//////////////////////////////////////////////////////////////////////////////

RefSymmSCMatrix
MBPT2_R12::density()
{
  ExEnv::out0() << "MBPT2_R12::density() is not implemented" << endl;
  abort();
  return 0;
}

//////////////////////////////////////////////////////////////////////////////

void
MBPT2_R12::compute()
{
  ExEnv::out0() << "switched off integral check in mbptr12.cc/compute" << endl;

  compute_energy_();

  if (r12eval()->compute_1rdm()) {
#if defined(HAVE_MPQC3_RUNTIME)
    ExEnv::out0() << indent << "Compute MP2-F12 one-electron properties" << std::endl;
    r12eval()->compute_TA_mp2f12_1rdm();
#else
    throw FeatureNotImplemented("F12 property code requires TiledArray. Recompile MPQC with TA support.",
                                __FILE__, __LINE__);
#endif
  }

#define TESTING_PSV 1
#if TESTING_PSV
  if (uocc_orbs_ == "pno" ||
      uocc_orbs_ == "pno-r12") {

    std::vector<double> r12_pair_energies[NSpinCases2];
    std::vector<double> mp2_pair_energies[NSpinCases2];

//    r12world_->initialize();
//    if (r12eval_.null()) {
//      // since r12intevalinfo uses this class' KeyVal to initialize, dynamic is set automatically
//      r12world_->world()->print_percent(print_percent_);
//      r12eval_ = new R12IntEval(r12world_);
//      r12eval_->debug(debug_);
//    }

    for(int s=0; s<2; ++s) {
      SpinCase2 S = static_cast<SpinCase2>(s);
      SpinCase1 spin1 = case1(S);
      SpinCase1 spin2 = case2(S);
      Ref<OrbitalSpace> occ1 = r12world()->refwfn()->occ(spin1);
      Ref<OrbitalSpace> occ2 = r12world()->refwfn()->occ(spin2);
      Ref<OrbitalSpace> occ1_act = r12eval()->ggspace(spin1);
      Ref<OrbitalSpace> occ2_act = r12eval()->ggspace(spin2);
      const unsigned int nocc = occ1->rank();  // assuming RHF at the moment

      // make PSVs
      const bool geminal_project = (uocc_orbs_ == "pno-r12");
      std::vector< std::pair<RefSCMatrix,RefSCMatrix> > psv_set = this->mp1_pno(S, geminal_project, uocc_orbs_pno_truncate_threshold_);

      // for each pair feed the PSV orbitals to R12WavefunctionWorld, compute, and extract only this pair's contribution
      SpinMOPairIter iter_ij(occ1_act->rank(),occ2_act->rank(),S);
      const unsigned int nij = iter_ij.nij();
      unsigned int ij = 0;
      for(iter_ij.start(); iter_ij; iter_ij.next(), ++ij) {

        const unsigned int i = iter_ij.i();
        const unsigned int j = iter_ij.j();

        // orbital set for this pair = occ orbitals + PSVs
        RefSCMatrix C_occ = occ1->coefs();
        RefSCMatrix C_psv = psv_set[ij].first;
        //RefSCMatrix C_psv = r12world()->refwfn()->uocc(spin1)->coefs();
        RefSCDimension aodim = C_occ.rowdim();
        RefSCDimension modim = new SCDimension(C_occ.ncol() + C_psv.ncol(), 1);
        modim->blocks()->set_subdim(0, new SCDimension(modim.n()));
        RefSCMatrix C = C_occ.kit()->matrix(C_occ.rowdim(),
                                            modim);
        C.assign_subblock(C_occ,
                          0, aodim.n()-1,
                          0, C_occ.ncol()-1);
        C.assign_subblock(C_psv,
                          0, aodim.n()-1,
                          C_occ.ncol(), modim.n()-1);
        std::vector<unsigned int> orbsym(C.ncol(), 0); // assuming C1 symmetry!
        RefSymmSCMatrix Pa = C.kit()->symmmatrix(modim); Pa.assign(0.0);
        for(unsigned int o=0; o<nocc; ++o)
          Pa.set_element(o, o, 1.0);
        RefSymmSCMatrix Pb = Pa; // RHF only at the moment

        // reinitialize R12WavefunctionWorld
        this->obsolete();
        r12world_->world()->initialize_ao_spaces();
        Ref<RefWavefunction> refwfn = new Extern_RefWavefunction(this->r12world()->world(), this->basis(), this->integral(),
                                                                 C, orbsym, Pa, Pb, nocc,
                                                                 this->nfzcore(), this->nfzvirt(), false);
        r12world_->refwfn(refwfn);

        // compute
        this->compute_energy_();

        // extract particular pair energy
        RefSCVector er12 = r12c_energy_->ef12(S);
        RefSCVector emp2r12 = r12c_energy_->emp2f12(S);
        const double er12_ij = er12(ij);
        const double emp2_ij = emp2r12(ij) - er12_ij;
        r12_pair_energies[S].push_back(er12_ij);
        mp2_pair_energies[S].push_back(emp2_ij);

      } // end of ij loop
    } // end of spincase2 loop


    std::ostream& so = ExEnv::out0();
    so << indent << "DONE w PSV MP2-R12 !!!!!" << std::endl;
    // print out pair energies
    for(int s=0; s<2; ++s) {
      SpinCase2 S = static_cast<SpinCase2>(s);
      SpinCase1 spin1 = case1(S);
      SpinCase1 spin2 = case2(S);
      Ref<OrbitalSpace> occ1_act = r12eval()->ggspace(spin1);
      Ref<OrbitalSpace> occ2_act = r12eval()->ggspace(spin2);

      so << std::endl << indent << prepend_spincase(S,"MBPT2-R12/C") << " pair energies:" << std::endl;
      so << indent << scprintf("    i       j        mp2(ij)        f12(ij)      mp2-f12(ij)") << endl;
      so << indent << scprintf("  -----   -----   ------------   ------------   ------------") << endl;
      SpinMOPairIter ij_iter(occ1_act->rank(),occ2_act->rank(),S);
      for (ij_iter.start(); ij_iter; ij_iter.next()) {
        const int i = ij_iter.i();
        const int j = ij_iter.j();
        const int ij = ij_iter.ij();
        const double ep_f12 = r12_pair_energies[S][ij];
        const double ep_mp2 = mp2_pair_energies[S][ij];
        const double ep_mp2f12 = ep_f12 + ep_mp2;
        so
            << indent
            << scprintf("  %3d     %3d     %12.9lf   %12.9lf   %12.9lf", i + 1,
                        j + 1, ep_mp2, ep_f12, ep_mp2f12) << std::endl;
      }

    }
    // print out total energies
    double er12 = 2.0 * std::accumulate(r12_pair_energies[AlphaAlpha].begin(),
                                        r12_pair_energies[AlphaAlpha].end(),
                                        0.0);
    er12 = std::accumulate(r12_pair_energies[AlphaBeta].begin(),
                           r12_pair_energies[AlphaBeta].end(),
                           er12);
    double emp2 = 2.0 * std::accumulate(mp2_pair_energies[AlphaAlpha].begin(),
                                        mp2_pair_energies[AlphaAlpha].end(),
                                        0.0);
    emp2 = std::accumulate(mp2_pair_energies[AlphaBeta].begin(),
                           mp2_pair_energies[AlphaBeta].end(),
                           emp2);
    const double etotal = emp2 + er12 + this->ref_energy();
    so <<endl<<indent
       <<scprintf("RHF energy [au]:                               %17.12lf\n", this->ref_energy());
    so <<indent
       <<scprintf("MP2 correlation energy [au]:                   %17.12lf\n", emp2);
    so <<indent
       <<scprintf("(MBPT2)-R12/C   correlation energy [au]:       %17.12lf\n", er12);
    so <<indent
       <<scprintf("MBPT2-R12/C   correlation energy [au]:         %17.12lf\n", emp2+er12);
    so <<indent
       <<scprintf("MBPT2-R12/C   energy [au]:                     %17.12lf\n", etotal) << endl;

  }
#endif

}

//////////////////////////////////////////////////////////////////////////////

void
MBPT2_R12::obsolete()
{
  r12eval_ = 0;
  r12a_energy_ = 0;
  r12ap_energy_ = 0;
  r12app_energy_ = 0;
  r12b_energy_ = 0;
  r12c_energy_ = 0;
  mp2_corr_energy_ = 0.0;
  cabs_singles_energy_ = 0.0;
  r12world_->world()->obsolete();
  r12world_->obsolete();
  MBPT2::obsolete();
}

//////////////////////////////////////////////////////////////////////////////

bool
MBPT2_R12::analytic_gradient_implemented() const
{
  return false;
}

//////////////////////////////////////////////////////////////////////////////

int
MBPT2_R12::value_implemented() const
{
  return 1;
}

////////////////////////////////////////////////////////////////////////////

void
MBPT2_R12::corrfactor(const Ref<R12Technology::CorrelationFactor>& cf)
{
    if (!r12world()->r12tech()->corrfactor()->equiv(cf)) {
      r12world()->r12tech()->corrfactor(cf);
      obsolete();
  }
}

/////////////////////////////////////////////////////////////////////////////

double
MBPT2_R12::corr_energy()
{
  energy();
  return energy() - ref_energy() - cabs_singles_energy();
}

/////////////////////////////////////////////////////////////////////////////

double
MBPT2_R12::r12_corr_energy()
{
  energy();
  return energy() - mp2_corr_energy_ - ref_energy() - cabs_singles_energy();
}

/////////////////////////////////////////////////////////////////////////////

double
MBPT2_R12::cabs_singles_energy()
{
  energy();
  return cabs_singles_energy_;
}

/////////////////////////////////////////////////////////////////////////////

namespace {

  int triang_half_INDEX_ordered(int i, int j) {
    return(i*(i+1)/2+j);
  }

  int triang_half_INDEX(int i, int j) {
    return((i>j) ? triang_half_INDEX_ordered(i,j) : triang_half_INDEX_ordered(j,i));
  }

  int ordinary_INDEX(int i, int j, int coldim) {
    return(i*coldim+j);
  }

  /// tpdm_index and init_ioff: for indexing of density matrices.
  int tpdm_index(int i, int j, int k, int l,int coldim){
    //int ind_half1=triang_half_INDEX(i,j);
    //int ind_half2=triang_half_INDEX(k,l);
    int ind_half1=ordinary_INDEX(i,j,coldim);
    int ind_half2=ordinary_INDEX(k,l,coldim);
    return(triang_half_INDEX(ind_half1,ind_half2));
  }

  void vector_to_symmmatrix(RefSymmSCMatrix &matrix, const RefSCVector &vector) {
    int dim = matrix.dim().n();
    for(int i=0; i<dim; i++){
      for(int j=0; j<=i; j++) {
        matrix.set_element(i,j,vector.get_element(triang_half_INDEX(i,j)));
      }
    }
  }

  void symmmatrix_to_vector(RefSCVector &vector, const RefSymmSCMatrix &matrix) {
    int dim = matrix.dim().n();
    for(int i=0; i<dim; i++){
      for(int j=0; j<=i; j++) {
        vector.set_element(triang_half_INDEX(i,j),matrix.get_element(i,j));
      }
    }
  }

  void vector_to_matrix(RefSCMatrix &matrix, const RefSCVector &vector) {
    int dim1 = matrix.rowdim().n();
    int dim2 = matrix.coldim().n();
    for(int i=0; i<dim1; i++) {
      for(int j=0; j<dim2; j++) {
        matrix.set_element(i,j,vector.get_element(ordinary_INDEX(i,j,dim2)));
      }
    }
  }

  int lowerupper_index(int p, int q);

  void vector_to_matrix(RefSCMatrix &matrix,const RefSCVector &vector,const SpinCase2 &pairspin) {
    int dim1 = matrix.rowdim().n();
    int dim2 = matrix.coldim().n();
    if(pairspin==AlphaBeta) {
      for(int i=0; i<dim1; i++) {
        for(int j=0; j<dim2; j++) {
          matrix.set_element(i,j,vector.get_element(ordinary_INDEX(i,j,dim2)));
        }
      }
    }
    else {  // pairspin==AlphaAlpha || pairspin==BetaBeta
      matrix->assign(0.0);
      for(int i=0; i<dim1; i++) {
        for(int j=0; j<i; j++) {
          const double value = vector.get_element(lowerupper_index(i,j));
          matrix.set_element(i,j,value);
          matrix.set_element(j,i,-value);
        }
      }
    }
  }

  void matrix_to_vector(RefSCVector &vector, const RefSCMatrix &matrix) {
    int dim1 = matrix.rowdim().n();
    int dim2 = matrix.coldim().n();
    for(int i=0; i<dim1; i++) {
      for(int j=0; j<dim2; j++) {
        vector.set_element(ordinary_INDEX(i,j,dim2),matrix.get_element(i,j));
      }
    }
  }

  void matrix_to_vector(RefSCVector &vector, const RefSCMatrix &matrix,const SpinCase2 &pairspin) {
    int dim1 = matrix.rowdim().n();
    int dim2 = matrix.coldim().n();
    if(pairspin==AlphaBeta) {
      for(int i=0; i<dim1; i++) {
        for(int j=0; j<dim2; j++) {
          vector.set_element(ordinary_INDEX(i,j,dim2),matrix.get_element(i,j));
        }
      }
    }
    else {  // pairspin==AlphaAlpha || pairspin==BetaBeta
      for(int i=0; i<dim1; i++) {
        for(int j=0; j<i; j++) {
          vector.set_element(lowerupper_index(i,j),matrix.get_element(i,j));
        }
      }
    }
  }

  int lowertriang_index(int p,int q) {
    if(q>=p){
      throw ProgrammingError("lowertriang_index(p,q) -- q must be smaller than p.",__FILE__,__LINE__);
    }
    int index=p*(p+1)/2+q-p;
    return(index);
  }

  int lowerupper_index(int p, int q) {
    if(p>q) {
      return(lowertriang_index(p,q));
    }
    else if(q>p) {
      return(lowertriang_index(q,p));
    }
    else {
      throw ProgrammingError("lowerupper_index(p,q) -- p and q are not allowed to be equal.",__FILE__,__LINE__);
    }
  }
  double indexsizeorder_sign(int p,int q) {
    if(p>q) {
      return(1.0);
    }
    else if(q>p) {
      return(-1.0);
    }
    else {
      return(0.0);
    }
  }

  int antisym_pairindex(int i, int j) {
    int max_ij = std::max(i, j);
    int min_ij = std::min(i, j);
    return (max_ij -1)* max_ij/2 + min_ij;
  }

}

std::vector< std::pair<RefSCMatrix,RefSCMatrix> >
MBPT2_R12::mp1_pno(SpinCase2 spin,
                   bool deflate_geminal,
                   double truncate_threshold)
{
  typedef std::vector< std::pair<RefSCMatrix,RefSCMatrix> > result_type;
  result_type result;

  ExEnv::out0() << indent << "Computing PNOs (spin=" << to_string(spin) << ",r12=" << (deflate_geminal?"yes":"no") << ")"
                          << endl;

  Ref<R12Amplitudes> amps = r12eval()->amps();
  {
    const SpinCase2 spincase2 = static_cast<SpinCase2>(spin);
    const SpinCase1 spin1 = case1(spincase2);
    const SpinCase1 spin2 = case2(spincase2);

    const Ref<OrbitalSpace>& occ1_act = r12eval()->occ_act(spin1);
    const Ref<OrbitalSpace>& occ2_act = r12eval()->occ_act(spin2);
    const Ref<OrbitalSpace>& occ1 = r12eval()->occ(spin1);
    const Ref<OrbitalSpace>& occ2 = r12eval()->occ(spin2);
    const Ref<OrbitalSpace>& cabs1 = r12world()->cabs_space(spin1);
    const Ref<OrbitalSpace>& cabs2 = r12world()->cabs_space(spin2);
    const Ref<OrbitalSpace>& uocc1_act = r12world()->refwfn()->uocc_act(spin1);
    const Ref<OrbitalSpace>& uocc2_act = r12world()->refwfn()->uocc_act(spin2);
    MPQC_ASSERT(uocc1_act == uocc2_act); // not ready to handle UHF yet

    // compute the vv block of MP1 (not MP2-R12) 1-RDM for each ij
    // \gamma_a^b = T^ij_ac T_ij^bc

    SpinMOPairIter ij_iter( occ1_act->rank(), occ2_act->rank(), spincase2);
    SpinMOPairIter ab_iter(uocc1_act->rank(),uocc2_act->rank(), spincase2);

    RefSCMatrix T2 = amps->T2(spincase2);
    //T2.print("MP1 T2 amplitudes");

    RefSCMatrix Fvv = amps->Fvv(spincase2);
    RefSCMatrix Fox = amps->Fox(spincase2);
    RefSCMatrix Fxo = amps->Fxo(spincase2);
    MPQC_ASSERT(this->r12world()->r12tech()->ansatz()->diag() == true);  // assuming diagonal ansatz
    RefSCMatrix Tg = Fvv.kit()->matrix(Fvv.rowdim(),Fvv.rowdim()); Tg.assign(0.0);
    firstorder_cusp_coefficients(spincase2, Tg, occ1_act, occ2_act, this->r12world()->r12tech()->corrfactor());
    Fvv = Tg * Fvv;   // to obtain amplitudes we need to contract in cusp coefficients
    Fox = Tg * Fox;   // to obtain amplitudes we need to contract in cusp coefficients
    Fxo = Tg * Fxo;   // to obtain amplitudes we need to contract in cusp coefficients
    //Fvv.print("R12 geminal amplitudes in vv space");

    // compute the AO basis Fock matrix so that PSV sets can be canonicalized
    RefSymmSCMatrix F_ao;
    {
      Ref<OrbitalSpace> obs_ao_space = this->r12world()->world()->tfactory()->orbital_registry()->value("mu");
      RefSCMatrix F_ao_rect = r12eval()->fock(obs_ao_space, obs_ao_space, AnySpinCase1);
      F_ao = F_ao_rect.kit()->symmmatrix(F_ao_rect.rowdim()); F_ao.assign(0.0);
      F_ao.accumulate_symmetric_sum(F_ao_rect); F_ao.scale(0.5);
    }

    for(ij_iter.start(); ij_iter; ij_iter.next()) {
      const int ij = ij_iter.ij();
      const unsigned int i = ij_iter.i();
      const unsigned int j = ij_iter.j();
      std::ostringstream oss; oss << "act occ pair " << i << "," << j;

      std::pair<RefSCMatrix, RefSCMatrix> result_ij; // pno,cabs coefficient pair

      RefSCVector T2_ij_vec = T2.get_row(ij);
      RefSCMatrix T2_ij = T2_ij_vec.kit()->matrix(uocc1_act->dim(), uocc2_act->dim());
      vector_to_matrix(T2_ij, T2_ij_vec, spincase2);
      //T2_ij.print((std::string("MP1 T2 amplitudes ") + oss.str()).c_str());

      RefSCMatrix T2_ij_eff = T2_ij;

      // compute the pair density and eigendecompose it
      RefSymmSCMatrix gamma1ab_ij = T2_ij.kit()->symmmatrix(T2_ij.rowdim()); gamma1ab_ij.assign(0.0);
      RefDiagSCMatrix gamma1ab_ij_evals;
      RefSCMatrix gamma1ab_ij_evecs;
      RefDiagSCMatrix gamma1_xx_evals;
      RefSCMatrix gamma1_xx_evecs;
      if (deflate_geminal == false) { // use standard MP1 pair densities
        gamma1ab_ij.accumulate_symmetric_product(T2_ij);
        gamma1ab_ij.accumulate_symmetric_product(T2_ij.t());
        //gamma1ab_ij.print("pair density");

//        RefSCMatrix gamma1ab_ij_rect = T2_ij * T2_ij.t() + T2_ij.t() * T2_ij;
//        gamma1ab_ij.assign(0.0);  gamma1ab_ij.accumulate_symmetric_sum(gamma1ab_ij_rect); gamma1ab_ij.scale(0.5);
//        gamma1ab_ij_rect.print("pair density?");
//        gamma1ab_ij.print("pair density???");

        gamma1ab_ij_evals = gamma1ab_ij.eigvals();
        gamma1ab_ij_evals.print("pair density eigenvalues");
      }
      else { // use densities deflated by the geminal contributions
        RefSCVector Fvv_ij_vec = Fvv.get_row(ij);
        RefSCMatrix Fvv_ij = Fvv_ij_vec.kit()->matrix(uocc1_act->dim(), uocc2_act->dim());
        vector_to_matrix(Fvv_ij, Fvv_ij_vec, spincase2);
        //Fvv_ij.print((std::string("R12 geminal amplitudes in vv space ") + oss.str()).c_str());

        // deflate T2 amplitudes by the geminal amplitudes in vv block
        RefSymmSCMatrix gamma1ab_ij_sansr12 = T2_ij.kit()->symmmatrix(T2_ij.rowdim()); gamma1ab_ij_sansr12.assign(0.0);
        T2_ij_eff = T2_ij - Fvv_ij;
#define DEFLATE_AMPLITUDES 1
#if DEFLATE_AMPLITUDES
        gamma1ab_ij_sansr12.accumulate_symmetric_product(T2_ij_eff);
        gamma1ab_ij_sansr12.accumulate_symmetric_product(T2_ij_eff.t());
        //gamma1ab_ij_sansr12.eigvals().print("pair density (sans R12 vv) eigenvalues");
#else
        gamma1ab_ij_sansr12.accumulate_symmetric_product(Fvv_ij);
        gamma1ab_ij_sansr12.accumulate_symmetric_product(Fvv_ij.t());
        gamma1ab_ij_sansr12.scale(-1.0);
        gamma1ab_ij_sansr12.accumulate_symmetric_product(T2_ij);
        gamma1ab_ij_sansr12.accumulate_symmetric_product(T2_ij.t());
        gamma1ab_ij_sansr12.eigvals().print("pair density (sans R12 vv) eigenvalues");
#endif
        // pass it outside this scope
        gamma1ab_ij = gamma1ab_ij_sansr12;

        gamma1ab_ij_evals = gamma1ab_ij_sansr12.eigvals();
        gamma1ab_ij_evals.print("pair density (sans R12) eigenvalues");

        // compute the "natural" CABS space
        if (1) {
          RefSCVector Fox_ij_vec = Fox.get_row(ij);
          RefSCMatrix Fox_ij = Fox_ij_vec.kit()->matrix(occ1->dim(),
                                                        cabs2->dim());
          vector_to_matrix(Fox_ij, Fox_ij_vec, AlphaBeta);
          RefSCVector Fxo_ij_vec = Fxo.get_row(ij);
          RefSCMatrix Fxo_ij = Fxo_ij_vec.kit()->matrix(cabs1->dim(),
                                                        occ2->dim());
          vector_to_matrix(Fxo_ij, Fxo_ij_vec, AlphaBeta);

          RefSymmSCMatrix gamma1_xx = gamma1ab_ij.kit()->symmmatrix(cabs1->dim());
          gamma1_xx.assign(0.0);
          gamma1_xx.accumulate_symmetric_product(Fxo_ij);
          gamma1_xx.accumulate_symmetric_product(Fox_ij.t());
          gamma1_xx_evals = gamma1_xx.eigvals();
          gamma1_xx_evals.print("eigenvalues of cabs/cabs pair density");
          gamma1_xx_evecs = gamma1_xx.eigvecs();
        }

      }
      gamma1ab_ij_evecs = gamma1ab_ij.eigvecs();

      // coefficient matrices should use blocked dimensions for using with OrbitalSpace
      {
        RefSCMatrix gamma1ab_ij_evecs_blk = uocc1_act->coefs()->kit()->matrix(uocc1_act->dim(),
                                                                              uocc1_act->dim());
        gamma1ab_ij_evecs_blk.assign(0.0);
        gamma1ab_ij_evecs_blk->convert_accumulate(gamma1ab_ij_evecs);
        gamma1ab_ij_evecs = gamma1ab_ij_evecs_blk;
      }
      {
        RefSCMatrix gamma1_xx_evecs_blk = uocc1_act->coefs()->kit()->matrix(cabs1->dim(),
                                                                            cabs1->dim());
        gamma1_xx_evecs_blk.assign(0.0);
        gamma1_xx_evecs_blk->convert_accumulate(gamma1_xx_evecs);
        gamma1_xx_evecs = gamma1_xx_evecs_blk;
      }

      // SVD the amplitudes
      if (0) {
        RefDiagSCMatrix T2_ij_svals = T2_ij.kit()->diagmatrix(T2_ij.rowdim());
        RefSCMatrix T2_ij_svecs_U = T2_ij.clone();
        RefSCMatrix T2_ij_svecs_V = T2_ij.clone();
        T2_ij_eff.svd(T2_ij_svecs_U, T2_ij_svals, T2_ij_svecs_V);
        T2_ij_svals.print("T2(eff) singular values");
//        T2_ij_svecs_U.print("T2(eff) singular vectors U");
//        T2_ij_svecs_V.print("T2(eff) singular vectors V");
        RefSCMatrix T2_ij_svecs_U_blk = uocc1_act->coefs()->kit()->matrix(T2_ij_svecs_U.rowdim(),
                                                                          T2_ij_svecs_U.coldim());
        T2_ij_svecs_U_blk.assign(0.0);
        T2_ij_svecs_U_blk->convert_accumulate(T2_ij_svecs_U);
        RefSCMatrix T2_ij_svecs_V_blk = uocc1_act->coefs()->kit()->matrix(T2_ij_svecs_V.rowdim(),
                                                                          T2_ij_svecs_V.coldim());
        T2_ij_svecs_V_blk.assign(0.0);
        T2_ij_svecs_V_blk->convert_accumulate(T2_ij_svecs_V);
//        (uocc1_act->coefs() * T2_ij_svecs_U_blk).print("T2(eff) singular vectors U in AO basis");
//        (uocc1_act->coefs() * T2_ij_svecs_V_blk).print("T2(eff) singular vectors V in AO basis");
      }

//      // canonicalize the PSV set before truncation, purely to be able to compare with the "after" scenario
//      {
//        RefSCMatrix psv_coefs_ao = uocc1_act->coefs() * gamma1ab_ij_evecs;
//        RefSymmSCMatrix F_pno = F_ao.kit()->symmmatrix(pno_coefs_ao.coldim()); F_pno.assign(0.0);
//        F_pno.accumulate_transform(pno_coefs_ao, F_ao, SCMatrix::TransposeTransform);
//        RefSCMatrix U = F_pno.eigvecs();
//        F_pno.eigvals().print("PNO orbital energies (before truncation)");
//      }
//
      // truncate the PNO sets according to their eigenvalues ("occupation numbers")
      RefSCMatrix pno_coefs_ao;
      {
        const unsigned int nuocc_act = uocc1_act->rank();
        long double sum_of_omitted_evals = 0.0;
        unsigned int a;
        for (a = 0; a < nuocc_act; ++a) {
          const double eval = std::fabs((long double) gamma1ab_ij_evals(a));
          if (eval + DBL_EPSILON < truncate_threshold)
            sum_of_omitted_evals += eval;
          else
            break;
        }
        const unsigned int nomitted = a;
        ExEnv::out0() << oss.str() << ": omitted " << nomitted
            << " eigenpairs (sum of omitted evals = " << sum_of_omitted_evals
            << ")" << std::endl;
        RefSCDimension pno_screened_dim = new SCDimension(nuocc_act - nomitted);
        pno_screened_dim->blocks()->set_subdim(
            0, new SCDimension(pno_screened_dim.n()));
        RefSCMatrix gamma1ab_ij_evecs_screened =
            gamma1ab_ij_evecs.kit()->matrix(gamma1ab_ij_evecs.rowdim(),
                                            pno_screened_dim);
        gamma1ab_ij_evecs_screened.assign_subblock(gamma1ab_ij_evecs,
                                                   0, gamma1ab_ij_evecs.nrow() - 1,
                                                   0, pno_screened_dim.n() - 1,
                                                   0, nomitted);
        pno_coefs_ao = uocc1_act->coefs() * gamma1ab_ij_evecs_screened;
      }

      // truncate CABS
      RefSCMatrix cabs_coefs_ao;
      {
        const unsigned int ncabs = cabs1->rank();
        long double sum_of_omitted_evals = 0.0;
        unsigned int a;
        for (a = 0; a < ncabs; ++a) {
          const double eval = std::fabs((long double) gamma1_xx_evals(a));
          if (eval + DBL_EPSILON < 1e-12)
            sum_of_omitted_evals += eval;
          else
            break;
        }
        const unsigned int nomitted = a;
        ExEnv::out0() << oss.str() << ": omitted " << nomitted
            << " CABS vectors (sum of omitted evals = " << sum_of_omitted_evals
            << ")" << std::endl;
        if (ncabs > nomitted) {
          RefSCDimension cabs_screened_dim = new SCDimension(ncabs - nomitted);
          cabs_screened_dim->blocks()->set_subdim(
              0, new SCDimension(cabs_screened_dim.n()));
          RefSCMatrix gamma1_xx_evecs_screened = gamma1_xx_evecs.kit()->matrix(
              gamma1_xx_evecs.rowdim(), cabs_screened_dim);
          gamma1_xx_evecs_screened.assign_subblock(gamma1_xx_evecs, 0,
                                                   gamma1_xx_evecs.nrow() - 1,
                                                   0, cabs_screened_dim.n() - 1,
                                                   0, nomitted);
          cabs_coefs_ao = cabs1->coefs() * gamma1_xx_evecs_screened;
        }
      }
      result_ij.second = cabs_coefs_ao;

      // canonicalize the resulting PNOs
      {
        RefSymmSCMatrix F_pno = F_ao.kit()->symmmatrix(pno_coefs_ao.coldim()); F_pno.assign(0.0);
        F_pno.accumulate_transform(pno_coefs_ao, F_ao, SCMatrix::TransposeTransform);
        RefSCMatrix U = F_pno.eigvecs();
        F_pno.eigvals().print("PNO orbital energies (after truncation)");
        RefSCMatrix pno_coefs_canonical_ao = pno_coefs_ao * U;

        // and finally, save the result
        result_ij.first = pno_coefs_canonical_ao;
      }

      result.push_back(result_ij);
    } // loop over act occ pairs

  } // spincase block

  return result;
}

//////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "ETS"
// End:
