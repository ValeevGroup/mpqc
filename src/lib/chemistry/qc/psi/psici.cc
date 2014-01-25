//
// psici.cc
//
// Copyright (C) 2008 Martin Torheyden
//
// Author: Martin Torheyden <mtorhey@vt.edu>
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


#include <cassert>
#include <psifiles.h>
#include <ccfiles.h>
#include <cmath>
#include <numeric>

#include <math/mmisc/pairiter.impl.h>
#include <util/misc/print.h>
#include <chemistry/qc/wfn/orbitalspace_utils.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/psi/psiwfn.h>
#include <chemistry/qc/psi/psici.h>
#include <chemistry/qc/psi/psiqtorder.h>

#include <psifiles.h>
#include <libiwl/iwl.h>
#include <libciomr/libciomr.h>

using namespace std;

namespace sc {

  namespace detail {

    vector<int> get_block_dims(const RefSCDimension &dim) {
      Ref<SCBlockInfo> blocks = dim->blocks();
      int nblocks = blocks->nblock();
      vector<int> block_dims(nblocks);
      for(int i=0; i<nblocks; i++) {
        block_dims[i] = blocks->subdim(i);
      }

      return(block_dims);
    }

    vector<int> get_block_dim_offsets(const RefSCDimension &dim) {
      vector<int> block_dims = get_block_dims(dim);
      const int ndim = block_dims.size();
      vector<int> block_dim_offsets(ndim);

      block_dim_offsets[0] = 0;
      //block_dim_offsets[1] = block_dims[0];
      for(int i=1; i<ndim; i++) {
        block_dim_offsets[i] = block_dim_offsets[i-1] + block_dims[i-1];
      }
      return(block_dim_offsets);
    }


    void print_block_dims(const RefSCDimension &dim) {
      vector<int> dims = get_block_dims(dim);
      ExEnv::out0() << indent;
      for(int i=0; i<dims.size(); i++) {
        ExEnv::out0() << setw(5) << dims[i];
      }
      ExEnv::out0() << endl;

      vector<int> dim_offsets = get_block_dim_offsets(dim);
      ExEnv::out0() << indent;
      for(int i=0; i<dim_offsets.size(); i++) {
        ExEnv::out0() << setw(5) << dim_offsets[i];
      }
      ExEnv::out0() << endl;
    }

    template <typename Int> void print_vector_int(const vector<Int> &vec_int) {
      ExEnv::out0() << indent;
      for(int i=0; i<vec_int.size(); i++) {
        ExEnv::out0() << setw(5) << vec_int[i];
      }
      ExEnv::out0() << endl;
    }

    void print_vector_vector_int(const vector<vector<int> > &vec_vec_int) {
      int nrows = vec_vec_int.size();
      for(int i=0; i<nrows; i++) {
        print_vector_int(vec_vec_int[i]);
      }
    }

    void print_vector_vector_vector_int(const vector<vector<vector<int> > > &vec_vec_vec_int) {
      int nblks = vec_vec_vec_int.size();
      ExEnv::out0() << indent;
      for(int i=0; i<nblks; i++) {
        ExEnv::out0() << "block " << (i+1) << ":" << endl;
        print_vector_vector_int(vec_vec_vec_int[i]);
      }
    }

    vector<double> getrow(int index,const RefSCMatrix &Matrix) {
      int ndim = Matrix.coldim().n();
      vector<double> rowvec(ndim);
      for(int i=0; i<ndim; i++) {
        rowvec[i] = Matrix.get_element(index,i);
      }

      return(rowvec);
    }

    /** a matrix element
     */
    class Element {
      public:
        int ind1;        // MO index in reference wave function
        int ind2;        // MO index in wave function
        double value;    // overlap of the ind1 MO -- ind2 MO overlap
        ~Element() {}
        Element(int Ind1,int Ind2,double Value) : ind1(Ind1), ind2(Ind2), value(Value) {}
        Element(const Element &other) {
          ind1=other.ind1; ind2=other.ind2; value=other.value;
        }
        Element &operator=(const Element &other) {
          if(this==&other) return(*this);
          ind1=other.ind1; ind2=other.ind2; value=other.value;
          return(*this);
        }
        bool operator>(const Element &other) const {return(fabs(value)>fabs(other.value));}
        bool operator<(const Element &other) const {return(fabs(value)<fabs(other.value));}
        void print() const {ExEnv::out0() << setw(7) << ind1 << setw(7) << ind2
                                          << setw(18) << setprecision(10) << value << endl; }
    };

    /** Class for symmetrry indices with energies for MOs. */
    class Energy_SymmIndex {
      public:
        int symmindex;  /// index number of symmetry block of MO symmetry
        double energy;  /// MO energy
        ~Energy_SymmIndex(){}
        Energy_SymmIndex(){}
        Energy_SymmIndex(int Symmindex,double Energy) : symmindex(Symmindex), energy(Energy) {}
        Energy_SymmIndex(const Energy_SymmIndex &other) {
          symmindex = other.symmindex; energy = other.energy;
        }
        Energy_SymmIndex &operator=(const Energy_SymmIndex &other) {
          if(this==&other) return(*this);
          symmindex = other.symmindex; energy = other.energy;
          return(*this);
        }
        bool operator>(const Energy_SymmIndex &other) const {return(energy>other.energy);}
        bool operator<(const Energy_SymmIndex &other) const {return(energy<other.energy);}
    };

    /** Given the MO overlap matrix between valence_obwfn_ (row index) and
     *  the wave function, create a list of objects of the class Element.
     */
    list<Element> matrix_to_element_list(const RefSCMatrix &value_mat) {
      RefSCDimension coldim = value_mat->coldim();
      RefSCDimension rowdim = value_mat->rowdim();
      int ncol_blks = coldim->blocks()->nblock();
      int nrow_blks = rowdim->blocks()->nblock();

      if(ncol_blks != nrow_blks)
        throw ProgrammingError("Error in detail::get_index_mat_from_value_mat -- number of block in rowdim and coldim does not coincide.",__FILE__,__LINE__);
      int nblks = nrow_blks;

      int ncols = coldim.n();
      int nrows = rowdim.n();
      //Ref<SCBlockInfo> colblocks = coldim->blocks();
      //Ref<SCBlockInfo> rowblocks = rowdim->blocks();
      vector<int> ncol_blk_dim = get_block_dims(coldim);
      vector<int> nrow_blk_dim = get_block_dims(rowdim);
      vector<int> ncol_blk_dim_offsets = get_block_dim_offsets(coldim);
      vector<int> nrow_blk_dim_offsets = get_block_dim_offsets(rowdim);

      list<Element> elem_list;

      for(int i=0; i<nblks; i++) {
        for(int j=0; j<nrow_blk_dim[i]; j++) {
          int rowind = nrow_blk_dim_offsets[i] + j;
          for(int k=0; k<ncol_blk_dim[i]; k++) {
            int colind = ncol_blk_dim_offsets[i] + k;
            double value = value_mat.get_element(rowind,colind);
            Element elem(rowind,colind,value);
            //elem.print();
            elem_list.push_back(elem);
          }
        }
      }

      return(elem_list);
    }

    void print_elem_list(const list<Element> &Elem_list) {
      for(list<Element>::const_iterator it=Elem_list.begin(); it!=Elem_list.end(); ++it) {
        it->print();
      }
    }

    /** Returns a map between the MOs of the reference wave function (first index) and
     *  the MOs of the wave function between which the overlap is the largest. The Element
     *  list has to be sorted.
     */
    std::map<unsigned int, unsigned int> reference_val_map(const list<detail::Element> &sorted_elem_list) {
      // the dimension of the reference_wfn and the wave function is the maximum of
      // element->ind1 and element->ind respectively plus one for all elements element
      // of sorted_elem_list.
      int dim_ref=0;
      int dim_wfn=0;
      for(list<Element>::const_iterator elem=sorted_elem_list.begin(); elem!=sorted_elem_list.end(); ++elem) {
        if(elem->ind1 > dim_ref) {
          dim_ref = elem->ind1;
        }
        if(elem->ind2 > dim_wfn) {
          dim_wfn = elem->ind2;
        }
      }
      dim_ref += 1;
      dim_wfn +=1;

      std::map<unsigned int, unsigned int> index_map;
      vector<int> chk_ref(dim_ref);
      vector<int> chk_wfn(dim_wfn);
      chk_ref.assign(dim_ref,0);
      chk_wfn.assign(dim_wfn,0);
      /// loop starting with element of largest overlap
      for(list<Element>::const_iterator elem=sorted_elem_list.begin(); elem!=sorted_elem_list.end(); ++elem) {
        int ind_ref = elem->ind1;
        int ind_wfn = elem->ind2;
        if((chk_ref[ind_ref]==0) && (chk_wfn[ind_wfn]==0)) {  // each index of each dimension can be selected only once.
          chk_ref[ind_ref] = 1; chk_wfn[ind_wfn] = 1;
          index_map.insert(pair<unsigned int, unsigned int>(ind_ref,ind_wfn));  // insert pair of indices whose overlap is the largest one.
        }
      }
      return(index_map);
    }

    template <typename Int>
    void print_int_index_map(const std::map<Int,Int> &Map) {
      for(typename std::map<Int,Int>::const_iterator elem = Map.begin(); elem != Map.end(); ++elem) {
        ExEnv::out0() << setw(7) << elem->first << setw(7) << elem->second << endl;
      }
    }

    template <typename Int> void print_blocks(const char *orbs_name,const vector<Int> &orbs_blocks,ostream &o) {
      o << indent << setw(12) << orbs_name << " = [";
      for(int i=0; i<orbs_blocks.size(); i++) {
        o << setw(5) << orbs_blocks[i];
      }
      o << "]" << endl;
    }

    RefSCMatrix obwfn_overlap(const Ref<PsiSCF>& psiwfn,
                              const Ref<OneBodyWavefunction>& mpqcwfn) {
      RefSCMatrix eigenvec_wfn_ao = psiwfn->coefs(); // these are in AO basis!
      RefSCMatrix eigenvec_reference_wfn = mpqcwfn->eigenvectors(); // these are in SO basis!
      Ref<GaussianBasisSet> basis_wfn = psiwfn->basis();
      Ref<GaussianBasisSet> basis_reference_wfn = mpqcwfn->basis();
      Ref<Integral> integral_wfn = psiwfn->integral();
      Ref<Integral> integral_reference_wfn = mpqcwfn->integral();
      Ref<PetiteList> plist_wfn = integral_wfn->petite_list();
      Ref<PetiteList> plist_reference_wfn = integral_reference_wfn->petite_list();
      RefSCMatrix eigenvec_valence_obwfn_ao = plist_reference_wfn->evecs_to_AO_basis(eigenvec_reference_wfn);

      Ref<OrbitalSpace> space_wfn = new OrbitalSpace("WFN","SCF MO space",eigenvec_wfn_ao,basis_wfn,integral_wfn);
      Ref<OrbitalSpace> space_reference_wfn = new OrbitalSpace("REFERENCE","SCF reference MO space",eigenvec_valence_obwfn_ao,basis_reference_wfn,integral_reference_wfn);

      RefSCMatrix overlap = compute_overlap_ints(space_reference_wfn,space_wfn);

      return overlap;
    }

    /** Returns a map between the MOs of the reference wave function (first index) and
     *  the MOs of the wave function between which the overlap is the largest. */
    std::map<unsigned int, unsigned int> reference_index_map(const RefSCMatrix &overlap) {
      /// create list of overlap matrix elements of type 'Element'
      list<detail::Element> elem_list = detail::matrix_to_element_list(overlap);

      /// sort the list of overlap matrix elements
      elem_list.sort(greater<detail::Element>());

      std::map<unsigned int,unsigned int> index_map = detail::reference_val_map(elem_list);

      return(index_map);
    }

    /** Returns a map between the valence_obwfn_->nmo MOs and the lowest MO's of the
     *  corresponding symmetries of the wave function in symmetry blocked MO order. */
    std::map<unsigned int,unsigned int> standard_index_map(const RefSCMatrix &overlap) {
      list<detail::Element> elem_list = detail::matrix_to_element_list(overlap);
      std::map<unsigned int, unsigned int> index_map = detail::reference_val_map(elem_list);

      return(index_map);
    }

    /** Index map from standard symmetry blocked order (as it is used in Psi3) to
     *  a symmetry blocked order where the lowest orbitals correspond to the valence orbitals
     *  from valence_obwfn_. */
    std::map<unsigned int, unsigned int>
    psi_index_map(const RefSCMatrix &overlap) {
      std::map<unsigned int, unsigned int> std_index_map = standard_index_map(overlap);
      std::map<unsigned int, unsigned int> ref_index_map = reference_index_map(overlap);
      std::map<unsigned int, unsigned int> new_index_map;
      for(std::map<unsigned int, unsigned int>::const_iterator elem1 = std_index_map.begin(); elem1 != std_index_map.end(); ++elem1) {
        for(std::map<unsigned int, unsigned int>::const_iterator elem2 = ref_index_map.begin(); elem2 != ref_index_map.end(); ++elem2) {
          if((elem1->first==elem2->first) && (elem1->second!=elem2->second)) {
            new_index_map.insert(pair<int,int>(elem1->second,elem2->second));
          }
        }
      }

      return(new_index_map);
    }


    /** Returns a vector of orbital indices in symmetry blocked order where the lowest
     *  MOs correspond to the respective MOs in the valence wfn. This vector
     *  can be used as the value of the "reorder" keyword in Psi to make the lowest-energy orbitals
     *  maximally similar to the valence orbitals. This is needed, for example, to avoid the inclusion
     *  of diffuse orbitals into the active space of RASSCF calculations.
     *
     *  \param overlap is the overlap between a valence wfn (e.g. minimal basis HF) and a given wfn
     */
    vector<unsigned int>
    ref_to_valence_moorder(const Ref<PsiSCF>& psiwfn,
                           const Ref<OneBodyWavefunction>& mpqcwfn) {
      const int nmo = psiwfn->nmo();
      RefSCMatrix overlap = obwfn_overlap(psiwfn, mpqcwfn);
      std::map<unsigned int, unsigned int> index_map = psi_index_map(overlap);
      vector<unsigned int> orbital_order;
      if(!index_map.empty()) {
        orbital_order.resize(nmo);
        for(int i=0; i<nmo; i++) orbital_order[i] = i;
        for(std::map<unsigned int, unsigned int>::const_iterator elem = index_map.begin(); elem != index_map.end(); ++elem) {
          const unsigned int ind1 = elem->first;
          const unsigned int ind2 = elem->second;
          const unsigned int tmp = orbital_order[ind2];
          orbital_order[ind2] = orbital_order[ind1];
          orbital_order[ind1] = tmp;
        }
      }

      return orbital_order;
    }

  } // end of namespace detail

  static ClassDesc PsiRASCI_cd(typeid(PsiRASCI), "PsiRASCI", 1,
                                   "public PsiCorrWavefunction",  0, create<PsiRASCI>,
                                   create<PsiRASCI>);

  PsiRASCI::PsiRASCI(const Ref<KeyVal> &keyval)
    : PsiCorrWavefunction(keyval) {

    opdm_print_ = keyval->booleanvalue("opdm_print", KeyValValueboolean(false));
    tpdm_print_ = keyval->booleanvalue("tpdm_print", KeyValValueboolean(false));

    root_ = keyval->intvalue("root", KeyValValueint(1));
    multiplicity_ = keyval->intvalue("multiplicity", KeyValValueint( this->reference()->multiplicity() ));
    nroots_ = keyval->intvalue("nroots", KeyValValueint(root_));
    if(nroots_<root_) {
      throw InputError("PsiRASCI::PsiRASCI(const Ref<KeyVal> &) -- nroots should not be smaller than root.",__FILE__,__LINE__);
    }
    // the default = -1 is to leave it out
    target_sym_ = keyval->intvalue("target_sym", KeyValValueint(-1));
    if (keyval->exists("target_sym") && (target_sym_ < 0 || target_sym_ >= this->molecule()->point_group()->char_table().nirrep()))
      throw InputError("PsiRASCI::PsiRASCI(const Ref<KeyVal> &) -- target_sym must be between 0 and number of irreps - 1",__FILE__,__LINE__);

    h0_blocksize_ = keyval->intvalue("h0_blocksize", KeyValValueint(40));
    repl_otf_ = keyval->booleanvalue("repl_otf", KeyValValueboolean(false));

    /// construction of valence_obwfn_
    valence_obwfn_ << keyval->describedclassvalue("valence_obwfn");
    if (valence_obwfn_) {
      double valence_obwfn_energy = valence_obwfn_->energy();//compute CLHF guess wavefunc
      this->reference()->compute();// compute reference, e.g. PsiCLHF
      moorder_ = detail::ref_to_valence_moorder(this->reference(),
                                                this->valence_obwfn_);
    }
    else {
      if(keyval->exists("moorder")) {
        const int nmo = this->reference()->nmo();
        if(keyval->count() != nmo)
          throw InputError("Input error for PsiRASCI -- moorder must have nmo elements.",__FILE__,__LINE__);
        moorder_ = this->read_occ(keyval, "moorder", nmo);
      }
    }

    ras1_ = read_occ(keyval,"ras1",nirrep_);
    ras2_ = read_occ(keyval,"ras2",nirrep_);
    ras3_ = read_occ(keyval,"ras3",nirrep_);
    ras1_max_ = keyval->intvalue("ras1_max", KeyValValueint(2));
    ras3_max_ = keyval->intvalue("ras3_max",KeyValValueint(2));

    vector<unsigned int> docc_act = this->docc_act();
    vector<unsigned int> socc = this->socc();
    vector<unsigned int> uocc_act = this->uocc_act();//these 3 lines are redundant.
    vector<unsigned int> mos = reference()->mopi();
    vector<unsigned int> frozen_docc = this->frozen_docc();
    vector<unsigned int> frozen_uocc = this->frozen_uocc();
    const unsigned int nirrep = reference()->nirrep();

    /// if ras1 is not specified, it is empty
    if(ras1_.empty()) {
      ras1_.resize(nirrep);
      std::fill(ras1_.begin(), ras1_.end(), 0);
    }

    /// if ras2 is not specified, look for valence_obwfn ...
    if(ras2_.empty()) {
      if (valence_obwfn_) { // if given, compute the valence space
        ras2_.resize(nirrep);
        // ras2_[i] = valence_obwfn_.orbs[i] - frozen_docc[i] - ras1_[i];
        for(unsigned int h=0; h<nirrep; ++h)
          ras2_[h] = valence_obwfn_->oso_dimension()->blocks()->size(h) - frozen_docc[h] - ras1_[h];
      }
      else
        InputError("valence_obwfn keyword is not given, thus ras2 vector must be specified",
                   __FILE__, __LINE__, "ras2");
    }

    // if ras3 is not specified, ras3_[i] = mos[i] - frozen_docc[i] - ras1_[i] - ras2_[i] - frozen_uocc[i];
    if(ras3_.empty()) {
      ras3_.resize(nirrep);
      std::transform(mos.begin(), mos.end(), frozen_docc.begin(), ras3_.begin(),
                     std::minus<unsigned int>());
      std::transform(ras3_.begin(), ras3_.end(), ras1_.begin(), ras3_.begin(),
                     std::minus<unsigned int>());
      std::transform(ras3_.begin(), ras3_.end(), ras2_.begin(), ras3_.begin(),
                     std::minus<unsigned int>());
      std::transform(ras3_.begin(), ras3_.end(), frozen_uocc.begin(), ras3_.begin(),
                     std::minus<unsigned int>());
    }

    const int energ_acc = -log10(desired_value_accuracy());
    energy_convergence_ = keyval->intvalue("energy_convergence",KeyValValueint(energ_acc));
    convergence_ = keyval->intvalue("convergence",KeyValValueint(energy_convergence_-2));
    maxiter_ = keyval->intvalue("maxiter",KeyValValueint(100));

    scf_levelshift_ = keyval->doublevalue("scf_levelshift",KeyValValuedouble(1.0));
    scf_stop_levelshift_ = keyval->intvalue("scf_stop_levelshift",KeyValValueint(10));

  }

  PsiRASCI::PsiRASCI(StateIn &s)
    : PsiCorrWavefunction(s) {
    s.get(opdm_print_);
    s.get(tpdm_print_);
    s.get(root_);
    s.get(multiplicity_);
    s.get(nroots_);
    s.get(h0_blocksize_);
    s.get(repl_otf_);
    s.get(energy_convergence_);
    s.get(convergence_);
    s.get(maxiter_);
    s.get(ras1_);
    s.get(ras2_);
    s.get(ras3_);
    s.get(ras1_max_);
    s.get(ras3_max_);
    s.get(target_sym_);

    s.get(scf_levelshift_);
    s.get(scf_stop_levelshift_);

    valence_obwfn_ << SavableState::restore_state(s);

    s.get(moorder_);
  }

  void PsiRASCI::save_data_state(StateOut &s) {
    PsiCorrWavefunction::save_state(s);
    s.put(opdm_print_);
    s.put(tpdm_print_);
    s.put(root_);
    s.put(multiplicity_);
    s.put(nroots_);
    s.put(h0_blocksize_);
    s.put(repl_otf_);
    s.put(energy_convergence_);
    s.put(convergence_);
    s.put(maxiter_);
    s.put(ras1_);
    s.put(ras2_);
    s.put(ras3_);
    s.put(ras1_max_);
    s.put(ras3_max_);
    s.put(target_sym_);

    s.put(scf_levelshift_);
    s.put(scf_stop_levelshift_);

    SavableState::save_state(valence_obwfn_.pointer(), s);

    s.put(moorder_);
  }

  PsiRASCI::~PsiRASCI() {}

  void PsiRASCI::write_input(int convergence) {
    Ref<PsiInput> input = get_psi_input();
    input->open();

    PsiCorrWavefunction::write_input(convergence);
    input->write_keyword("psi:wfn", "detci");
    // no need to run cscf -- the orbitals have already been determined
    input->write_keyword("psi:exec", "( \"cints\" \"transqt2\" \"detci\")");

    const bool rasscf = false;
    write_rasci_input(convergence, rasscf);

    input->close();
  }

  void PsiRASCI::write_rasci_input(int convergence, bool rasscf) {

    Ref<PsiInput> input = get_psi_input();

    input->write_keyword("psi:jobtype","sp");
    input->write_keyword("detci:opdm","true");
    input->write_keyword("detci:tpdm","true");

    input->write_keyword("detci:opdm_print",opdm_print_);
    input->write_keyword("detci:tpdm_print",tpdm_print_);

    input->write_keyword("detci:root",root_);
    input->write_keyword("detci:s",(multiplicity_ - 1)/2);
    // If computing a non-singlet state using closed-shell as a reference need to tell detci to solve for M_S = 0 state
    const bool m_s_eq_0 = this->reference()->spin_polarized() == false;
    input->write_keyword("detci:ms0",m_s_eq_0);
    input->write_keyword("detci:calc_ssq",true);
    input->write_keyword("detci:num_roots",nroots_);
    if (target_sym_ != -1 && !rasscf)
      input->write_keyword("detci:ref_sym",target_sym_);

    if(h0_blocksize_ > 0) { /// use "h0_blocksize" keyword only if it is >0
      input->write_keyword("detci:h0_blocksize",h0_blocksize_);
    }

    if((repl_otf_==true) && !rasscf) {  /// don't use "repl_otf" keyword for RASSCF calculations.
      input->write_keyword("detci:repl_otf","true");
    }

    input->write_keyword("detci:guess_vector", "h0_block");
    input->write_keyword("detci:export_vector", "true");
    input->write_keyword("detci:energy_convergence", energy_convergence_);
    input->write_keyword("detci:convergence", convergence_);
    input->write_keyword("detci:maxiter", maxiter_);

    if (!ras1_.empty()) input->write_keyword_array("psi:ras1",ras1_);
    if (!ras2_.empty()) input->write_keyword_array("psi:ras2",ras2_);
    if (!ras3_.empty()) input->write_keyword_array("psi:ras3",ras3_);
    input->write_keyword("detci:ex_lvl",ras1_max_);
    input->write_keyword("psi:ras3_max",ras3_max_);

    input->write_keyword("scf:levelshift",scf_levelshift_);
    input->write_keyword("scf:stop_levelshift",scf_stop_levelshift_);

    if(!moorder_.empty()) {
      input->write_keyword("scf:reorder", "after");
      input->write_keyword_array("scf:moorder", moorder_);
    }
  }

  void PsiRASCI::compute() {

    // compute reference NOW so that orbital info is available for making the correlated wfn input
    reference_->compute();
    const double energy_scf = reference_energy();
    PsiWavefunction::compute();

    const double energy_ci = exenv()->chkpt().rd_etot();
    const double energy_rasscf = 0.0;

    ExEnv::out0() << indent << "SCF energy: " << setprecision(12) << energy_scf << endl;
    ExEnv::out0() << indent << "correlation (CI-SCF) energy: " << setprecision(12) << (energy_ci-energy_scf) << endl;
    ExEnv::out0() << indent << "total CI energy: " << setprecision(12) << energy_ci << endl;

    set_energy(energy_ci);
  }

  void PsiRASCI::print(std::ostream& os) const {
    os << indent << "PsiRASCI:" << endl;
    os << incindent;
    detail::print_blocks("frzocc",this->frozen_docc(),os);
    detail::print_blocks("ras1",ras1_,os);
    detail::print_blocks("ras2",ras2_,os);
    detail::print_blocks("ras3",ras3_,os);
    os << indent << "ras1_max = " << ras1_max_ << std::endl;
    os << indent << "ras3_max = " << ras3_max_ << std::endl;
    PsiCorrWavefunction::print(os);
    os << decindent;
  }

  std::vector<unsigned int>
  PsiRASCI::map_density_to_sb() {
    // maps symm -> RAS
    std::vector<unsigned int> fmap = index_map_symmtocorrorder(frozen_docc(),
                                                             ras1(),
                                                             ras2(),
                                                             ras3(),
                                                             frozen_uocc());
    return index_map_inverse( fmap );
  }

  const Ref<OrbitalSpace>&  PsiRASCI::orbs_sb(SpinCase1 spin) {
    if(orbs_sb_[spin]) {
      return orbs_sb_[spin];
    }
    if (spin==Beta) {
      return this->orbs_sb(Alpha);
    }

    // reread orbital coefficient and energies in case this is an MCSCF
    PsiChkpt chkpt(exenv(), integral(), debug());
    const std::string id = ParsedOrbitalSpaceKey::key(std::string("pp(sym)"),spin);
    orbs_sb_[spin] = new OrbitalSpace(id,prepend_spincase(spin,"MOs (Psi3)"),
                                          chkpt.coefs(spin),basis(),integral(),
                                          chkpt.evals(spin),0,0,OrbitalSpace::symmetry);

    return orbs_sb_[spin];
  }

  const Ref<OrbitalSpace>& PsiRASCI::occ(SpinCase1 spin) {
    const Ref<OrbitalSpace>& orbs = this->orbs_sb(spin);
    if (ras3_max() > 0)
      return orbs;
    else { // ras3 is empty -> select frozen_docc + ras1 + ras2
      if(occ_[spin]) {
        return occ_[spin];
      }
      if (spin==Beta) {
        return this->occ(Alpha);
      }

      std::vector<bool> occ_mask(orbs->rank(), false);
      const int nirrep = this->nirrep();
      std::vector<unsigned int> mopi = this->reference()->mopi();
      std::vector<bool>::iterator mo = occ_mask.begin();
      for(int h=0; h<nirrep; ++h)
      {
        const int nocc = frozen_docc_[h] + ras1_[h] + ras2_[h];
        std::fill(mo, mo+nocc, true);
        mo += mopi[h];
      }
      std::string label = prepend_spincase(spin, "PsiRASCI occupied MOs");
      std::string id = ParsedOrbitalSpaceKey::key(std::string("x(sym)"),spin);
      occ_[spin] =
          new MaskedOrbitalSpace(id, label, orbs, occ_mask);

      return occ_[spin];
    }
  }

  RefSymmSCMatrix PsiRASCI::mo_density(SpinCase1 spin) {
    if (ras3_max_ > 0)
      return PsiCorrWavefunction::mo_density(spin);
    else {
      RefSymmSCMatrix opdm_occ = this->onepdm_occ(spin);
      // if density is already in required spaces?
      Ref<OrbitalSpace> orbs_sb = this->orbs_sb(spin);
      Ref<OrbitalSpace> occ = this->occ(spin);
      if (*occ == *orbs_sb)
        return opdm_occ;

      RefSymmSCMatrix result = this->basis_matrixkit()->symmmatrix(orbs_sb->dim());
      result.assign(0.0);
      std::vector<int> omap = map(*occ, *orbs_sb);
      const int nmo = orbs_sb->rank();

      for(int R=0; R<nmo; ++R)
        for(int C=0; C<=R; ++C) { //lower triangle
          const int rr = omap[R];
          const int cc = omap[C];
          if (rr == -1 || cc == -1) continue;
          const double rdm_R_C = opdm_occ.get_element(rr, cc);
          result.set_element(R, C, rdm_R_C);
        }

      return result;
    }
    MPQC_ASSERT(false); // unreachable
  }

  RefSymmSCMatrix PsiRASCI::onepdm_occ(SpinCase1 spin) {

    if (ras3_max_ > 0 && this->nfzv() == 0) return this->mo_density(spin);

    if ((spin == Beta || spin == AnySpinCase1) && !this->spin_polarized())
      return onepdm_occ(Alpha);

    if (spin == AnySpinCase1 && this->spin_polarized())
      ProgrammingError("asked for any spin density but the density is spin-polarized");

    if (onepdm_occ_[spin])
      return onepdm_occ_[spin];

    // ensure that this has been computed
    { const double energy = this->value(); }

    // What is about to happen:
    // 1-rdm reported by Psi is in RAS order, hence need to map it to symmetry-blocked order
    // Psi reports density in terms of occupied + CI orbitals defined as frozen_docc + ras1 + ras2 + ras3
    // this does not depend on whether ras3_max_ > 0 or not! Since we are here this means ras3_max_ = 0.
    // thus, will first extract density in FRZOCC+RAS1+RAS2+RAS3, then filter them out RAS3 if necessary

    // here need to map RAS1+RAS2+RAS3 to full SB
    std::vector<unsigned int> dmap1;
    {
      std::vector<unsigned int> empty(nirrep_, 0);
      std::vector<unsigned int> fmap = index_map_symmtocorrorder(frozen_docc(),
                                                               ras1(),
                                                               ras2(),
                                                               ras3(),
                                                               empty);
      dmap1 = index_map_inverse( fmap );
    }
    std::vector<unsigned int> opdm_orbs_pi(nirrep_);
    for(int i=0; i<nirrep_; i++) {
      opdm_orbs_pi[i] = frozen_docc_[i] + ras1_[i] + ras2_[i] + ras3_[i];
    }
    onepdm_occ_[spin] = detail::rdopdm(spin,
                                       opdm_orbs_pi,
                                       dmap1,
                                       this->basis_matrixkit());

    const unsigned int norbs_in_ras3 = std::accumulate(ras3_.begin(), ras3_.end(), 0u);
    if (norbs_in_ras3 > 0) { // filter out RAS3 orbitals since ras3_max is 0
      std::vector<unsigned int> opdm_occ_pi(nirrep_);
      for(int i=0; i<nirrep_; i++) {
        opdm_occ_pi[i] = frozen_docc_[i] + ras1_[i] + ras2_[i];
      }
      const unsigned int nocc = std::accumulate(opdm_occ_pi.begin(), opdm_occ_pi.end(), 0);
      RefSCDimension occdim = new SCDimension(nocc,
                                              nirrep_,
                                              reinterpret_cast<int*>(const_cast<unsigned int*>(&(opdm_occ_pi[0])))
                                            );
      for (unsigned int h=0; h<nirrep_; ++h)
        occdim->blocks()->set_subdim(h, new SCDimension(opdm_occ_pi[h]));
      RefSymmSCMatrix opdm_occ = this->basis_matrixkit()->symmmatrix(occdim);

      unsigned int occ_offset = 0;
      unsigned int orbs_offset = 0;
      for(unsigned int h=0; h<nirrep_; h++) {
        const unsigned int nocc = opdm_occ_pi[h];
        const unsigned int norbs = opdm_orbs_pi[h];
        for(unsigned int r=0; r<nocc; ++r) {
          for(unsigned int c=0; c<=r; ++c) {
            const double value = onepdm_occ_[spin].get_element( orbs_offset + r, orbs_offset + c );
            opdm_occ.set_element( occ_offset + r, occ_offset + c, value );
          }
        }
        occ_offset += nocc;
        orbs_offset += norbs;
      }
      // tada!
      onepdm_occ_[spin] = opdm_occ;
    }
    // else this is what we need

    if(debug_>=DefaultPrintThresholds::mostN2) {
      onepdm_occ_[spin].print(prepend_spincase(spin,"Psi opdm (occupied part)").c_str());
    }

    return onepdm_occ_[spin];
  }

  RefSymmSCMatrix PsiRASCI::twopdm_occ(SpinCase2 spin) {

    if (ras3_max_ > 0 && this->nfzv() == 0) return this->twopdm_dirac(spin);

    if (spin == BetaBeta && !this->spin_polarized())
      return twopdm_occ(AlphaAlpha);

    if (twopdm_occ_[spin])
      return twopdm_occ_[spin];

    // 2-rdm reported by Psi is in RAS order, hence need to map it to symmetry-blocked order
    // map_density_to_sb() reports mapping from full RAS to full SB orders
    // here need to map RAS1+RAS2 (since ras3_max = 0) to full SB
    std::vector<unsigned int> dmap;
    {
      std::vector<unsigned int> empty(nirrep_, 0);
      std::vector<unsigned int> fmap = index_map_symmtocorrorder(frozen_docc(),
                                                               ras1(),
                                                               ras2(),
                                                               empty,
                                                               empty);
      dmap = index_map_inverse( fmap );
    }
    twopdm_occ_[spin] = detail::rdtpdm(spin, dmap);

    if(debug_>=DefaultPrintThresholds::mostN4) {
      twopdm_occ_[spin].print(prepend_spincase(spin,"Psi 2-RDM").c_str());
    }

    return twopdm_occ_[spin];
  }

  RefSymmSCMatrix PsiRASCI::twopdm_occ() {

    if (ras3_max_ > 0 && this->nfzv() == 0) return this->twopdm_dirac();

    if (twopdm_sf_occ_)
      return twopdm_sf_occ_;

    // 2-rdm reported by Psi is in RAS order, hence need to map it to symmetry-blocked order
    // map_density_to_sb() reports mapping from full RAS to full SB orders
    // here need to map RAS1+RAS2 (since ras3_max = 0) to full SB
    std::vector<unsigned int> dmap;
    {
      std::vector<unsigned int> empty(nirrep_, 0);
      std::vector<unsigned int> fmap = index_map_symmtocorrorder(frozen_docc(),
                                                               ras1(),
                                                               ras2(),
                                                               empty,
                                                               empty);
      dmap = index_map_inverse( fmap );
    }
    const bool spinfree_is_true = true;
    twopdm_sf_occ_ = detail::rdtpdm(AlphaBeta, dmap, spinfree_is_true);

    if(debug_>=DefaultPrintThresholds::mostN4) {
      twopdm_sf_occ_.print("Psi spin-free 2-RDM");
    }

    return twopdm_sf_occ_;
  }

  double
  PsiRASCI::magnetic_moment() const {
    return static_cast<double>(multiplicity_ - 1);
  }

//////////////////////////////////////////////////////////

  ClassDesc
  PsiRASSCF::class_desc_(typeid(PsiRASSCF),
                         "PsiRASSCF",
                         1,               // version
                         "public PsiRASCI", // must match parent
                         0,               // change to create<PsiRASSCF> if this class is DefaultConstructible
                         create<PsiRASSCF>, // change to 0 if this class is not KeyValConstructible
                         create<PsiRASSCF>  // change to 0 if this class is not StateInConstructible
  );

  PsiRASSCF::PsiRASSCF(const Ref<KeyVal>& keyval) : PsiRASCI(keyval) {

    relax_core_ = keyval->booleanvalue("relax_core",KeyValValueboolean(false));
    state_average_ = keyval->booleanvalue("state_average",KeyValValueboolean(false));

    const int energ_acc = -log10(desired_value_accuracy());
    rasscf_energy_convergence_ = keyval->intvalue("rasscf_energy_convergence",KeyValValueint(energ_acc));
    rasscf_convergence_ = keyval->intvalue("rasscf_convergence",KeyValValueint(rasscf_energy_convergence_-2));
    diis_start_ = keyval->intvalue("diis_start",KeyValValueint(0));
    rasscf_maxiter_ = keyval->intvalue("rasscf_maxiter",KeyValValueint(200));

    // the default = -1 is to leave it out
    rasscf_target_sym_ = keyval->intvalue("rasscf_target_sym", KeyValValueint(-1));
    if (keyval->exists("rasscf_target_sym") && (rasscf_target_sym_ < 0 || rasscf_target_sym_ >= this->molecule()->point_group()->char_table().nirrep()))
      throw InputError("PsiRASSCF::PsiRASSCF(const Ref<KeyVal> &) -- rasscf_target_sym must be between 0 and number of irreps - 1",__FILE__,__LINE__);

    // reset some defaults for the RASCI class
    energy_convergence_ = keyval->intvalue("energy_convergence",KeyValValueint(rasscf_energy_convergence_+2));
    convergence_ = keyval->intvalue("convergence",KeyValValueint(rasscf_convergence_+2));
    max_step_ = keyval->doublevalue("max_step",KeyValValuedouble(0.30));
    if (keyval->exists("ras3") == false)
      std::fill(ras3_.begin(), ras3_.end(), 0);

    ras1_max_ = keyval->intvalue("ras1_max", KeyValValueint(0)); // overwrite ras1_max in PsiRASSCF to make the defaults proper
    ras3_max_ = keyval->intvalue("ras3_max", KeyValValueint(0)); // overwrite ras3_max in PsiRASSCF to make the defaults proper

  }

  PsiRASSCF::PsiRASSCF(StateIn& s) : PsiRASCI(s) {
    s.get(diis_start_);
    s.get(rasscf_maxiter_);
    s.get(relax_core_);
    s.get(state_average_);
}

  PsiRASSCF::~PsiRASSCF() {
    // this may be necessary if this is a templated class
    const bool make_sure_class_desc_initialized = (&class_desc_ != 0);
  }

  void
  PsiRASSCF::save_data_state(StateOut& s) {
    s.put(diis_start_);
    s.put(rasscf_maxiter_);
    s.put(relax_core_);
    s.put(state_average_);
  }

  void
  PsiRASSCF::compute() {

    run_detci_only_ = false;

    // compute reference NOW so that orbital info is available for making the correlated wfn input
    reference_->compute();
    PsiWavefunction::compute();
    const double energy_rasscf = exenv()->chkpt().rd_etot();
    ExEnv::out0() << indent << "RASSCF energy: " << setprecision(12) << energy_rasscf << endl;

    // if this is state-averaged RASSCF
    // run detci again as follows:
    // 1) clean, but keep the checkpoint file
    // 2) then run "cints", "transqt2", "detci"
    if (state_average_ == true) {
      // unset state_average to avoid producing densities averaged over several roots
      const bool state_average_orig = state_average_; state_average_ = false;
      run_detci_only_ = true;

      // false -> keep the checkpoint file!
      exenv()->run_psiclean(false);
      exenv()->keep_output();
      PsiWavefunction::compute();

      state_average_ = state_average_orig;
      run_detci_only_ = false;
    }

  }

  void
  PsiRASSCF::write_input(int convergence) {

    Ref<PsiInput> input = get_psi_input();
    input->open();

    PsiCorrWavefunction::write_input_frozen2restricted(convergence, relax_core_);
    input->write_keyword("psi:wfn", "rasscf");
    if (run_detci_only_)
      input->write_keyword("psi:exec", "( \"cints\" \"transqt2\" \"detci\")");

    if (nroots_ > 1 && state_average_ == true) {
      std::vector<unsigned int> average_states(nroots_);
      std::vector<double> average_weights(nroots_);
      const double weight = 1.0/nroots_;
      for(int s=0; s<nroots_; ++s) {
        average_states[s] = s+1;
        average_weights[s] = weight;
      }

      input->write_keyword_array("detci:average_states", average_states);
      input->write_keyword_array("detci:average_weights", average_weights);
    }


    input->write_keyword("detcas:convergence",rasscf_convergence_);
    input->write_keyword("detcas:energy_convergence",rasscf_energy_convergence_);
    input->write_keyword("default:ncasiter",rasscf_maxiter_);
    if(diis_start_ < 0) {
      input->write_keyword("detcas:diis_start",1000000);
    }
    else
      input->write_keyword("detcas:diis_start", diis_start_);
    if(max_step_ != 0.3)
      input->write_keyword("detcas:max_step", max_step_);

    if (rasscf_target_sym_ != -1)
      input->write_keyword("detci:ref_sym",rasscf_target_sym_);

    const bool rasscf = true;
    write_rasci_input(convergence, rasscf);

    // need to compute restricted_uocc as well
    // ruocc[i] = mos[i] - frozen_docc[i] - ras1_[i] - ras2_[i] - ras3_[i] - frozen_uocc[i];
    const unsigned int nirr = this->nirrep();
    std::vector<unsigned int> ruocc(nirr);
    std::vector<unsigned int> mos = reference()->mopi();
    vector<unsigned int> frozen_docc = this->frozen_docc();
    vector<unsigned int> frozen_uocc = this->frozen_uocc();
    for(unsigned int h=0; h<nirr; ++h)
      ruocc[h] =  mos[h] - frozen_docc[h] - ras1_[h] - ras2_[h] - ras3_[h] - frozen_uocc[h];
    input->write_keyword_array("psi:restricted_uocc",ruocc);

    input->close();

  }

  void PsiRASSCF::print(std::ostream& os) const {

    os << indent << "PsiRASSCF:" << endl;
    os << incindent;
    os << indent << "relax_core = " << (relax_core_ ? "true" : "false") << std::endl;
    PsiRASCI::print(os);
    os << decindent;

  }

} // end of namespace sc

