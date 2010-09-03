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


#ifdef __GNUC__
#pragma implementation
#endif

#include <assert.h>
#include <psifiles.h>
#include <ccfiles.h>
#include <cmath>

#include <chemistry/qc/mbptr12/pairiter.impl.h>
#include <chemistry/qc/mbptr12/print.h>
#include <chemistry/qc/mbptr12/orbitalspace_utils.h>
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
     *  of diffuse orbitals into the active space of CASSCF calculations.
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

  static ClassDesc PsiCI_cd(typeid(PsiCI), "PsiCI", 1,
                                   "public PsiCorrWavefunction",  0, create<PsiCI>,
                                   create<PsiCI>);

  PsiCI::PsiCI(const Ref<KeyVal> &keyval)
    : PsiCorrWavefunction(keyval) {

    rasscf_ = keyval->booleanvalue("rasscf",KeyValValueboolean(false));
    relax_core_ = keyval->booleanvalue("relax_core",KeyValValueboolean(false));
    state_average_ = false;
    if (rasscf_) state_average_ = keyval->booleanvalue("state_average",KeyValValueboolean(false));
    wfn_type_ = rasscf_ ? "detcas" : "detci";

    opdm_print_ = keyval->booleanvalue("opdm_print", KeyValValueboolean(false));
    tpdm_print_ = keyval->booleanvalue("tpdm_print", KeyValValueboolean(false));

    root_ = keyval->intvalue("root", KeyValValueint(1));
    detci_num_roots_ = keyval->intvalue("detci_num_roots", KeyValValueint(root_));
    if(detci_num_roots_<root_) {
      throw InputError("PsiCI::PsiCI(const Ref<KeyVal> &) -- detci_num_roots should not be smaller than root.",__FILE__,__LINE__);
    }
    detcas_detci_num_roots_ = keyval->intvalue("detcas_detci_num_roots", KeyValValueint(root_));
    if(detcas_detci_num_roots_<root_) {
      throw InputError("PsiCI::PsiCI(const Ref<KeyVal> &) -- detcas_detci_num_roots should not be smaller than root.",__FILE__,__LINE__);
    }
    // the default = -1 is to leave it out
    detci_ref_sym_ = keyval->intvalue("detci_ref_sym", KeyValValueint(-1));
    if (keyval->exists("detci_ref_sym") && (detci_ref_sym_ < 0 || detci_ref_sym_ >= this->molecule()->point_group()->char_table().nirrep()))
      throw InputError("PsiCI::PsiCI(const Ref<KeyVal> &) -- detci_ref_sym must be between 0 and number of irreps - 1",__FILE__,__LINE__);

    h0_blocksize_ = keyval->intvalue("h0_blocksize", KeyValValueint(40));
    ex_lvl_ = keyval->intvalue("ex_lvl", KeyValValueint(2));
    repl_otf_ = keyval->booleanvalue("repl_otf", KeyValValueboolean(false));

    ras1_ = read_occ(keyval,"ras1",nirrep_);
    ras2_ = read_occ(keyval,"ras2",nirrep_);
    ras3_ = read_occ(keyval,"ras3",nirrep_);
    ras3_max_ = keyval->intvalue("ras3_max",KeyValValueint(ex_lvl_));

    // Psi has not done a correlated calculation, these will fail
    //assert(false);
    vector<unsigned int> docc_act = this->docc_act();
    vector<unsigned int> socc = this->socc();
    vector<unsigned int> uocc_act = this->uocc_act();
    vector<unsigned int> mos = reference()->mopi();
    vector<unsigned int> frozen_docc = this->frozen_docc();
    vector<unsigned int> frozen_uocc = this->frozen_uocc();
    const unsigned int nirrep = reference()->nirrep();

    /// if ras1 is not given in input, it is docc_act + socc
    if(ras1_.empty()) {
      ras1_.resize(nirrep);
      std::transform(docc_act.begin(), docc_act.end(), socc.begin(), ras1_.begin(),
                     std::plus<unsigned int>());
    }

    /// if ras2 is not given in input, it is empty
    if(ras2_.empty()) {
      ras2_.resize(nirrep);
      std::fill(ras2_.begin(), ras2_.end(), 0);
    }

    // if ras3 empty, ras3_[i] = mos[i] - frozen_docc[i] - ras1_[i] - ras2_[i] - frozen_uocc[i];
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

    int energ_acc = -log10(desired_value_accuracy());
    detci_energy_convergence_ = keyval->intvalue("detci_energy_convergence",KeyValValueint(energ_acc));
    detci_convergence_ = keyval->intvalue("detci_convergence",KeyValValueint(9));
    detcas_energy_convergence_ = keyval->intvalue("detcas_energy_convergence",KeyValValueint(detci_energy_convergence_-2));
    detcas_convergence_ = keyval->intvalue("detcas_convergence",KeyValValueint(detci_convergence_-3));
    detcas_detci_energy_convergence_ = keyval->intvalue("detcas_detci_energy_convergence",KeyValValueint(detcas_energy_convergence_+2));
    detcas_detci_convergence_ = keyval->intvalue("detcas_detci_convergence",KeyValValueint(detcas_convergence_+2));
    if(keyval->exists("detcas_diis")) {
      detcas_diis_ = keyval->booleanvalue("detcas_diis");
      detcas_diis_start_ = keyval->intvalue("detcas_diis_start");
    }
    else {
      detcas_diis_ = false;
    }
    detcas_detci_maxiter_ = keyval->intvalue("detcas_detci_maxiter",KeyValValueint(5));
    detci_maxiter_ = keyval->intvalue("detci_maxiter",KeyValValueint(100));
    detcas_maxiter_ = keyval->intvalue("detcas_maxiter",KeyValValueint(200));

    reorder_ = keyval->stringvalue("reorder", KeyValValuestring("after"));

    scf_levelshift_ = keyval->doublevalue("scf_levelshift",KeyValValuedouble(0.1));
    scf_stop_levelshift_ = keyval->intvalue("scf_stop_levelshift",KeyValValueint(100));

    /// construction of valence_obwfn_
    valence_obwfn_ << keyval->describedclassvalue("valence_obwfn");
    if (valence_obwfn_.nonnull()) {
      if(reorder_!="after")
        throw InputError("PsiCI::PsiCI -- if valence_obwfn is specified, only reorder=after makes sense",__FILE__,__LINE__);

      if(!keyval->exists("moorder")) {
        double valence_obwfn_energy = valence_obwfn_->energy();
        moorder_ = detail::ref_to_valence_moorder(this->reference(),
                                                  this->valence_obwfn_);
      }

    }

    if(!reorder_.empty()) {
      if(keyval->exists("moorder")) {
        const int nmo = this->reference()->nmo();
        if(keyval->count() != nmo)
          throw InputError("Input error for PsiCI -- moorder must have nmo elements.",__FILE__,__LINE__);
        moorder_ = this->read_occ(keyval, "moorder", nmo);
      }
    }

    run_detci_only_ = false;
  }

  PsiCI::PsiCI(StateIn &s)
    : PsiCorrWavefunction(s) {
    int opdm_print; s.get(opdm_print); opdm_print_ = (bool)opdm_print;
    int tpdm_print; s.get(tpdm_print); tpdm_print_ = (bool)tpdm_print;
    s.get(root_);
    s.get(detci_num_roots_);
    s.get(detcas_detci_num_roots_);
    s.get(h0_blocksize_);
    s.get(ex_lvl_);
    int repl_otf; s.get(repl_otf); repl_otf_ = (bool)repl_otf;
    s.get(detci_energy_convergence_);
    s.get(detcas_energy_convergence_);
    s.get(detcas_detci_energy_convergence_);
    s.get(detci_convergence_);
    s.get(detcas_convergence_);
    s.get(detcas_detci_convergence_);
    int detcas_diis; s.get(detcas_diis); detcas_diis_ = (bool)detcas_diis;
    s.get(detcas_detci_maxiter_);
    s.get(detci_maxiter_);
    s.get(rasscf_);
    s.get(relax_core_);
    s.get(state_average_);
    s.get(wfn_type_);
    s.get(ras1_);
    s.get(ras2_);
    s.get(ras3_);
    s.get(ras3_max_);
    s.get(run_detci_only_);
    s.get(detci_ref_sym_);

    s.get(scf_levelshift_);
    s.get(scf_stop_levelshift_);

    valence_obwfn_ << SavableState::restore_state(s);

    s.get(reorder_);
    s.get(moorder_);
  }

  void PsiCI::save_data_state(StateOut &s) {
    PsiCorrWavefunction::save_state(s);
    s.put((int)opdm_print_);
    s.put((int)tpdm_print_);
    s.put(root_);
    s.put(detci_num_roots_);
    s.put(detcas_detci_num_roots_);
    s.put(h0_blocksize_);
    s.put(ex_lvl_);
    s.put((int)repl_otf_);
    s.put(detci_energy_convergence_);
    s.put(detcas_energy_convergence_);
    s.put(detcas_detci_energy_convergence_);
    s.put(detci_convergence_);
    s.put(detcas_convergence_);
    s.put(detcas_detci_convergence_);
    s.put((int)detcas_diis_);
    s.put(detcas_detci_maxiter_);
    s.put(detci_maxiter_);
    s.put(rasscf_);
    s.put(relax_core_);
    s.put(state_average_);
    s.put(wfn_type_);
    s.put(ras1_);
    s.put(ras2_);
    s.put(ras3_);
    s.put(ras3_max_);
    s.put(run_detci_only_);
    s.put(detci_ref_sym_);

    s.put(scf_levelshift_);
    s.put(scf_stop_levelshift_);

    SavableState::save_state(valence_obwfn_.pointer(), s);

    s.put(reorder_);
    s.put(moorder_);
  }

  PsiCI::~PsiCI() {}

  void PsiCI::write_input(int convergence) {
    Ref<PsiInput> input = get_psi_input();
    input->open();

    // are we running detcas? logic is impossibly complicated at the moment -- it sucks
    // 1) if rasscf = false -> no orbital optimization, no need to run detcas
    // 2) ras3_max > 0 -> probably doing MRCI from RAS orbitals, will not optimize the orbitals
    // 3) if this is extra CI step using RAS orbitals -> will not run detcas either
    const bool do_ci_only = (rasscf_==false) || (ras3_max_>0) || (run_detci_only_==true);

    if(do_ci_only) {  // running ci without orbital optimization of core orbitals
      PsiCorrWavefunction::write_input(convergence);
    }
    else { // cas with restricted orbitals
      PsiCorrWavefunction::write_input_frozen2restricted(convergence, true);
    }

    input->write_keyword("psi:wfn",wfn_type_.c_str());
    input->write_keyword("psi:jobtype","sp");

    /// one- and two-particle density matrices have to be evaluated in any case
    input->write_keyword("detci:opdm","true");
    input->write_keyword("detci:tpdm","true");

    if(opdm_print_==true) {
      input->write_keyword("detci:opdm_print","true");
    }
    if(tpdm_print_==true) {
      input->write_keyword("detci:tpdm_print","true");
    }

    if(do_ci_only) {  // in detci calculations
      input->write_keyword("detci:root",root_);
      input->write_keyword("detci:num_roots",detci_num_roots_);
      if (detci_ref_sym_ != -1) input->write_keyword("detci:ref_sym",detci_ref_sym_);
    }
    else {  // in detcas calculations
      input->write_keyword("detci:root",root_);
      input->write_keyword("detci:num_roots",detcas_detci_num_roots_);
      if (detcas_detci_num_roots_ > 1 && state_average_ == true) {
        std::vector<unsigned int> average_states(detcas_detci_num_roots_);
        std::vector<double> average_weights(detcas_detci_num_roots_);
        const double weight = 1.0/detcas_detci_num_roots_;
        for(int s=0; s<detcas_detci_num_roots_; ++s) {
          average_states[s] = s+1;
          average_weights[s] = weight;
        }

        input->write_keyword_array("detci:average_states", average_states);
        input->write_keyword_array("detci:average_weights", average_weights);
      }
    }


    if(h0_blocksize_ > 0) { /// use "h0_blocksize" keyword only if it is >0
      input->write_keyword("detci:h0_blocksize",h0_blocksize_);
    }

    input->write_keyword("detci:ex_lvl",ex_lvl_);

    if((repl_otf_==true) && do_ci_only) {  /// don't use "repl_otf" keyword for CASSCF calculations.
      input->write_keyword("detci:repl_otf","true");
    }

    input->write_keyword("detci:guess_vector","h0_block");

    input->write_keyword("detci:export_vector","true");

    // RAS stuff
    if(rasscf_ && (ras3_max_==0) && wfn_type_ == std::string("casscf")) {
      if (!ras2_.empty()) input->write_keyword_array("psi:active",ras2_);
      input->write_keyword("detci:energy_convergence",detcas_detci_energy_convergence_);
      input->write_keyword("detcas:energy_convergence",detcas_energy_convergence_);
      input->write_keyword("detci:convergence",detcas_detci_convergence_);
      input->write_keyword("detcas:convergence",detcas_convergence_);
      input->write_keyword("default:ncasiter",detcas_maxiter_);
      if(detcas_diis_==false) {
        input->write_keyword("detcas:diis_start",1000000);
      }
      else input->write_keyword("detcas:diis_start", detcas_diis_start_);
      input->write_keyword("detci:maxiter",detcas_detci_maxiter_);
    }
    else {
      if (!ras1_.empty()) input->write_keyword_array("psi:ras1",ras1_);
      if (!ras2_.empty()) input->write_keyword_array("psi:ras2",ras2_);
      if (!ras3_.empty()) input->write_keyword_array("psi:ras3",ras3_);
      input->write_keyword("psi:ras3_max",ras3_max_);
      input->write_keyword("detci:energy_convergence",detci_energy_convergence_);
      input->write_keyword("detci:convergence",detci_convergence_);
      input->write_keyword("detci:maxiter",detci_maxiter_);
    }

    // if doing rasscf but wfn_type_ == "detci" means we are using rasscf orbitals in checkpoint, hence don't run input or cscf
    if (do_ci_only) {
      input->write_keyword("psi:exec","(\"cints\" \"transqt2\" \"detci\")");
    }

    input->write_keyword("scf:levelshift",scf_levelshift_);
    input->write_keyword("scf:stop_levelshift",scf_stop_levelshift_);

    if(!moorder_.empty()) {
      input->write_keyword("scf:reorder",reorder_.c_str());
      input->write_keyword_array("scf:moorder",moorder_);
    }

    input->close();
  }

  void PsiCI::compute() {
    double energy_rasscf = 0.0;
    if (rasscf_) {

      // to compute RASSCF wavefunction temporarily customize the object so that the standard ::write_input works appropriately
      const int ras3_max_orig = ras3_max_;  ras3_max_ = 0;
      const std::string wfn_type_orig = wfn_type_;  wfn_type_ = "casscf";
      // ras1 orbitals in CAS should be treated same way as frozen orbitals
      const std::vector<unsigned int> frozen_docc_orig(frozen_docc_);
      const std::vector<unsigned int> ras1_orig(ras1_);
      {
        std::transform(frozen_docc_orig.begin(), frozen_docc_orig.end(), ras1_orig.begin(), frozen_docc_.begin(),
                       std::plus<unsigned int>());
        std::fill(ras1_.begin(), ras1_.end(), 0);
      }

      // compute RASSCF
      PsiWavefunction::compute();
      energy_rasscf = exenv()->chkpt().rd_etot();
      ExEnv::out0() << indent << "rasscf energy for the impatient: " << setprecision(12) << energy_rasscf << endl;

      // unwind the changes made previously
      ras3_max_ = ras3_max_orig;
      wfn_type_ = wfn_type_orig;
      frozen_docc_ = frozen_docc_orig;
      ras1_ = ras1_orig;

      // if ras3_max > 0 (i.e. user wants MRCI on top of RASSCF)
      //    OR
      //    this is state-averaged RASSCF
      // run detci again as follows:
      //
      // clean, but keep the checkpoint file
      // then run "cints", "transqt2", "detci"
      //
      // the difference between the 2 cases is wfn type:
      // detci in the former
      // casscf in the latter (because the densities must be expressed in active space only... ugh!)
      if (ras3_max_ > 0 ||
          (ras3_max_ == 0 && state_average_ == true)
         ) {
        // unset state_average to avoid producing densities averaged over several roots
        const bool state_average_orig = state_average_; state_average_ = false;
        const std::string wfn_type_orig = wfn_type_;
        if (ras3_max_ > 0) // need wfn = detci
          wfn_type_ = std::string("detci");
        else // this is rasscf -> wfn = casscf
          wfn_type_ = std::string("casscf");
        run_detci_only_ = true;

        // false -> keep the checkpoint file!
        exenv()->run_psiclean(false);
        PsiWavefunction::compute();

        run_detci_only_ = false;
        state_average_ = state_average_orig;
        wfn_type_ = wfn_type_orig;
      }
    }
    else {
      PsiWavefunction::compute();
    }

    const double energy_scf = reference_energy();
    const double energy_ci = exenv()->chkpt().rd_etot();

    ExEnv::out0() << indent << "SCF energy: " << setprecision(12) << energy_scf << endl;
    ExEnv::out0() << indent << "RASSCF energy: " << setprecision(12) << energy_rasscf << endl;
    ExEnv::out0() << indent << "correlation (CI-SCF) energy: " << setprecision(12) << (energy_ci-energy_scf) << endl;
    ExEnv::out0() << indent << "total CI energy: " << setprecision(12) << energy_ci << endl;

    set_energy(energy_ci);
  }

  void PsiCI::print(std::ostream& os) const {
    os << indent << "PsiCI:" << endl;
    os << incindent;
    detail::print_blocks("frzocc",this->frozen_docc(),os);
    detail::print_blocks("ras1",ras1_,os);
    detail::print_blocks("ras2",ras2_,os);
    detail::print_blocks("ras3",ras3_,os);
    os << indent << "ras3_max = " << ras3_max_ << std::endl;
    os << indent << "rasscf = " << (rasscf_ ? "true" : "false") << std::endl;
    os << indent << "relax_core = " << (relax_core_ ? "true" : "false") << std::endl;
    PsiCorrWavefunction::print(os);
    os << decindent;
  }

  std::vector<unsigned int>
  PsiCI::map_density_to_sb() {
    // maps symm -> RAS
    std::vector<unsigned int> fmap = index_map_symmtoqtorder(frozen_docc(),
                                                             ras1(),
                                                             ras2(),
                                                             ras3(),
                                                             frozen_uocc());
    return index_map_inverse( fmap );
  }

  const Ref<OrbitalSpace>&  PsiCI::orbs_sb(SpinCase1 spin) {
    if(orbs_sb_[spin].nonnull()) {
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

  const Ref<OrbitalSpace>& PsiCI::occ(SpinCase1 spin) {
    const Ref<OrbitalSpace>& orbs = this->orbs_sb(spin);
    if (ras3_max() > 0)
      return orbs;
    else { // ras3 is empty -> select frozen_docc + ras1 + ras2
      if(occ_[spin].nonnull()) {
        return occ_[spin];
      }
      if (spin==Beta) {
        return this->occ(Alpha);
      }

      std::vector<bool> occ_mask(orbs->rank(), false);
      const int nirrep = this->nirrep();
      std::vector<unsigned int> mopi = this->reference()->mopi();
      std::vector<bool>::iterator mo = occ_mask.begin();
      for(int h=0; h<nirrep; ++h) {
        const int nocc = frozen_docc_[h] + ras1_[h] + ras2_[h];
        std::fill(mo, mo+nocc, true);
        mo += mopi[h];
      }
      std::string label = prepend_spincase(spin, "PsiCI occupied MOs");
      std::string id = ParsedOrbitalSpaceKey::key(std::string("x(sym)"),spin);
      occ_[spin] =
          new MaskedOrbitalSpace(id, label, orbs, occ_mask);

      return occ_[spin];
    }
  }

  RefSymmSCMatrix PsiCI::mo_density(SpinCase1 spin) {
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
        for(int C=0; C<=R; ++C) {
          const int rr = omap[R];
          const int cc = omap[C];
          if (rr == -1 || cc == -1) continue;
          const double rdm_R_C = opdm_occ.get_element(rr, cc);
          result.set_element(R, C, rdm_R_C);
        }

      return result;
    }
    assert(false); // unreachable
  }

  RefSymmSCMatrix PsiCI::onepdm_occ(SpinCase1 spin) {

    if (ras3_max_ > 0) return this->mo_density(spin);

    if ((spin == Beta || spin == AnySpinCase1) && !this->spin_polarized())
      return onepdm_occ(Alpha);

    if (spin == AnySpinCase1 && this->spin_polarized())
      ProgrammingError("asked for any spin density but the density is spin-polarized");

    if (onepdm_occ_[spin].nonnull())
      return onepdm_occ_[spin];

    // ensure that this has been computed
    { const double energy = this->value(); }

    // 1-rdm reported by Psi is in RAS order, hence need to map it to symmetry-blocked order
    // map_density_to_sb() reports mapping from full RAS to full SB orders
    // here lets map RAS (excluding RAS3 and frozen virtuals) to occ()
    std::vector<unsigned int> dmap;
    {
      std::vector<unsigned int> empty(nirrep_, 0);
      std::vector<unsigned int> fmap = index_map_symmtoqtorder(frozen_docc(),
                                                               ras1(),
                                                               ras2(),
                                                               empty,
                                                               empty);
      dmap = index_map_inverse( fmap );
    }
    std::vector<unsigned int> occpi(nirrep_);
    for(int i=0; i<nirrep_; i++) {
      occpi[i] = frozen_docc_[i] + ras1_[i] + ras2_[i];
    }
    onepdm_occ_[spin] = detail::rdopdm(spin,
                                       occpi,
                                       dmap,
                                       this->basis_matrixkit());

    if(debug_>=DefaultPrintThresholds::mostN2) {
      onepdm_occ_[spin].print(prepend_spincase(spin,"Psi opdm").c_str());
    }

    return onepdm_occ_[spin];
  }

  RefSymmSCMatrix PsiCI::twopdm_occ(SpinCase2 spin) {

    if (ras3_max_ > 0) return this->twopdm_dirac(spin);

    if (spin == BetaBeta && !this->spin_polarized())
      return twopdm_occ(AlphaAlpha);

    if (twopdm_occ_[spin].nonnull())
      return twopdm_occ_[spin];

    // 2-rdm reported by Psi is in RAS order, hence need to map it to symmetry-blocked order
    // map_density_to_sb() reports mapping from full RAS to full SB orders
    // here lets map RAS (excluding RAS3 and frozen virtuals) to occ()
    std::vector<unsigned int> dmap;
    {
      std::vector<unsigned int> empty(nirrep_, 0);
      std::vector<unsigned int> fmap = index_map_symmtoqtorder(frozen_docc(),
                                                               ras1(),
                                                               ras2(),
                                                               empty,
                                                               empty);
      dmap = index_map_inverse( fmap );
    }
    twopdm_occ_[spin] = detail::rdtpdm(spin, dmap);

    if(debug_>=DefaultPrintThresholds::mostN4) {
      twopdm_occ_[spin].print(prepend_spincase(spin,"Psi tpdm").c_str());
    }

    return twopdm_occ_[spin];
  }

#if 0
  RefSymmSCMatrix
  PsiCI::rdm1(const SpinCase1 &spin) const {
    string opdm_label_str;
    FILE *outfile;
    if(spin==Alpha) {
      opdm_label_str = "MO-basis Alpha OPDM";
      outfile = fopen("psiout_opdm_Alpha","w");
    }
    else { // spin==Beta
      opdm_label_str = "MO-basis Beta OPDM";
      outfile = fopen("psiout_opdm_Beta","w");
    }
    char *opdm_label=(char *)opdm_label_str.c_str();
    const int nmo = reference()->nmo();
    double *opdm_arr = new double [nmo*nmo];
    detail::rdopdm(PSIF_MO_OPDM,opdm_label,opdm_arr,nmo*nmo,0,1,outfile);
    fclose(outfile);

    const RefDiagSCMatrix evals = reference()->evals(spin);
    const RefSCDimension modim = evals.dim();
    RefSymmSCMatrix opdm_mat = evals.kit()->symmmatrix(modim);
    opdm_mat.assign(opdm_arr);
    delete [] opdm_arr;

    if(debug_>=DefaultPrintThresholds::mostN2) {
      opdm_mat.print(prepend_spincase(spin,"Psi opdm").c_str());
    }

    return(opdm_mat);
  }
#endif

#if 0
  void PsiCI::print_all_blocks(ostream &o) {
    vector<unsigned int> mos = mo_blocks();
    detail::print_blocks("mos",mos,o);
    vector<unsigned int> frzc = frzcpi();
    detail::print_blocks("frzc",frzc,o);
    vector<unsigned int> docc = doccpi();
    detail::print_blocks("docc",docc,o);
    vector<unsigned int> docc_act = doccpi_act();
    detail::print_blocks("docc_act",docc_act,o);
    vector<unsigned int> socc = soccpi();
    detail::print_blocks("socc",socc,o);
    vector<unsigned int> uocc_act = uoccpi_act();
    detail::print_blocks("uocc_act",uocc_act,o);
    vector<unsigned int> uocc = uoccpi();
    detail::print_blocks("uocc",uocc,o);
    vector<unsigned int> frzv = frzvpi();
    detail::print_blocks("frzv",frzv,o);
  }
#endif

} // end of namespace sc

