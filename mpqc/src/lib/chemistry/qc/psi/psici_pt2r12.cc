//
// psici_pt2r12.cc
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

#include <math/scmat/local.h>
#include <chemistry/qc/mbptr12/pairiter.impl.h>
#include <chemistry/qc/mbptr12/print.h>
#include <chemistry/qc/psi/psiqtorder.h>
#include <chemistry/qc/psi/psiwfn.h>
#include <chemistry/qc/psi/psici_pt2r12.h>
#include <chemistry/qc/basis/split.h>
#include <chemistry/qc/scf/clhf.h>
#include <chemistry/qc/basis/petite.h>
#include <functional>
#include <algorithm>

using namespace std;
using std::map;

namespace sc {

  namespace psici_pt2r12 {

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
      for(int i=0; i<dims.size(); i++) {
        ExEnv::out0() << setw(5) << dims[i];
      }
      ExEnv::out0() << endl;

      vector<int> dim_offsets = get_block_dim_offsets(dim);
      for(int i=0; i<dim_offsets.size(); i++) {
        ExEnv::out0() << setw(5) << dim_offsets[i];
      }
      ExEnv::out0() << endl;
    }

    template <typename Int> void print_vector_int(const vector<Int> &vec_int) {
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

    /** class for an indexed overlap matrix element
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
    list<Element> create_indexed_list_from_value_matrix(const RefSCMatrix &value_mat) {
      RefSCDimension coldim = value_mat->coldim();
      RefSCDimension rowdim = value_mat->rowdim();
      int ncol_blks = coldim->blocks()->nblock();
      int nrow_blks = rowdim->blocks()->nblock();

      if(ncol_blks != nrow_blks)
        throw ProgrammingError("Error in psici_pt2r12::get_index_mat_from_value_mat -- number of block in rowdim and coldim does not coincide.",__FILE__,__LINE__);
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
    std::map<unsigned int, unsigned int> reference_val_map(const list<psici_pt2r12::Element> &sorted_elem_list) {
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

    void print_blocks(const char *orbs_name,const vector<unsigned int> &orbs_blocks,ostream &o) {
      o << setw(12) << orbs_name;
      for(int i=0; i<orbs_blocks.size(); i++) {
        o << setw(5) << orbs_blocks[i];
      }
      o << endl;
    }

    void print_blocks(const char *orbs_name,const vector<int> &orbs_blocks,ostream &o) {
      o << setw(12) << orbs_name;
      for(int i=0; i<orbs_blocks.size(); i++) {
        o << setw(5) << orbs_blocks[i];
      }
      o << endl;
    }

  }

  static ClassDesc PsiCI_PT2R12_cd(typeid(PsiCI_PT2R12), "PsiCI_PT2R12", 1,
                                   "public PsiCorrWavefunction_PT2R12",  0, create<PsiCI_PT2R12>,
                                   create<PsiCI_PT2R12>);

  PsiCI_PT2R12::PsiCI_PT2R12(const Ref<KeyVal> &keyval)
    : PsiCorrWavefunction_PT2R12(keyval) {

    rasscf_ = keyval->booleanvalue("rasscf",KeyValValueboolean(false));
    wfn_type_ = rasscf_ ? "detcas" : "detci";

    opdm_print_ = keyval->booleanvalue("opdm_print", KeyValValueboolean(false));
    tpdm_print_ = keyval->booleanvalue("tpdm_print", KeyValValueboolean(false));

    root_ = keyval->intvalue("root", KeyValValueint(1));
    detci_num_roots_ = keyval->intvalue("detci_num_roots", KeyValValueint(root_));
    if(detci_num_roots_<root_) {
      throw InputError("PsiCI_PT2R12::PsiCI_PT2R12(const Ref<KeyVal> &) -- detci_num_roots should not be smaller than root.",__FILE__,__LINE__);
    }
    detcas_detci_num_roots_ = keyval->intvalue("detcas_detci_num_roots", KeyValValueint(root_));
    if(detcas_detci_num_roots_<root_) {
      throw InputError("PsiCI_PT2R12::PsiCI_PT2R12(const Ref<KeyVal> &) -- detcas_detci_num_roots should not be smaller than root.",__FILE__,__LINE__);
    }

    h0_blocksize_ = keyval->intvalue("h0_blocksize", KeyValValueint(40));
    ex_lvl_ = keyval->intvalue("ex_lvl", KeyValValueint(2));
    repl_otf_ = keyval->booleanvalue("repl_otf", KeyValValueboolean(false));

    print_all_blocks();

    ras1_ = read_occ(keyval,"ras1",nirrep_);
    ras2_ = read_occ(keyval,"ras2",nirrep_);
    ras3_ = read_occ(keyval,"ras3",nirrep_);
    ras3_max_ = keyval->intvalue("ras3_max",KeyValValueint(ex_lvl_));

    vector<unsigned int> frozen_docc = frzc_blocks();
    vector<unsigned int> docc_act = docc_act_blocks();
    vector<unsigned int> socc = socc_blocks();
    vector<unsigned int> uocc_act = uocc_act_blocks();
    vector<unsigned int> frozen_uocc = frzv_blocks();
    vector<unsigned int> mos = mo_blocks();
    int length = frozen_docc.size();

    /// if ras1 is not given in input, it is docc_act + socc
    if(ras1_.empty()) {
      ras1_.resize(length);
      std::transform(docc_act.begin(), docc_act.end(), socc.begin(), ras1_.begin(),
                     std::plus<unsigned int>());
    }

    /// if ras2 is not given in input, it is uocc_act
    if(ras2_.empty()) {
      ras2_ = uocc_act;
    }

    // if ras3 empty, ras3_[i] = mos[i] - frozen_docc[i] - ras1_[i] - ras2_[i] - frozen_uocc[i];
    if(ras3_.empty()) {
      ras3_.resize(length);
      std::transform(mos.begin(), mos.end(), frozen_docc.begin(), ras3_.begin(),
                     std::minus<unsigned int>());
      std::transform(ras3_.begin(), ras3_.end(), ras1_.begin(), ras3_.begin(),
                     std::minus<unsigned int>());
      std::transform(ras3_.begin(), ras3_.end(), ras2_.begin(), ras3_.begin(),
                     std::minus<unsigned int>());
      std::transform(ras3_.begin(), ras3_.end(), frozen_uocc.begin(), ras3_.begin(),
                     std::minus<unsigned int>());
    }

    psici_pt2r12::print_blocks("ras1",ras1_,ExEnv::out0());
    psici_pt2r12::print_blocks("ras2",ras2_,ExEnv::out0());
    psici_pt2r12::print_blocks("ras3",ras3_,ExEnv::out0());

    int energ_acc = -log10(desired_value_accuracy());
    detci_energy_convergence_ = keyval->intvalue("detci_energy_convergence",KeyValValueint(energ_acc));
    detci_convergence_ = keyval->intvalue("detci_convergence",KeyValValueint(9));
    detcas_energy_convergence_ = keyval->intvalue("detcas_energy_convergence",KeyValValueint(detci_energy_convergence_-2));
    detcas_convergence_ = keyval->intvalue("detcas_convergence",KeyValValueint(detci_convergence_-3));
    detcas_detci_energy_convergence_ = keyval->intvalue("detcas_detci_energy_convergence",KeyValValueint(detcas_energy_convergence_+2));
    detcas_detci_convergence_ = keyval->intvalue("detcas_detci_convergence",KeyValValueint(detcas_convergence_+2));
    if(keyval->exists("detcas_diis")) {
      detcas_diis_ = keyval->booleanvalue("detcas_diis");
    }
    else {
      detcas_diis_ = false;
    }
    detcas_detci_maxiter_ = keyval->intvalue("detcas_detci_maxiter",KeyValValueint(5));
    detci_maxiter_ = keyval->intvalue("detci_maxiter",KeyValValueint(100));
    if(keyval->exists("detcas_detci_average_states")) {
      const int naverage_states = keyval->count("detcas_detci_average_states");
      detcas_detci_average_states_ = this->read_occ(keyval,"detcas_detci_average_states",naverage_states);
    }

    reorder_ = keyval->stringvalue("reorder");

    scf_levelshift_ = keyval->doublevalue("scf_levelshift",KeyValValuedouble(0.1));
    scf_stop_levelshift_ = keyval->intvalue("scf_stop_levelshift",KeyValValueint(100));

    /// construction of valence_obwfn_
    valence_obwfn_ << keyval->describedclassvalue("valence_obwfn");
    if(valence_obwfn_.nonnull()) {
      double valence_obwfn_energy = valence_obwfn_->energy();
      ExEnv::out0() << "valence one-body wavefunction energy: " << setprecision(12) << valence_obwfn_energy << endl;

      RefSCMatrix overlap_reference_wfn = overlap_with_valence_obwfn();

      if(reorder_.empty()) {
        reorder_ = "after";
      }

      if(reorder_!="after")
        throw InputError("PsiCI_PT2R12::PsiCI_PT2R12 -- if reference_wfn is specified, only reorder=after makes sense",__FILE__,__LINE__);

      if(!keyval->exists("moorder")) {
        moorder_ = map_to_valence_order(overlap_reference_wfn);
      }

    }

    if(!reorder_.empty()) {
      if(keyval->exists("moorder")) {
        const int nmo = reference_mpqc_->oso_dimension().n();
        if(keyval->count() != nmo)
          throw InputError("Input error for PsiCI_PT2R12 -- moorder must have nmo elements.",__FILE__,__LINE__);
        moorder_ = this->read_occ(keyval, "moorder", nmo);
      }
    }

  }

  PsiCI_PT2R12::PsiCI_PT2R12(StateIn &s)
    : PsiCorrWavefunction_PT2R12(s) {
    s.get(eci_);
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
    s.get(detcas_detci_average_states_);
    int rasscf; s.get(rasscf); rasscf_ = (bool)rasscf;
    s.get(wfn_type_);
    s.get(ras1_);
    s.get(ras2_);
    s.get(ras3_);
    s.get(ras3_max_);

    s.get(scf_levelshift_);
    s.get(scf_stop_levelshift_);

    valence_obwfn_ << SavableState::restore_state(s);

    s.get(reorder_);
    s.get(moorder_);
  }

  void PsiCI_PT2R12::save_data_state(StateOut &s) {
    PsiCorrWavefunction_PT2R12::save_state(s);
    s.put(eci_);
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
    s.put(detcas_detci_average_states_);
    s.put((int)rasscf_);
    s.put(wfn_type_);
    s.put(ras1_);
    s.put(ras2_);
    s.put(ras3_);
    s.put(ras3_max_);

    s.put(scf_levelshift_);
    s.put(scf_stop_levelshift_);

    SavableState::save_state(valence_obwfn_.pointer(), s);

    s.put(reorder_);
    s.put(moorder_);
  }

  PsiCI_PT2R12::~PsiCI_PT2R12() {}

  void PsiCI_PT2R12::write_input(int convergence) {
    Ref<PsiInput> input = get_psi_input();
    input->open();

    PsiCorrWavefunction_PT2R12::write_input(convergence);

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

    if((!rasscf_) || (ras3_max_>0)) {  // in detci calculations
      input->write_keyword("detci:root",root_);
      input->write_keyword("detci:num_roots",detci_num_roots_);
    }
    else {  // in detcas calculations
      input->write_keyword("detci:root",root_);
      input->write_keyword("detci:num_roots",detcas_detci_num_roots_);
    }


    if((h0_blocksize_ > 0) && ((!rasscf_) || (ras3_max_>0))) { /// use "h0_blocksize" keyword only if it is >0 and only for CI calculations, not CASSCF.
      input->write_keyword("detci:h0_blocksize",h0_blocksize_);
    }

    input->write_keyword("detci:ex_lvl",ex_lvl_);

    if((repl_otf_==true) && ((!rasscf_) || (ras3_max_>0))) {  /// don't use "repl_otf" keyword for CASSCF calculations.
      input->write_keyword("detci:repl_otf","true");
    }

    input->write_keyword("detci:guess_vector","h0_block");

    input->write_keyword("detci:export_vector","true");

    // RAS stuff
    ExEnv::out0() << "ras3_max_ = " << ras3_max_ << endl;
    ExEnv::out0() << "rasscf_ = " << ((rasscf_==true) ? "true" : "false") << endl;
    if(rasscf_ && (ras3_max_==0)) {
      if (!ras2_.empty()) input->write_keyword_array("psi:active",ras2_);
      input->write_keyword("detci:energy_convergence",detcas_detci_energy_convergence_);
      input->write_keyword("detcas:energy_convergence",detcas_energy_convergence_);
      input->write_keyword("detci:convergence",detcas_detci_convergence_);
      input->write_keyword("detcas:convergence",detcas_convergence_);
      input->write_keyword("default:ncasiter",200);
      if(detcas_diis_==false) {
        input->write_keyword("detcas:diis_start",1000000);
      }
      input->write_keyword("detci:maxiter",detcas_detci_maxiter_);
      if(!detcas_detci_average_states_.empty()) {
        input->write_keyword_array("detci:average_states",detcas_detci_average_states_);
      }
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
    if (rasscf_ && wfn_type_ == std::string("detci")) {
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

  RefSCMatrix PsiCI_PT2R12::MPQC2PSI_transform_matrix(SpinCase1 spin) {
    Ref<SCMatrixKit> localkit = new LocalSCMatrixKit;

    Ref<OrbitalSpace> orbs_mpqc = r12eval_->orbs(spin);   // MPQC symmetry-blocked orbitals
    Ref<OrbitalSpace> orbs_symm = orbs_sb(spin);          // Psi3 symmetry-blocked orbitals

    if(ras1_.empty() || ras2_.empty() || ras3_.empty())
      throw InputError("PsiCI_PT2R12::MPQC2PSI_transform_matrix -- ras1_, ras2_ and ras3_ all have to be specified.",__FILE__,__LINE__);

    std::vector<unsigned int> ras1(ras1_.size()); std::copy(ras1_.begin(), ras1_.end(), ras1.begin());
    std::vector<unsigned int> ras2(ras2_.size()); std::copy(ras2_.begin(), ras2_.end(), ras2.begin());
    std::vector<unsigned int> ras3(ras3_.size()); std::copy(ras3_.begin(), ras3_.end(), ras3.begin());
    RefSCMatrix coeffrasorder = coeffsymmtorasorder(orbs_symm->coefs(),
                                                    this->frozen_docc(),
                                                    ras1, ras2, ras3,
                                                    this->frozen_uocc());

    const int nmo = orbs_symm->coefs().coldim().n();
    RefSCDimension aodim = orbs_symm->coefs().rowdim();

    Ref<OrbitalSpace> orbs_ras = new OrbitalSpace(string("moras"),string("MOs in ras order"),orbs_symm,coeffrasorder,
                              orbs_symm->basis());

    RefSCMatrix coeffmpqc = orbs_mpqc->coefs();
    // TODO should use blocked dimensions here, but there are too many changes to be made downstream at the moment
    // see also psiqtorder.cc
#if 1
    RefSCDimension mpqcmodim = new SCDimension(nmo);
#else
    int* nfunc_per_block = new int[1];
    nfunc_per_block[0] = nmo;
    RefSCDimension mpqcmodim = new SCDimension(nmo, 1, nfunc_per_block, "MOs in MPQC order");
    if (nmo)
      mpqcmodim->blocks()->set_subdim(0, new SCDimension(nfunc_per_block[0]));
#endif
    RefSCMatrix coeffmpqc_nb = localkit->matrix(aodim,mpqcmodim);
    for(int i=0; i<aodim.n(); i++) {
      for(int j=0; j<nmo; j++) {
        coeffmpqc_nb.set_element(i,j,coeffmpqc.get_element(i,j));
      }
    }
    Ref<OrbitalSpace> orbs_mpqc_nb = new OrbitalSpace(string("mompqcnb"),string("MOs in MPQC non-blocked order"),
                                                      orbs_mpqc,coeffmpqc_nb,orbs_mpqc->basis());

    //RefSCMatrix MPQC2PSI_tform = transform(*orbs_ras,*orbs_mpqc_nb,localkit);
    RefSCMatrix MPQC2PSI_tform = transform(*orbs_ras,*orbs_mpqc_nb,localkit);

    return(MPQC2PSI_tform);
  }

  void PsiCI_PT2R12::compute() {
    double energy_rasscf = 0.0;
    if (rasscf_) {
      const int ras3_max_orig = ras3_max_;  ras3_max_ = 0;
      const std::string wfn_type_orig = wfn_type_;  wfn_type_ = "casscf";
      PsiWavefunction::compute();
      ras3_max_ = ras3_max_orig;
      wfn_type_ = wfn_type_orig;
      energy_rasscf = exenv()->chkpt().rd_etot();
      ExEnv::out0() << "rasscf energy for the impatient: " << setprecision(12) << energy_rasscf << endl;
      exenv()->run_psi_module("psiclean");
      /// if ras3_max > 0 (i.e. user wants MRCI on top of RASSCF), rebuild the input file, then run "cinst", "transqt2", "detci"
      if (ras3_max_ > 0) {
        // use wfn=detci now
        const std::string wfn_type_orig = wfn_type_; wfn_type_ = std::string("detci");
        PsiWavefunction::compute();
        wfn_type_ = wfn_type_orig;
      }
      else {
        PsiWavefunction::compute();
      }
    }
    else {
      PsiWavefunction::compute();
    }

    const double energy_scf = reference_energy();
    double energy_conv = exenv()->chkpt().rd_etot();
    double energy_from_densities;
    double value_acc = this->desired_value_accuracy();
    // value_acc is only for check against energy_from_densities;
    // should not be too tight avoid unintentional ProgrammingError to be thrown.
    value_acc = (value_acc>=1.0e-9) ? value_acc : 1.0e-9;
    if(r12evalinfo_->r12tech()->corrfactor()->geminaldescriptor()->type()!=string("invalid")) {
      ExEnv::out0() << "conventional total energy for the impatient: " << setprecision(12) << energy_conv << endl;
      energy_from_densities = energy_conventional();
      ExEnv::out0() << "conventional total energy computed from the densities: " << setprecision(12) << energy_from_densities << endl;
      ExEnv::out0() << "conventional total energy for comparison: " << setprecision(12) << energy_conv << endl;
      if(fabs(energy_from_densities-energy_conv) > 10.0*value_acc) {
        throw ProgrammingError("PsiCI_PT2R12::compute -- the used densities don't reproduce the conventional energy.",__FILE__,__LINE__);
      }
    }

    double energy_correction_r12 = 0.0;
    // compute R12 correction only if given non-null correlation factor
    if(r12evalinfo_->r12tech()->corrfactor()->geminaldescriptor()->type()!=string("invalid")) {
      const RefSymmSCMatrix opdm_alpha = onepdm_refmo(Alpha);
      const RefSymmSCMatrix opdm_beta = onepdm_refmo(Beta);
      r12eval_->r12info()->set_opdm(opdm_alpha,opdm_beta);

      double energy_pt2r12[NSpinCases2];
      if(r12evalinfo_->ansatz()->projector()==LinearR12::Projector_1) {
        for(int i=0; i<NSpinCases2; i++) {
          SpinCase2 pairspin = static_cast<SpinCase2>(i);
          energy_pt2r12[i] = energy_PT2R12_projector1(pairspin);
          energy_correction_r12 +=  energy_pt2r12[i];
        }
      }
      else {
        for(int i=0; i<NSpinCases2; i++) {
          SpinCase2 pairspin = static_cast<SpinCase2>(i);
          energy_pt2r12[i] = energy_PT2R12_projector2(pairspin);
          energy_correction_r12 +=  energy_pt2r12[i];
        }
      }
    }

    ExEnv::out0() << "SCF energy: " << setprecision(12) << energy_scf << endl;
    ExEnv::out0() << "rasscf energy: " << setprecision(12) << energy_rasscf << endl;
    ExEnv::out0() << "conventional total energy: " << setprecision(12) << energy_conv << endl;
    ExEnv::out0() << "conventional correlation energy: " << setprecision(12) << (energy_conv-energy_scf) << endl;
    const double energy = exenv()->chkpt().rd_etot() + energy_correction_r12;
    ExEnv::out0() << "PT2R12 energy correction: " << setprecision(12) << energy_correction_r12 << endl;
    ExEnv::out0() << "PT2R12 corrected total energy: " << setprecision(12) << energy << endl;
    ExEnv::out0() << "PT2R12 corrected correlation energy: " << setprecision(12) << (energy-energy_scf) << endl;

    set_energy(energy);
  }

  RefSCMatrix PsiCI_PT2R12::overlap_with_valence_obwfn() {
    RefSCMatrix eigenvec_wfn = reference_mpqc_->eigenvectors();
    RefSCMatrix eigenvec_reference_wfn = valence_obwfn_->eigenvectors();
    Ref<GaussianBasisSet> basis_wfn = reference_mpqc_->basis();
    Ref<GaussianBasisSet> basis_reference_wfn = valence_obwfn_->basis();
    Ref<Integral> integral_wfn = reference_mpqc_->integral();
    Ref<Integral> integral_reference_wfn = valence_obwfn_->integral();
    Ref<PetiteList> plist_wfn = integral_wfn->petite_list();
    Ref<PetiteList> plist_reference_wfn = integral_reference_wfn->petite_list();
    RefSCMatrix eigenvec_wfn_ao = plist_wfn->evecs_to_AO_basis(eigenvec_wfn);
    RefSCMatrix eigenvec_valence_obwfn_ao = plist_reference_wfn->evecs_to_AO_basis(eigenvec_reference_wfn);

    Ref<OrbitalSpace> space_wfn = new OrbitalSpace("WFN","SCF MO space",eigenvec_wfn_ao,basis_wfn,integral_wfn);
    Ref<OrbitalSpace> space_reference_wfn = new OrbitalSpace("REFERENCE","SCF reference MO space",eigenvec_valence_obwfn_ao,basis_reference_wfn,integral_reference_wfn);

    RefSCMatrix overlap;
    R12IntEvalInfo::compute_overlap_ints(space_reference_wfn,space_wfn,overlap);

    return(overlap);
  }

  std::map<unsigned int, unsigned int> PsiCI_PT2R12::reference_index_map(const RefSCMatrix &valuematrix) {
    /// create list of overlap matrix elements of type 'Element'
    list<psici_pt2r12::Element> elem_list = psici_pt2r12::create_indexed_list_from_value_matrix(valuematrix);

//#define PRINT_DIMS
#ifdef PRINT_DIMS
    ExEnv::out0() << "reference wave function block dimensions:" << endl;
    psici_pt2r12::print_block_dims(valuematrix.rowdim());
    ExEnv::out0() << "wave function block dimensions:" << endl;
    psici_pt2r12::print_block_dims(valuematrix.coldim());
#endif

    /// sort the list of overlap matrix elements
    elem_list.sort(greater<psici_pt2r12::Element>());

    std::map<unsigned int,unsigned int> index_map = psici_pt2r12::reference_val_map(elem_list);

    return(index_map);
  }

  std::map<unsigned int,unsigned int> PsiCI_PT2R12::standard_index_map(const RefSCMatrix &valuematrix) {
    list<psici_pt2r12::Element> elem_list = psici_pt2r12::create_indexed_list_from_value_matrix(valuematrix);
    std::map<unsigned int, unsigned int> index_map = psici_pt2r12::reference_val_map(elem_list);

    return(index_map);
  }

  std::map<unsigned int, unsigned int> PsiCI_PT2R12::psi_index_map(const RefSCMatrix &valuematrix) {
    std::map<unsigned int, unsigned int> std_index_map = standard_index_map(valuematrix);
    std::map<unsigned int, unsigned int> ref_index_map = reference_index_map(valuematrix);
    std::map<unsigned int, unsigned int> new_index_map;
    for(std::map<unsigned int, unsigned int>::const_iterator elem1 = std_index_map.begin(); elem1 != std_index_map.end(); ++elem1) {
      for(std::map<unsigned int, unsigned int>::const_iterator elem2 = ref_index_map.begin(); elem2 != ref_index_map.end(); ++elem2) {
        if((elem1->first==elem2->first) && (elem1->second!=elem2->second)) {
          new_index_map.insert(pair<int,int>(elem1->second,elem2->second));
        }
      }
    }
    if(!new_index_map.empty()) {
      ExEnv::out0() << "MO indices to be exchanged in symmetry blocked (Psi3) order" << endl;\
      ExEnv::out0() << "to make the lowest orbitals correspond to the reference wave function:" << endl;
      psici_pt2r12::print_int_index_map(new_index_map);
    }

    return(new_index_map);
  }

  vector<unsigned int> PsiCI_PT2R12::map_to_valence_order(const RefSCMatrix &overlap) {
    const int nmo = reference_mpqc_->oso_dimension().n();
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

    ExEnv::out0() << "Symmetry blocked MO index order where the lowest MOs correspond" << endl;
    ExEnv::out0() << "to the respective orbitals of the reference wave function." << endl;
    ExEnv::out0() << "This vector can be used as it is for the Psi3 cscf reorder keyword." << endl;
    psici_pt2r12::print_vector_int(orbital_order);

    return(orbital_order);
  }

  vector<unsigned int> PsiCI_PT2R12::mo_symmetries_in_energetic_order() {
    Ref<OrbitalSpace> orbs_sb = r12evalinfo_->refinfo()->orbs_sb();
    RefDiagSCMatrix evals_sb = orbs_sb->evals();
    const int nmo = evals_sb.dim().n();
    Ref<SCBlockInfo> blocks_sb = evals_sb.dim()->blocks();

    std::vector<psici_pt2r12::Energy_SymmIndex> energy_symmindex_vec(nmo);
    int ind = 0;
    for(int i=0; i<blocks_sb->nblock(); i++) {
      const int blocksize = blocks_sb->subdim(i).n();
      for(int j=0; j<blocksize; j++) {
        energy_symmindex_vec[ind] = psici_pt2r12::Energy_SymmIndex(i,evals_sb.get_element(ind));

        ind += 1;
      }
    }

    sort(energy_symmindex_vec.begin(),energy_symmindex_vec.end());

    vector<unsigned int> symmindex(nmo);
    for(int i=0; i<nmo; i++) {
      symmindex[i] = energy_symmindex_vec[i].symmindex;
    }

    return(symmindex);
  }

  vector<unsigned int> PsiCI_PT2R12::mo_blocks() {
    Ref<SCBlockInfo> blocks_sb = r12evalinfo_->refinfo()->orbs_sb()->evals().dim()->blocks();
    int nirrep = blocks_sb->nblock();
    vector<unsigned int> mos(nirrep);
    mos.assign(nirrep,0);
    for(int i=0; i<nirrep; i++) {
      mos[i] = blocks_sb->subdim(i);
    }
    return(mos);
  }

  vector<unsigned int> PsiCI_PT2R12::frzc_blocks() {
    Ref<SCBlockInfo> blocks_sb = r12evalinfo_->refinfo()->orbs_sb()->evals().dim()->blocks();
    int nfzc = r12evalinfo_->refinfo()->nfzc();
    int nirrep = blocks_sb->nblock();
    vector<unsigned int> frzc_vec(nirrep);
    frzc_vec.assign(nirrep,0);
    vector<unsigned int> symmindex = mo_symmetries_in_energetic_order();
    for(int i=0; i<nfzc; i++) {
      frzc_vec[symmindex[i]] += 1;
    }
    return(frzc_vec);
  }

  vector<unsigned int> PsiCI_PT2R12::docc_blocks() {
    Ref<SCBlockInfo> blocks_sb = r12evalinfo_->refinfo()->docc_sb()->evals().dim()->blocks();
    int nirrep = blocks_sb->nblock();
    vector<unsigned int> docc(nirrep);
    docc.assign(nirrep,0);
    for(int i=0; i<nirrep; i++) {
      docc[i] = blocks_sb->subdim(i);
    }
    return(docc);
  }

  vector<unsigned int> PsiCI_PT2R12::docc_act_blocks() {
    vector<unsigned int> docc = docc_blocks();
    vector<unsigned int> frzc = frzc_blocks();
    int nirrep = docc.size();
    vector<unsigned int> docc_act(nirrep);
    for(int i=0; i<nirrep; i++) {
      docc_act[i] = docc[i] - frzc[i];
    }
    return(docc_act);
  }

  vector<unsigned int> PsiCI_PT2R12::socc_blocks() {
    Ref<SCBlockInfo> blocks_sb = r12evalinfo_->refinfo()->orbs_sb()->evals().dim()->blocks();
    int nirrep = blocks_sb->nblock();
    vector<unsigned int> socc(nirrep);
    socc.assign(nirrep,0);
    vector<unsigned int> symmindex = mo_symmetries_in_energetic_order();
    const int nmo = symmindex.size();
    double tol = 1.0e-6;
    for(int i=0; i<nmo; i++) {
      if(fabs(r12evalinfo_->refinfo()->ref()->occupation(i)-1.0)<tol) {
        socc[symmindex[i]] += 1;
      }
    }
    return(socc);
  }

  vector<unsigned int> PsiCI_PT2R12::uocc_act_blocks() {
    vector<unsigned int> uocc = uocc_blocks();
    vector<unsigned int> frzv = frzv_blocks();
    int nirrep = uocc.size();
    vector<unsigned int> uocc_act(nirrep);
    for(int i=0; i<nirrep; i++) {
      uocc_act[i] = uocc[i] - frzv[i];
    }
    return(uocc_act);
  }

  vector<unsigned int> PsiCI_PT2R12::uocc_blocks() {
    Ref<SCBlockInfo> blocks_sb = r12evalinfo_->refinfo()->uocc_sb()->evals().dim()->blocks();
    int nirrep = blocks_sb->nblock();
    int nirrep_docc = r12evalinfo_->refinfo()->docc_sb()->evals().dim()->blocks()->nblock();
    vector<unsigned int> uocc(nirrep_docc);
    if(nirrep!=nirrep_docc){
      uocc.assign(uocc.size(),0.0);
      return(uocc);
    }
    for(int i=0; i<nirrep; i++) {
      uocc[i] = blocks_sb->subdim(i);
    }
    return(uocc);
  }

  vector<unsigned int> PsiCI_PT2R12::frzv_blocks() {
    Ref<SCBlockInfo> blocks_sb = r12evalinfo_->refinfo()->orbs_sb()->evals().dim()->blocks();
    int nfzv = r12evalinfo_->refinfo()->nfzv();
    int nirrep = blocks_sb->nblock();
    int nmo = blocks_sb->nelem();
    vector<unsigned int> frzv_vec(nirrep);
    frzv_vec.assign(nirrep,0);
    vector<unsigned int> symmindex = mo_symmetries_in_energetic_order();
    for(int i=nmo-1,j=0; (j<nfzv) && (i>=0); j++,i--) {
      frzv_vec[symmindex[i]] += 1;
    }
    return(frzv_vec);
  }

  void PsiCI_PT2R12::print_all_blocks(ostream &o) {
    vector<unsigned int> mos = mo_blocks();
    psici_pt2r12::print_blocks("mos",mos,o);
    vector<unsigned int> frzc = frzc_blocks();
    psici_pt2r12::print_blocks("frzc",frzc,o);
    vector<unsigned int> docc = docc_blocks();
    psici_pt2r12::print_blocks("docc",docc,o);
    vector<unsigned int> docc_act = docc_act_blocks();
    psici_pt2r12::print_blocks("docc_act",docc_act,o);
    vector<unsigned int> socc = socc_blocks();
    psici_pt2r12::print_blocks("socc",socc,o);
    vector<unsigned int> uocc_act = uocc_act_blocks();
    psici_pt2r12::print_blocks("uocc_act",uocc_act,o);
    vector<unsigned int> uocc = uocc_blocks();
    psici_pt2r12::print_blocks("uocc",uocc,o);
    vector<unsigned int> frzv = frzv_blocks();
    psici_pt2r12::print_blocks("frzv",frzv,o);
  }

}
