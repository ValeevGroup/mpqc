//
// ref.h
//
// Copyright (C) 2005 Edward Valeev
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

#ifndef _chemistry_qc_mbptr12_ref_h
#define _chemistry_qc_mbptr12_ref_h

#include <chemistry/qc/wfn/wfn.h>
#include <chemistry/qc/scf/scf.h>
#include <chemistry/qc/lcao/wfnworld.h>
#include <chemistry/qc/wfn/orbitalspace.h>
#include <chemistry/qc/wfn/spin.h>

namespace sc {

  class Wavefunction;

  /// PopulatedOrbitalSpace is an OrbitalSpace populated with a density.
  /// It holds OrbitalSpaces representing subsets of the OrbitalSpace,
  /// for example, corresponding to various occupancies.
  class PopulatedOrbitalSpace : virtual public SavableState {
    public:
    /**
     * @param oreg an OrbitalSpaceRegistry object that will know of the computed spaces
     * @param spin spin-case that will be used to compute labels of OrbitalSpace objects
     * @param bs basis set
     * @param integral Integral factory used to support coefficients
     * @param coefs coefficients of orbitals expanded in basis (AO by MO matrix)
     * @param occs occupation vector
     * @param active the mask used to freeze orbitals (active[i] == false means i will be frozen)
     * @param energies orbital energies.
     * @param eorder_increasing if true, energy-ordered spaces will order orbitals in the order of increasing energies
     * @param vbs OrbitalSpace that represents the unoccupied orbitals.
     *            The default is 0, which means to use empty orbitals from coefs.
     * @param fbrun the FockBuildRuntime object used to compute Fock matrices. if vbs != 0, fbrun must be specified.
     */
      PopulatedOrbitalSpace(const Ref<OrbitalSpaceRegistry>& oreg,
                            SpinCase1 spin, const Ref<GaussianBasisSet>& bs,
                            const Ref<Integral>& integral,
                            const RefSCMatrix& coefs,
                            const std::vector<double>& occs,
                            const std::vector<bool>& active,
                            const RefDiagSCMatrix& energies,
                            bool eorder_increasing = true,
                            Ref<OrbitalSpace> vbs = 0,
                            Ref<FockBuildRuntime> fbrun = 0,
                            const std::vector<double>& rasscf_occs=std::vector<double>(1,-1.0)
                           ); //if rasscf_occs is specified, uses rasscf_occs instead of occs to construct occ_act_mask,
                              //eventually to construct occ_act_sb_ differently->force GG/gg space to be ras1+ras2 orbs
      PopulatedOrbitalSpace(const bool doscreen, const double occ_thres, RefSymmSCMatrix OBS_mo_ordm, const Ref<OrbitalSpaceRegistry>& oreg,
                                  SpinCase1 spin, const Ref<GaussianBasisSet>& bs,
                                  const Ref<Integral>& integral,
                                  RefSCMatrix& coefs,
                                  const std::vector<double>& occs,
                                  std::vector<bool>& active,
                                  const RefDiagSCMatrix& energies,
                                  bool eorder_increasing = true,
                                  Ref<OrbitalSpace> vbs = 0,
                                  Ref<FockBuildRuntime> fbrun = 0); // this is to construct screened orb spaces
      PopulatedOrbitalSpace(StateIn& si);
      ~PopulatedOrbitalSpace();
      void save_data_state(StateOut& so);

      /// an orbital is occupied if its occupancy is greater than this
      static double zero_occupancy() { return 1e-12; }

      const Ref<OrbitalSpaceRegistry>& orbital_registry() const { return oreg_; }

      const Ref<OrbitalSpace>& orbs_sb() const { return orbs_sb_; }
      const Ref<OrbitalSpace>& orbs() const { return orbs_; }
      const Ref<OrbitalSpace>& occ_sb() const { return occ_sb_; }
      const Ref<OrbitalSpace>& occ_act_sb() const { return occ_act_sb_; }
      const Ref<OrbitalSpace>& unscreen_occ_act_sb() const { return unscreen_occ_act_sb_; }
      const Ref<OrbitalSpace>& occ() const { return occ_; }
      const Ref<OrbitalSpace>& occ_act() const { return occ_act_; }
      const Ref<OrbitalSpace>& uocc_sb() const { return uocc_sb_; }
      const Ref<OrbitalSpace>& conv_uocc_sb() const { return conv_uocc_sb_; }
      const Ref<OrbitalSpace>& conv_occ_sb() const { return conv_occ_sb_; }
      const Ref<OrbitalSpace>& uocc_act_sb() const { return uocc_act_sb_; }
      const Ref<OrbitalSpace>& uocc() const { return uocc_; }
      const Ref<OrbitalSpace>& uocc_act() const { return uocc_act_; }

    private:
      Ref<OrbitalSpaceRegistry> oreg_;
      Ref<OrbitalSpace> orbs_sb_;
      Ref<OrbitalSpace> orbs_;
      Ref<OrbitalSpace> occ_sb_;
      Ref<OrbitalSpace> occ_act_sb_;
      Ref<OrbitalSpace> unscreen_occ_act_sb_; // we keep the unscreened (but rotated) orbitals to use when transformation RDM
      Ref<OrbitalSpace> occ_;
      Ref<OrbitalSpace> occ_act_;
      Ref<OrbitalSpace> uocc_sb_;
      Ref<OrbitalSpace> uocc_act_sb_;
      Ref<OrbitalSpace> uocc_;
      Ref<OrbitalSpace> uocc_act_;
      Ref<OrbitalSpace> conv_uocc_;//this denotes the virtual orb space in the conventional sense,e.g. in MRCI, they are orbs other than
                                   // core and active occupied ones
      Ref<OrbitalSpace> conv_uocc_sb_;
      Ref<OrbitalSpace> conv_occ_sb_;

      /// purges the spaces from the registry
      void purge();
  };

  /**
     RefWavefunction represents the reference wave function used in correlated calculations.

     Single-determinantal and multi-determinantal Wavefunction can be used as a reference.
     However since Wavefunction does not have proper constructors, it's not implemented
     as a proper Adapter to Wavefunction and thus implements many member functions of Wavefunction.
     The main content is a set of OrbitalSpace objects.
  */
  class RefWavefunction : virtual public SavableState {
    protected:

      /// return true to override density fitting settings of RefWavefunction objects to
      /// those of WavefunctionWorld in which they live. This should only be used for testing
      /// as it may produce Fock matrices that will differ from those of the original reference.
      static bool use_world_dfinfo() {
        return false;
      }

      RefWavefunction(StateIn&);
      /// @param world WavefunctionWorld to which this object belongs
      /// @param basis The basis set supporting the reference wave function
      /// @param integral The integral factory used to compute the reference wavefunction
      RefWavefunction(const Ref<WavefunctionWorld>& world,
                      const Ref<GaussianBasisSet>& basis,
                      const Ref<Integral>& integral);

    public:
    ~RefWavefunction();
    void save_data_state(StateOut&);

    /// obsoletes this object
    /// @sa R12WavefunctionWorld::obsolete()
    virtual void obsolete();

    const Ref<WavefunctionWorld>& world() const { return world_; }
    const Ref<GaussianBasisSet>& basis() const { return basis_; }
    const Ref<Integral>& integral() const { return integral_; }
    /// returns the basis supporting unoccupied orbitals. The defauls is same as returned by basis().
    virtual const Ref<GaussianBasisSet>& uocc_basis() const { return basis(); }

    /// if true, force alpha and beta 1-rdm the same for spin-free algorithms
    void set_spinfree(bool TrueOrFalse);
    /// @sa MolecularEnergy::energy()
    virtual double energy() =0;
    /// Set the accuracy to which the value is to be computed. @sa Function::set_desired_value_accuracy()
    void set_desired_value_accuracy(double);
    /// Return the accuracy with which the value has been computed. @sa Function::actual_value_accuracy()
    virtual double actual_value_accuracy () const =0;
    /// @sa Return the accuracy with which the value is to be computed. @sa Function::desired_value_accuracy()
    virtual double desired_value_accuracy() const =0;
    /// @sa Return true if the accuracy was set to default. @sa Function::desired_value_accuracy_set_to_default()
    virtual bool desired_value_accuracy_set_to_default() const;
    /// @sa Wavefunction::nelectron()
    virtual int nelectron() const =0;
    /// @sa Wavefunction::spin_polarized()
    virtual bool spin_polarized() const =0;
    /// @sa Wavefunction::dk()
    virtual int dk() const =0;
    /// @sa Wavefunction::momentum_basis()
    virtual Ref<GaussianBasisSet> momentum_basis() const =0;
    /** Returns the SO core Hamiltonian in the given basis and momentum
        basis.  The momentum basis is not needed if no Douglas-Kroll
        correction is being performed.
    @sa Wavefunction::core_hamiltonian_for_basis()  */
    virtual RefSymmSCMatrix core_hamiltonian_for_basis(const Ref<GaussianBasisSet> &basis,
                                                       const Ref<GaussianBasisSet> &p_basis) =0;
    /// return the AO basis density
    virtual RefSymmSCMatrix ordm(SpinCase1 spin) const =0;
    /// return the MO basis density (MOs are given by orbs_sb())
    virtual RefSymmSCMatrix ordm_orbs_sb(SpinCase1 spin) const;
    virtual RefSymmSCMatrix ordm_occ_sb(SpinCase1 spin) const;

    /// is this a single-determinantal reference?
    virtual bool sdref() const =0;

    /// Returns the space of symmetry-blocked orthogonal SOs (spans the entire space of the basis)
    const Ref<OrbitalSpace>& oso_space() const;
    /// Return the space of symmetry-blocked MOs of the given spin
    const Ref<OrbitalSpace>& orbs_sb(SpinCase1 spin = AnySpinCase1) const;
    /// Return the space of energy-sorted MOs of the given spin
    const Ref<OrbitalSpace>& orbs(SpinCase1 spin = AnySpinCase1) const;
    /// Return the space of symmery-blocked occupied MOs of the given spin
    const Ref<OrbitalSpace>& occ_sb(SpinCase1 spin = AnySpinCase1) const;
    /// Return the space of symmery-blocked active occupied MOs of the given spin
    const Ref<OrbitalSpace>& occ_act_sb(SpinCase1 spin = AnySpinCase1) const;
    /// Return the space of symmery-blocked active occupied MOs of the given spin, (orbs rotated but not screened)
    const Ref<OrbitalSpace>& unscreen_occ_act_sb(SpinCase1 spin = AnySpinCase1) const;
    /// Return the space of occupied MOs of the given spin
    const Ref<OrbitalSpace>& occ(SpinCase1 spin = AnySpinCase1) const;
    /// Return the space of active occupied MOs of the given spin
    const Ref<OrbitalSpace>& occ_act(SpinCase1 spin = AnySpinCase1) const;
    /// Return the space of symmetry-blocked unoccupied (virtual) MOs of the given spin
    const Ref<OrbitalSpace>& uocc_sb(SpinCase1 spin = AnySpinCase1) const;
    const Ref<OrbitalSpace>& conv_uocc_sb(SpinCase1 spin = AnySpinCase1) const;
    const Ref<OrbitalSpace>& conv_occ_sb(SpinCase1 spin = AnySpinCase1) const;
    /// Return the space of symmetry-blocked active unoccupied (virtual) MOs of the given spin
    const Ref<OrbitalSpace>& uocc_act_sb(SpinCase1 spin = AnySpinCase1) const;
    /// Return the space of unoccupied (virtual) MOs of the given spin
    const Ref<OrbitalSpace>& uocc(SpinCase1 spin = AnySpinCase1) const;
    /// Return the space of active unoccupied (virtual) MOs of the given spin
    const Ref<OrbitalSpace>& uocc_act(SpinCase1 spin = AnySpinCase1) const;

    // 'original' orbs_sb: without transformation to the screened orbital space
    const Ref<OrbitalSpace>& orig_orbs_sb(SpinCase1 spin = AnySpinCase1) const;

    double occ_thres() {return occ_thres_;}
    void set_occ_thres(double vv) {occ_thres_ = vv;}
    bool force_rasscf(){return force_correlate_rasscf_;}
    void set_force_correlate_rasscf(const bool force) {force_correlate_rasscf_ = force;}
    bool do_screen() {return do_screen_;}
    void set_do_screen(bool screenornot) {do_screen_ = screenornot;}
    Ref<PopulatedOrbitalSpace> & get_poporbspace(SpinCase1 spin = Alpha) {return spinspaces_[spin];}
    Ref<PopulatedOrbitalSpace>& get_screened_poporbspace(SpinCase1 spin = Alpha) {return screened_spinspaces_[spin];}
    /// which DensityFittingRuntime used to compute this reference wave function
    virtual Ref<DensityFittingInfo> dfinfo() const =0;

    protected:
    /// For spin-free algorithms, if this is true, we would set both alpha/beta 1-rdm to the average of them; defaults to false
    bool force_average_AB_rdm1_;

    private:
    Ref<WavefunctionWorld> world_;   // who owns this?
    Ref<GaussianBasisSet> basis_;
    Ref<Integral> integral_;
    bool omit_uocc_;
    bool force_correlate_rasscf_;
    double occ_thres_; // this parameter is in fact assigned to the value of correlate_min_occ_ of R12WavefunctionWorld, where is a logical
                        // place to define such a variable, since it more belong to the R12 world. However, we need it to control RefWavefunction too.
    bool do_screen_;


    /// used to implement set_desired_value_accuracy()
    virtual void _set_desired_value_accuracy(double eps) =0;

    protected:
    /// initializes the object
    void init() const;
    /// calling this will cause the object to be re-initialized next time it is used
    virtual void reset();
    mutable Ref<PopulatedOrbitalSpace> spinspaces_[NSpinCases1];
    mutable Ref<PopulatedOrbitalSpace> screened_spinspaces_[NSpinCases1]; // this is used to construct the orbital spaces which are screened based on occ number using natural orbitals.
    bool screened_space_init_ed_; // to avoid potential complexity and redundancy, we use this parameter to monitor whether "screened_spincases_" is initialized or not,
                                  // so that we only construct it once.
    bool orig_space_init_ed_; // created to potentially monitor the initialization state of spincases_; I find the space is initialized again and again now by calling init_spaces(), confusing.

    /// initialize OrbitalSpace objects
    virtual void init_spaces() = 0;
  };

  /// RefWavefunction specialization for a single-determinant wave function
  class SD_RefWavefunction : public RefWavefunction {
    public:
      /// construct from a OneBodyWavefunction object
      /// @param world The WavefunctionWorld in which this objects lives.
      /// @param obwfn The OneBodyWavefunction object that specifies the orbitals
      /// @param spin_restricted If false and obwfn is an open-shell spin-restricted OneBodyWavefunction,
      ///        will use semicanonical orbitals. The value of this parameter will be ignored for closed-shell
      ///        and spin-unrestricted open-shell references.
      /// @param nfzc The number of lowest-energy occupied orbitals to be kept inactive
      /// @param nfzv The number of highest-energy unoccupied orbitals to be kept inactive
      /// @param vir_space The space describing the unoccupied orbitals. Default is 0, which
      ///        means use unoccupied orbitals from obwfn.
      /// @param occ_orbs The method used to provide the occupied orbitals. This can
      /// be "pipek-mezey or "canonical".
      ///
      /// N.B. This will feed the FockBuildRuntime in world with the density matrices from obwfn!
      SD_RefWavefunction(const Ref<WavefunctionWorld>& world,
                            const Ref<OneBodyWavefunction>& obwfn,
                            bool spin_restricted = true,
                            unsigned int nfzc = 0,
                            unsigned int nfzv = 0,
                            Ref<OrbitalSpace> vir_space = 0,
                            std::string occ_orbs = std::string("canonical"));
      SD_RefWavefunction(StateIn&);
      ~SD_RefWavefunction();
      void save_data_state(StateOut&);

      void obsolete();

      bool sdref() const { return true; }
      const Ref<OneBodyWavefunction>& obwfn() const { return obwfn_; }
      const Ref<OrbitalSpace>& vir_space() const { return vir_space_; }
      const Ref<GaussianBasisSet>& uocc_basis() const {
        if (vir_space_.nonnull()) return vir_space_->basis();
        else return this->basis();
      }

      double energy() { return obwfn()->energy(); }
      double actual_value_accuracy () const { return obwfn()->actual_value_accuracy(); }
      double desired_value_accuracy() const { return obwfn()->desired_value_accuracy(); }
      bool desired_value_accuracy_set_to_default() const { return obwfn()->desired_value_accuracy_set_to_default(); }
      int nelectron() const { return obwfn()->nelectron(); }
      bool spin_polarized() const { return obwfn_->spin_polarized(); }
      bool spin_restricted() const { return spin_restricted_; }
      int dk() const { return obwfn()->dk(); }
      Ref<GaussianBasisSet> momentum_basis() const { return obwfn()->momentum_basis(); }
      RefSymmSCMatrix core_hamiltonian_for_basis(const Ref<GaussianBasisSet> &basis,
                                                 const Ref<GaussianBasisSet> &p_basis);
      unsigned int nfzc() const { return nfzc_; }
      unsigned int nfzv() const { return nfzv_; }
      RefSymmSCMatrix ordm(SpinCase1 spin) const;
      Ref<DensityFittingInfo> dfinfo() const;
      const std::string& occ_orbitals() const { return occ_orbitals_; }
    private:
      Ref<OneBodyWavefunction> obwfn_;
      Ref<OrbitalSpace> vir_space_;
      bool spin_restricted_;
      std::string occ_orbitals_;
      unsigned int nfzc_;
      unsigned int nfzv_;
      void init_spaces();
      void init_spaces_restricted();
      void init_spaces_unrestricted();
      void _set_desired_value_accuracy(double eps) { obwfn_->set_desired_value_accuracy(eps); }
  };

  /// RefWavefunction specialization for a general wave function specified by its orbitals and rank-1 reduced density matrices
  class Extern_RefWavefunction : public RefWavefunction {
    public:
      /// Constructs Extern_RefWavefunction using the MO-basis 1-RDMs + MO coefficients (same for alpha and beta spincase)
      /// @param world The WavefunctionWorld in which this objects lives.
      /// @param basis The basis set
      /// @param integral The integral object that determines the ordering of basis functions in shells
      /// @param orbs  The MO coefficient matrix
      /// @param symm  Irreps of MOs
      /// @param alpha_1rdm The alpha-spin density matrix in MO basis
      /// @param beta_1rdm The beta-spin density matrix in MO basis (assuming if alpha_1rdm and beta_1rdm point to the SAME object
      ///   then if alpha and beta densities are identical.
      /// @param nocc orbitals [0,nocc) will be occupied
      /// @param nfzc orbitals [0,nfzc) will be inactive
      /// @param nfzv orbitals [nmo-nfzv,nmo) will be inactive
      /// @param omit_uocc If true, omit all unoccupied orbitals (i.e. make the unoccupied space empty). N.B. This is
      ///                      not the same as "freezing" the unoccupieds.
      Extern_RefWavefunction(const Ref<WavefunctionWorld>& world,
                  const Ref<GaussianBasisSet>& basis,
                  const Ref<Integral>& integral,
                  const RefSCMatrix& orbs,
                  const std::vector<unsigned int>& orbsymm,
                  const RefSymmSCMatrix& alpha_1rdm,
                  const RefSymmSCMatrix& beta_1rdm,
                  unsigned int nocc,
                  unsigned int nfzc = 0,
                  unsigned int nfzv = 0,
                  bool omit_uocc = false);
      Extern_RefWavefunction(const Ref<WavefunctionWorld>& world,
                  const Ref<GaussianBasisSet>& basis,
                  const Ref<Integral>& integral,
                  const RefSCMatrix& orbs,
                  const std::vector<unsigned int>& orbsymm,
                  const RefSymmSCMatrix& alpha_1rdm,
                  const RefSymmSCMatrix& beta_1rdm,
                  std::vector<unsigned int> mopi,
                  std::vector<unsigned int> occpi,
                  std::vector<unsigned int> corrpi,
                  std::vector<unsigned int> fzcpi,
                  std::vector<unsigned int> fzvpi,
                  bool force_correlate_rasscf = false,
                  bool omit_uocc = false);
      Extern_RefWavefunction(StateIn&);
      virtual ~Extern_RefWavefunction();
      void save_data_state(StateOut&);
      RefSymmSCMatrix ordm(SpinCase1 spin) const { return rdm_[spin]; }

      void obsolete() {
//        throw FeatureNotImplemented("cannot obsolete Extern_R12RefWavefunction",
//                                                    __FILE__, __LINE__);
      }

      bool sdref() const;
      double energy() { return 0.0; }
      double actual_value_accuracy () const { return DBL_EPSILON; }
      double desired_value_accuracy() const { return DBL_EPSILON; }
      int nelectron() const { return nelectron_; }
      bool spin_polarized() const { return rdm_[Alpha] != rdm_[Beta]; }
      bool spin_restricted() const { return true; }
      /// reimplements RefWavefunction::dk(). Currently only nonrelativistic references are supported.
      int dk() const { return 0; }
      Ref<GaussianBasisSet> momentum_basis() const { return 0; }
      RefSymmSCMatrix core_hamiltonian_for_basis(const Ref<GaussianBasisSet> &basis,
                                                 const Ref<GaussianBasisSet> &p_basis);
      unsigned int nfzc() const { return nfzc_; }
      unsigned int nfzv() const { return nfzv_; }
      bool omit_uocc() const { return omit_uocc_; }
      bool ordm_idempotent() const { return ordm_idempotent_; }
    private:
      RefSymmSCMatrix rdm_[NSpinCases1];
      unsigned int nfzc_;
      unsigned int nfzv_;
      unsigned int nelectron_;
      bool omit_uocc_;
      bool ordm_idempotent_;


      void pre_init(std::vector<unsigned int> mopi,
                    std::vector<unsigned int> occpi,
                    std::vector<unsigned int> corrpi,
                    std::vector<unsigned int> fzcpi,
                    std::vector<unsigned int> fzvpi,
                    const RefSCMatrix& orbs,
                const std::vector<unsigned int>& orbsym);
      void init_spaces() {throw sc::ProgrammingError("For Extern_RefWavefunction, spaces must be init-ed in constructor");}
      void init_spaces(unsigned int nocc, const RefSCMatrix& orbs,
                       const std::vector<unsigned int>& orbsym);
      void init_spaces(std::vector<unsigned int> mopi,
                       std::vector<unsigned int> occpi,
                       std::vector<unsigned int> corrpi,
                       std::vector<unsigned int> fzcpi,
                       std::vector<unsigned int> fzvpi,
                       const RefSCMatrix& orbs,
                   const std::vector<unsigned int>& orbsym);
      void _set_desired_value_accuracy(double eps) {
        // do nothing
      }

      Ref<DensityFittingInfo> dfinfo() const;
  };

  /// This factory produces the RefWavefunction that corresponds to the type of ref object
  struct RefWavefunctionFactory {
      static Ref<RefWavefunction>
        make(const Ref<WavefunctionWorld>& world,
             const Ref<Wavefunction>& ref,
             bool spin_restricted = true,
             unsigned int nfzc = 0,
             unsigned int nfzv = 0,
             Ref<OrbitalSpace> vir_space = 0);
  };

  /// canonicalize A
  Ref<OrbitalSpace>
  compute_canonvir_space(const Ref<FockBuildRuntime>& fb_rtime,
                         const Ref<OrbitalSpace>& A,
                         SpinCase1 spin);
};

#endif

