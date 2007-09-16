#ifdef __GNUC__
#pragma interface
#endif

#ifndef _chemistry_qc_psi_psiwfn_h
#define _chemistry_qc_psi_psiwfn_h

#include <chemistry/qc/wfn/wfn.h>
#include <chemistry/qc/psi/psiexenv.h>
#include <chemistry/qc/mbptr12/mbptr12.h>
#include <chemistry/qc/mbptr12/spin.h>
#include <chemistry/qc/mbptr12/moindexspace.h>

namespace sc {
  
  ///////////////////////////////////////////////////////////////////
  /** PsiWavefunction is an abstract base for all Psi wave functions.
   Its KeyVal constructor is invoked by all KeyVal constructors of
   concrete implementations of PsiWavefunction.
   */

  class PsiWavefunction : public Wavefunction {
      
      Ref<PsiExEnv> exenv_;
      /// All Psi wave functions can at least compute the energy
      int value_implemented() const {
        return 1;
      }
      
    protected:
      int nirrep_;
      std::vector<int> docc_;
      std::vector<int> socc_;
      int multp_;
      int charge_;
      char *memory_;
      /// Prepares a complete Psi input file. The input file is assumed to have been opened.
      virtual void write_input(int conv) =0;

      std::vector<int> read_occ(const Ref<KeyVal> &keyval, const char *name,
                                int nirrep);

      // return the debug level
      int debug() const;
      
    public:
      /** The KeyVal constructor.

       <dl>

       <dt><tt>psienv</tt><dd> Specifies a PsiExEnv object.  There
       is no default.

       <dt><tt>memory</tt><dd> This integer specifies the amount of memory
       (in bytes) for Psi to use. The default is 2000000.

       <dt><tt>debug</tt><dd> This integer can be used to produce output
       for debugging.  The default is 0.

       </dl> */
      PsiWavefunction(const Ref<KeyVal>&);
      PsiWavefunction(StateIn&);
      ~PsiWavefunction();

      void save_data_state(StateOut&);

      /** Writes out Psi input file entries specific to this PsiWavefunction.
       The input file is assumed to have been opened. */
      virtual void write_basic_input(int conv);
      void compute();
      void print(std::ostream&o=ExEnv::out0()) const;
      RefSymmSCMatrix density();
      int nelectron();

      /// Return an associated PsiExEnv object
      Ref<PsiExEnv> exenv() const {
        return exenv_;
      }
      ;
      /// Return an associated PsiInput object
      Ref<PsiInput> get_psi_input() const {
        return exenv_->get_psi_input();
      }
      ;
  };
  
  ///////////////////////////////////////////////////////////////////
  /// PsiSCF is an abstract base for all Psi SCF wave functions

  class PsiSCF : public PsiWavefunction {
      RefDiagSCMatrix evals_;
      RefSCMatrix coefs_;
    public:
      PsiSCF(const Ref<KeyVal>&);
      PsiSCF(StateIn&);
      ~PsiSCF();
      void save_data_state(StateOut&);

      enum RefType {rhf, hsoshf, uhf};
      /// Returns the PsiSCF::RefType of this particular Psi SCF wave function
      virtual PsiSCF::RefType reftype() const =0;
      /// Returns the eigenvalues matrix
      virtual const RefDiagSCMatrix& evals();
      /// Returns the coefficient matrix
      virtual const RefSCMatrix& coefs();

      /// number of MOs
      unsigned int nmo();
      /// number of occupied MOs of spin
      unsigned int nocc(SpinCase1 spin);
  };
  
  ///////////////////////////////////////////////////////////////////
  /// PsiCLHF is a concrete implementation of Psi RHF wave function

  class PsiCLHF : public PsiSCF {
    protected:
      void write_input(int conv);
    public:
      PsiCLHF(const Ref<KeyVal>&);
      PsiCLHF(StateIn&);
      ~PsiCLHF();

      void write_basic_input(int conv);
      int spin_polarized() {
        return 0;
      }
      ;
      int gradient_implemented() const {
        return 1;
      }
      ;
      PsiSCF::RefType reftype() const {
        return rhf;
      }
      ;
  };
  
  ///////////////////////////////////////////////////////////////////
  /// PsiHSOSHF is a concrete implementation of Psi ROHF wave function

  class PsiHSOSHF : public PsiSCF {
    protected:
      void write_input(int conv);
    public:
      PsiHSOSHF(const Ref<KeyVal>&);
      PsiHSOSHF(StateIn&);
      ~PsiHSOSHF();

      void write_basic_input(int conv);
      int spin_polarized() {
        return 0;
      }
      ;
      int gradient_implemented() const {
        return 1;
      }
      ;
      PsiSCF::RefType reftype() const {
        return hsoshf;
      }
      ;
  };
  
  ///////////////////////////////////////////////////////////////////
  /// PsiUHF is a concrete implementation of Psi UHF wave function

  class PsiUHF : public PsiSCF {
    protected:
      void write_input(int conv);
    public:
      PsiUHF(const Ref<KeyVal>&);
      PsiUHF(StateIn&);
      ~PsiUHF();

      void write_basic_input(int conv);
      int spin_polarized() {
        return 1;
      }
      ;
      int gradient_implemented() const {
        return 1;
      }
      ;
      PsiSCF::RefType reftype() const {
        return uhf;
      }
      ;
  };
  
  ///////////////////////////////////////////////////////////////////
  /// PsiCorrWavefunction is a Psi correlated wave function

  class PsiCorrWavefunction : public PsiWavefunction {
    protected:
      Ref<PsiSCF> reference_;
      Ref<MOIndexSpace> occ_act_sb_;
      Ref<MOIndexSpace> vir_act_sb_;
      std::vector<int> frozen_docc_;
      std::vector<int> frozen_uocc_;
      void write_input(int conv);
    public:
      PsiCorrWavefunction(const Ref<KeyVal>&);
      PsiCorrWavefunction(StateIn&);
      ~PsiCorrWavefunction();
      void save_data_state(StateOut&);
      int spin_polarized() {
        return reference_->spin_polarized();
      }
      ;
      /// symmetry-blocked space of active occupied orbitals from Psi3
      const Ref<MOIndexSpace>& occ_act_sb(SpinCase1);
      /// symmetry-blocked space of active virtual orbitals from Psi3
      const Ref<MOIndexSpace>& vir_act_sb(SpinCase1);
  };
  
  ///////////////////////////////////////////////////////////////////
  /// PsiCC is a Psi coupled cluster wave function

  class PsiCC : public PsiCorrWavefunction {
      std::vector<RefSCMatrix> T_;
      RefSCMatrix Tau2_;
      std::vector<RefSCMatrix> Lambda_;

    protected:
      // set to true if want to run only if Psi3 and MPQC orbitals match exactly up to a phase
      static const bool use_sparsemap_only_ = false;

      /// set to true to test whether T2 transform from Psi3 to MPQC orbitals works correctly (causes Psi3 to do MP1 instead of CCSD)
      static bool test_t2_phases_;
      static void do_test_t2_phases() {
        test_t2_phases_ = true;
      }
      
      /// read in T1-like quantity using DPD label L
      RefSCMatrix T1(const std::string& L);
      /// read in T2-like quantity using DPD label L
      RefSCMatrix T2(const std::string& L);

      /// transform T1 to the new basis using sparse maps
      RefSCMatrix
          transform_T1(
                       const SparseMOIndexMap& occ_act_map,
                       const SparseMOIndexMap& vir_act_map,
                       const RefSCMatrix& T1,
                       const Ref<SCMatrixKit>& kit = SCMatrixKit::default_matrixkit()) const;
      /// transform T2 to the new basis using sparse maps
      RefSCMatrix
          transform_T2(
                       const SparseMOIndexMap& occ1_act_map,
                       const SparseMOIndexMap& occ2_act_map,
                       const SparseMOIndexMap& vir1_act_map,
                       const SparseMOIndexMap& vir2_act_map,
                       const RefSCMatrix& T2,
                       const Ref<SCMatrixKit>& kit = SCMatrixKit::default_matrixkit()) const;
      /// transform T1 to the new basis using dense transforms
      RefSCMatrix
          transform_T1(
                       const RefSCMatrix& occ_act_tform,
                       const RefSCMatrix& vir_act_tform,
                       const RefSCMatrix& T1,
                       const Ref<SCMatrixKit>& kit = SCMatrixKit::default_matrixkit()) const;
      /// transform T2 to the new basis using dense transforms
      RefSCMatrix
          transform_T2(
                       const RefSCMatrix& occ1_act_tform,
                       const RefSCMatrix& occ2_act_tform,
                       const RefSCMatrix& vir1_act_tform,
                       const RefSCMatrix& vir2_act_tform,
                       const RefSCMatrix& T2,
                       const Ref<SCMatrixKit>& kit = SCMatrixKit::default_matrixkit()) const;
      /// compare T2 and T2_ref (check that elements < zero are in the same place and elements > soft_zero have the same sign)
      void compare_T2(const RefSCMatrix& T2, const RefSCMatrix& T2_ref,
                      unsigned int no1, unsigned int no2, unsigned int nv1,
                      unsigned int nv2, double zero = 1e-8) const;

    public:
      PsiCC(const Ref<KeyVal>&);
      PsiCC(StateIn&);
      ~PsiCC();
      void save_data_state(StateOut&);

      /// return T amplitudes of rank i (i >= 1). The amplitudes are expressed in terms of Psi3 orbitals (symmetry-blocked).
      virtual const RefSCMatrix& T(unsigned int i);
      /// return Tau2 amplitudes. The amplitudes are expressed in terms of Psi3 orbitals (symmetry-blocked).
      virtual const RefSCMatrix& Tau2();
      /// return Lambda amplitudes of rank i (i >= 1)
      virtual const RefSCMatrix& Lambda(unsigned int i);
  };
  
  ///////////////////////////////////////////////////////////////////
  /// PsiCCSD is a concrete implementation of Psi CCSD wave function

  class PsiCCSD : public PsiCC {
    protected:
      void write_input(int conv);
    public:
      PsiCCSD(const Ref<KeyVal>&);
      PsiCCSD(StateIn&);
      ~PsiCCSD();
      void save_data_state(StateOut&);
      int gradient_implemented() const;
  };
  
  ///////////////////////////////////////////////////////////////////
  /// PsiCCSD_T is a concrete implementation of Psi CCSD(T) wave function

  class PsiCCSD_T : public PsiCC {
    protected:
      void write_input(int conv);
    public:
      PsiCCSD_T(const Ref<KeyVal>&);
      PsiCCSD_T(StateIn&);
      ~PsiCCSD_T();

      void save_data_state(StateOut&);
      int gradient_implemented() const;
  };
  
  ///////////////////////////////////////////////////////////////////
  /// PsiCCSD-PT2R12 is a concrete implementation of CCSD-PT2R12 wave function

  class MBPT2_R12;
  
  class PsiCCSD_PT2R12 : public PsiCC {
      /// set to true to test the code against MP2-R12. Will also test the T2 Psi3->MPQC transform
      static const bool mp2_only_ = false;
      /// set to true to use Ts instead of Lambdas
      static const bool replace_Lambda_with_T_ = true;
      /// default is to include up to 3rd-order terms in the energy
      static const unsigned int completeness_order_for_energy_ = 3;
      /// the max order for the intermediates is one less
      static const unsigned int
          completeness_order_for_intermediates_ = completeness_order_for_energy_
              - 1;

      double eccsd_;

      Ref<MBPT2_R12> mbptr12_;
    protected:
      void write_input(int conv);
    public:
      PsiCCSD_PT2R12(const Ref<KeyVal>&);
      PsiCCSD_PT2R12(StateIn&);
      ~PsiCCSD_PT2R12();
      void save_data_state(StateOut&);
      int gradient_implemented() const;
      /// Psi only produces the CCSD wave function
      void compute();
  };

}
#endif
