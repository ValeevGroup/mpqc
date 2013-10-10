#ifndef _chemistry_qc_ci_ci_h
#define _chemistry_qc_ci_ci_h

#include <chemistry/qc/nbody/nbwfn.h>

//#define MPQC_PROFILE_ENABLE
#include "mpqc/utility/profile.hpp"
#include "mpqc/ci/ci.hpp"

namespace sc {

  /// @addtogroup ChemistryElectronicStructureNBody
  /// @{

  /**
   * CI is a configuration interaction ManyBodyWavefunction. Currently only full CI is supported.
   */
  class CI: public ManyBodyWavefunction {
    public:

      /** A KeyVal constructor is used to generate a CI
          object from the input. This constructor accepts all keywords
          of the KeyVal constructor of the ManyBodyWavefunction class, plus the additional
          keywords listed below.

          <table border="1">

          <tr><td><b>%Keyword</b><td><b>Type</b><td><b>Default</b><td><b>Description</b>

           <tr><td><tt>total_charge</tt><td>int<td>see description<td>Specifies the total charge
           of the system. This charge is defined without taking into account custom nuclear
           charges or classical charges, i.e. total charge = sum of atomic numbers of nuclei
           - number of electrons.
           The default is to assume the same charge as that of the RefWavefunction object, i.e.
           @c molecule()->total_Z()-refwfn()->nelectron()

           <tr><td><tt>magnetic_moment</tt><td>int<td>see description<td>The S (or J) magnetic moment
           of the target state(s), in units of \f$ \hbar/2 \f$.
           The default is the magnetic moment of the @c reference RefWavefunction.

           <tr><td><tt>max_ex_rank</tt><td>unsigned int<td>0<td>The maximum excitation rank. The default
           is zero, which denotes full CI. This is equivalent to setting rank=number of electrons
           in active orbitals.

           </table>
       */
      CI(const Ref<KeyVal> &kv);
      CI(StateIn&);
      void save_data_state(StateOut &so);
      virtual ~CI();

      void print(std::ostream& os=ExEnv::out0()) const;

      int value_implemented() const;

      void compute();

      int nelectron();
      double magnetic_moment() const;
      RefSymmSCMatrix density();

    private:
      mpqc::ci::Config config_;
      std::vector<double> E_;
      static std::vector<double> compute(const Ref<RefWavefunction> &ref,
                                         const mpqc::ci::Config& config);
  };

  /// @}
  // end of addtogroup ChemistryElectronicStructureNBody

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
