#ifndef _chemistry_qc_ci_ci_h
#define _chemistry_qc_ci_ci_h

#include <chemistry/qc/nbody/nbwfn.h>

namespace sc {

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

          <tr><td><tt>key1</tt><td>string<td>false<td>Keyword1 description

           <tr><td><tt>key2</tt><td>string<td>false<td>Keyword2 description

           </table>
       */
      CI(const Ref<KeyVal> &kv);
      CI(StateIn&);
      void save_data_state(StateOut &so);
      virtual ~CI();

      int value_implemented() const;

      void compute();

      RefSymmSCMatrix density();
    
    private:
      Ref<KeyVal> kv_;
      std::vector<double> E_;
      static std::vector<double> compute(const Ref<RefWavefunction> &ref,
                                         const Ref<KeyVal> &kv);
  };

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
