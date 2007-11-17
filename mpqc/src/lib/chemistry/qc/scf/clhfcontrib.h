
#include <scconfig.h>
#include <util/misc/exenv.h>
#include <chemistry/qc/scf/fockbuild.h>

namespace sc {

class CLHFContribution: public GenericFockContribution {
  public:
    CLHFContribution(Ref<GaussianBasisSet> &f_b1,
                     Ref<GaussianBasisSet> &f_b2,
                     Ref<GaussianBasisSet> &p_b,
                     const std::string &fockbuildmatrixtype);
    CLHFContribution(const CLHFContribution &);
    ~CLHFContribution();

    Ref<FockContribution> clone();

    /** Compute the Coulomb contribution applying no two electron integral
        permutations. The integrals in buf are the full redundant set
        (nI*nJ*nK*nL integrals). The computes only the Coulomb contribution
        to the Fock matrix. */
    void contrib_e_J(double factor,
                     int I, int J, int K, int L,
                     int nI, int nJ, int nK, int nL,
                     const double * restrictxx buf);

    /** Compute the exchange contribution applying no two electron integral
        permutations. The integrals in buf are the full redundant set
        (nI*nJ*nK*nL integrals). The computes only the Coulomb contribution
        to the Fock matrix. */
    void contrib_e_K(double factor,
                     int I, int J, int K, int L,
                     int nI, int nJ, int nK, int nL,
                     const double * restrictxx buf);

    void contrib_p12_p13p24_J(double factor,
                              int I, int J, int K, int L,
                              int nI, int nJ, int nK, int nL,
                              const double * restrictxx buf);
    void contrib_p12_p13p24_K(double factor,
                              int I, int J, int K, int L,
                              int nI, int nJ, int nK, int nL,
                              const double * restrictxx buf);
    void contrib_p34_p13p24_J(double factor,
                              int I, int J, int K, int L,
                              int nI, int nJ, int nK, int nL,
                              const double * restrictxx buf);
    void contrib_p34_p13p24_K(double factor,
                              int I, int J, int K, int L,
                              int nI, int nJ, int nK, int nL,
                              const double * restrictxx buf);
    void contrib_p12_p34_J(double factor,
                           int I, int J, int K, int L,
                           int nI, int nJ, int nK, int nL,
                           const double * restrictxx buf);
    void contrib_p12_p34_K(double factor,
                           int I, int J, int K, int L,
                           int nI, int nJ, int nK, int nL,
                           const double * restrictxx buf);
    void contrib_p34_J(double factor,
                       int I, int J, int K, int L,
                       int nI, int nJ, int nK, int nL,
                       const double * restrictxx buf);
    void contrib_p34_K(double factor,
                       int I, int J, int K, int L,
                       int nI, int nJ, int nK, int nL,
                       const double * restrictxx buf);
    void contrib_p13p24_J(double factor,
                          int I, int J, int K, int L,
                          int nI, int nJ, int nK, int nL,
                          const double * restrictxx buf);
    void contrib_p13p24_K(double factor,
                          int I, int J, int K, int L,
                          int nI, int nJ, int nK, int nL,
                          const double * restrictxx buf);

    /** Compute the Coulomb contribution applying all two electron integral
        permutations. I, J, K, and L indices must all be unique.
    */
    void contrib_all_J(double factor,
                       int I, int J, int K, int L,
                       int nI, int nJ, int nK, int nL,
                       const double * restrictxx buf);

    /** Compute the exchange contribution applying all two electron integral
        permutations. I, J, K, and L indices must all be unique.
    */
    void contrib_all_K(double factor,
                       int I, int J, int K, int L,
                       int nI, int nJ, int nK, int nL,
                       const double * restrictxx buf);
};

}
