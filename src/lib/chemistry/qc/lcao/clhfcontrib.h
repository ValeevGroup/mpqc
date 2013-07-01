
#ifndef _chemistry_qc_lcao_clhfcontrib_h
#define _chemistry_qc_lcao_clhfcontrib_h

#include <mpqc_config.h>
#include <util/misc/exenv.h>
#include <chemistry/qc/lcao/fockbuild.h>

namespace sc {

/**
     Computes components of the Fock matrix necessary for closed-shell
     calculations (i.e. CLSCF). Requires as input the total AO density matrix P(0).
     Output matrices are in AO basis. If f_b1 == f_b2 then the output is the skeleton AO matrix
     that needs to be symmetrized with PetiteList, else the output is the full AO matrix.
 */
class CLHFContribution: public GenericFockContribution {
  public:
    CLHFContribution(const Ref<GaussianBasisSet> &f_b1,
                     const Ref<GaussianBasisSet> &f_b2,
                     const Ref<GaussianBasisSet> &p_b,
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
                     const double * RESTRICT buf);

    /** Compute the exchange contribution applying no two electron integral
        permutations. The integrals in buf are the full redundant set
        (nI*nJ*nK*nL integrals). The computes only the Coulomb contribution
        to the Fock matrix. */
    void contrib_e_K(double factor,
                     int I, int J, int K, int L,
                     int nI, int nJ, int nK, int nL,
                     const double * RESTRICT buf);

    void contrib_p12_p13p24_J(double factor,
                              int I, int J, int K, int L,
                              int nI, int nJ, int nK, int nL,
                              const double * RESTRICT buf);
    void contrib_p12_p13p24_K(double factor,
                              int I, int J, int K, int L,
                              int nI, int nJ, int nK, int nL,
                              const double * RESTRICT buf);
    void contrib_p34_p13p24_J(double factor,
                              int I, int J, int K, int L,
                              int nI, int nJ, int nK, int nL,
                              const double * RESTRICT buf);
    void contrib_p34_p13p24_K(double factor,
                              int I, int J, int K, int L,
                              int nI, int nJ, int nK, int nL,
                              const double * RESTRICT buf);
    void contrib_p12_p34_J(double factor,
                           int I, int J, int K, int L,
                           int nI, int nJ, int nK, int nL,
                           const double * RESTRICT buf);
    void contrib_p12_p34_K(double factor,
                           int I, int J, int K, int L,
                           int nI, int nJ, int nK, int nL,
                           const double * RESTRICT buf);
    void contrib_p34_J(double factor,
                       int I, int J, int K, int L,
                       int nI, int nJ, int nK, int nL,
                       const double * RESTRICT buf);
    void contrib_p34_K(double factor,
                       int I, int J, int K, int L,
                       int nI, int nJ, int nK, int nL,
                       const double * RESTRICT buf);
    void contrib_p13p24_J(double factor,
                          int I, int J, int K, int L,
                          int nI, int nJ, int nK, int nL,
                          const double * RESTRICT buf);
    void contrib_p13p24_K(double factor,
                          int I, int J, int K, int L,
                          int nI, int nJ, int nK, int nL,
                          const double * RESTRICT buf);

    /** Compute the Coulomb contribution applying all two electron integral
        permutations. I, J, K, and L indices must all be unique.
    */
    void contrib_all_J(double factor,
                       int I, int J, int K, int L,
                       int nI, int nJ, int nK, int nL,
                       const double * RESTRICT buf);

    /** Compute the exchange contribution applying all two electron integral
        permutations. I, J, K, and L indices must all be unique.
    */
    void contrib_all_K(double factor,
                       int I, int J, int K, int L,
                       int nI, int nJ, int nK, int nL,
                       const double * RESTRICT buf);
};

}

#endif
