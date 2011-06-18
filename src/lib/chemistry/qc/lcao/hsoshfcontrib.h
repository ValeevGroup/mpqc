
#ifndef _chemistry_qc_lcao_hsoshfcontrib_h
#define _chemistry_qc_lcao_hsoshfcontrib_h

#include <scconfig.h>
#include <util/misc/exenv.h>
#include <chemistry/qc/lcao/fockbuild.h>

namespace sc {

/**
   Computes components of the Fock matrix necessary for high-spin open-shell
   calculations (e.g. HSOSSCF and UnrestrictedSCF). Requires as input
   two density matrices P: P(0) is the total density and P(1) is the spin density
   (P(alpha) - P(beta)). Spin density is only used to compute a contribution to the exchange matrix.
   The resulting matrices can then be combined to produce alpha and beta components of J, K, and F matrices
   matrices as follows: J(alpha) = J(beta) = J(0), K(alpha) = K(0) + K(1), K(beta) = K(0) - K(1),
   F(alpha) = F(0) + F(1), F(beta) = F(0) - F(1). The corresponding closed-shell (F(c)) and open-shell (F(o))
   matrices used in HSOSSCF can be obtained as follows: F(o) = F(alpha), F(c) = (F(alpha) + F(beta))/2.

   Output matrices are in AO basis. If f_b1 == f_b2 then the output is the skeleton AO matrix
   that needs to be symmetrized with PetiteList, else the output is the full AO matrix.
 */
class HSOSHFContribution: public GenericFockContribution {
  public:
    HSOSHFContribution(const Ref<GaussianBasisSet> &f_b1,
                     const Ref<GaussianBasisSet> &f_b2,
                     const Ref<GaussianBasisSet> &p_b,
                     const std::string &fockbuildmatrixtype);
    HSOSHFContribution(const HSOSHFContribution &);
    ~HSOSHFContribution();

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

#endif
