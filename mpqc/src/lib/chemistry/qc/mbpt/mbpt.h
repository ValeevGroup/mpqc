
#ifndef _chemistry_qc_mbpt_mbpt_h
#define _chemistry_qc_mbpt_mbpt_h

#ifdef __GNUC__
#pragma interface
#endif

#include <util/group/memory.h>
#include <util/group/message.h>
#include <chemistry/qc/basis/obint.h>
#include <chemistry/qc/basis/tbint.h>
#include <chemistry/qc/scf/scf.h>

////////////////////////////////////////////////////////////////////////////

class MBPT2: public Wavefunction {
#   define CLASSNAME MBPT2
#   define HAVE_KEYVAL_CTOR
#   define HAVE_STATEIN_CTOR
#   include <util/state/stated.h>
#   include <util/class/classd.h>
  protected:
    int debug_;
    RefSCF reference_;
    RefMemoryGrp mem;
    int nfzc, nfzv;
    int mem_alloc;

    int eliminate_in_gmat_;
    const double *intbuf_;
    RefTwoBodyInt tbint_;
    RefTwoBodyDerivInt tbintder_;
    int nbasis;
    RefMessageGrp msg_;
    int nvir, nocc, nsocc;

    char *method_;

  protected:
    void init_variables();

    // implement the Compute::compute() function
    void compute();

    // Fill in the eigenvectors and eigenvalues (Guest & Saunders general
    // form is used for the Fock matrix in the open shell case).
    void eigen(RefDiagSCMatrix &vals, RefSCMatrix &vecs);

    // calculate the opt2 energy using algorithm v1
    void compute_hsos_v1();

    // calculate the opt2 energy using algorithm v2
    void compute_hsos_v2();

    // calculate the opt2 energy using the load balanced version of v2
    //void compute_hsos_v2_lb();

    // calculate the closed shell mp2 energy and gradient
    int compute_cs_batchsize(int mem_static, int nocc_act);
    int make_cs_gmat(RefSymmSCMatrix& Gmat, double *DPmat);
    void form_max_dens(double *DPmat, signed char *maxp);
    int init_cs_gmat();
    void done_cs_gmat();
    int make_g_d_nor(RefSymmSCMatrix& Gmat,
                     double *DPmat, const double *mgdbuff);
    void cs_cphf(double **scf_vector,
                 double *Laj, double *eigval, RefSCMatrix& P2aj);
    void s2pdm_contrib(const double *intderbuf, double *PHF,
                       double *P2AO, double **ginter);
    void hcore_cs_grad(double *PMP2, double **ginter);
    void overlap_cs_grad(double *WMP2, double **ginter);
    void compute_cs_grad();
  public:
    MBPT2(StateIn&);
    MBPT2(const RefKeyVal&);
    ~MBPT2();

    void save_data_state(StateOut&);

    RefSymmSCMatrix density();

    int gradient_implemented();
    int value_implemented();

    void print(ostream&o=cout);
};
SavableState_REF_dec(MBPT2);

#endif

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
// End:
