#ifdef __GNUG__
#pragma interface
#endif

#ifndef _chemistry_qc_int1e_h
#define _chemistry_qc_int1e_h

#include <util/ref/ref.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/intv3/fjt.h>

extern "C" {
#include <chemistry/qc/intv2/atoms.h>
}

class Int1eV3: public VRefCount {
  protected:
    centers_t *cs1;
    centers_t *cs2;
    RefGaussianBasisSet bs1_;
    RefGaussianBasisSet bs2_;
    double *fjttable_;
    RefFJT fjt_;

    // statics from comp_1e.c:
  protected:
    double oo2zeta_a;
    double oo2zeta_b;
    double sMus;
    double sTs;
    double xi;
    double *A;
    double *B;
    double *C;
    double ss;
    double PmA[3];
    double PmB[3];
    double PmC[3];
    double zeta;
    double oo2zeta;
    shell_t *shell1, *shell2;
    int exponent_weighted;
    int scale_shell_result;
    double result_scale_factor;
    int three_center;
    centers_t *third_centers;
    int third_centernum;
    int init_order;
    double *buff;
    double *cartesianbuffer;
    int mu;

  protected:
    void accum_shell_1der(
        double *buff, int ish, int jsh,
        centers_t *dercs, int centernum,
        double (Int1eV3::*)(int,int,int,int,int,int,int,int)
        );
    double comp_shell_overlap(int gc1, int i1, int j1, int k1,
                              int gc2, int i2, int j2, int k2);
    double comp_prim_overlap(int i1, int j1, int k1,
                             int i2, int j2, int k2);
    double comp_shell_kinetic(int gc1, int i1, int j1, int k1,
                              int gc2, int i2, int j2, int k2);
    double comp_prim_kinetic(int i1, int j1, int k1,
                             int i2, int j2, int k2);
    double comp_shell_nuclear(int gc1, int i1, int j1, int k1,
                              int gc2, int i2, int j2, int k2);
    void accum_shell_efield(double *buff, int ish, int jsh);
    double comp_prim_nuclear(int i1, int j1, int k1,
                             int i2, int j2, int k2, int m);
    void comp_shell_efield(double *efield,
                           int gc1, int i1, int j1, int k1,
                           int gc2, int i2, int j2, int k2);
    double comp_prim_efield(int xyz, int i1, int j1, int k1,
                            int i2, int j2, int k2, int m);
    void comp_shell_dipole(double* dipole,
                           int gc1, int i1, int j1, int k1,
                           int gc2, int i2, int j2, int k2);
    double comp_prim_dipole(int im, int jm, int km,
                            int i1, int j1, int k1,
                            int i2, int j2, int k2);
    // routines from comp_1e:
  protected:
    void int_accum_shell_overlap_1der(int ish, int jsh,
                                      centers_t *dercs, int centernum);
    void int_done_1e();
    void int_initialize_1e(int flags, int order);
    double int_prim_overlap(shell_t *pshell1, shell_t *pshell2,
                            double *pA, double *pB,
                            int prim1, int prim2,
                            int i1, int j1, int k1,
                            int i2, int j2, int k2);
    void int_accum_shell_kinetic(int ish, int jsh);
    void int_accum_shell_kinetic_1der(int ish, int jsh,
                                      centers_t *dercs, int centernum);
    void int_accum_shell_nuclear_1der(int ish, int jsh,
                                      centers_t *dercs, int centernum);
    void int_accum_shell_nuclear_hfc_1der(int ish, int jsh,
                                          centers_t *dercs, int centernum);
    void int_accum_shell_nuclear_hf_1der(int ish, int jsh,
                                         centers_t *dercs, int centernum);
    void int_accum_shell_nuclear_nonhf_1der(int ish, int jsh,
                                            centers_t *dercs, int centernum);
    void int_accum_shell_efield(int ish, int jsh,
                                double *position);
    void int_accum_shell_point_charge(int ish, int jsh,
                                      int ncharge, double* charge,
                                      double** position);
    void int_shell_nuclear_hf_1der(int ish, int jsh,
                                   centers_t *dercs, int centernum);
    void int_shell_nuclear_nonhf_1der(int ish, int jsh,
                                      centers_t *dercs, int centernum);
    void int_accum_shell_dipole(int ish, int jsh,
                                double *com);

    // from offsets.cc
  protected:
    void int_initialize_offsets1(centers_t *cs1, centers_t *cs2);
    void int_done_offsets1(centers_t *cs1, centers_t *cs2);

  public:
    Int1eV3(const RefGaussianBasisSet&,
            const RefGaussianBasisSet&,
            int order);
    ~Int1eV3();

    double *buffer() { return buff; }
    RefGaussianBasisSet basis() { if (bs1_==bs2_) return bs1_; return 0; }
    RefGaussianBasisSet basis1() { return bs1_; }
    RefGaussianBasisSet basis2() { return bs2_; }

    void kinetic(int ish, int jsh);
    void nuclear(int ish, int jsh);
    void overlap(int ish, int jsh);
    void hcore(int ish, int jsh);
    void efield(int ish, int jsh, double position[3]);
    void point_charge(int ish, int jsh,
                      int ncharge, double* charge,
                      double** position);
    void dipole(int ish, int jsh,
                double *com);

    void hcore_1der(int ish, int jsh,
                    int dercs, int centernum);
    void kinetic_1der(int ish, int jsh,
                      int dercs, int centernum);
    void nuclear_1der(int ish, int jsh,
                      int dercs, int centernum);
    void overlap_1der(int ish, int jsh,
                      int dercs, int centernum);
};
REF_dec(Int1eV3);


#endif
