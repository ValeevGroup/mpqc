#ifdef __GNUG__
#pragma interface
#endif

#ifndef _chemistry_qc_int2e_h
#define _chemistry_qc_int2e_h

#include <util/ref/ref.h>
#include <chemistry/qc/basis/basis.h>
#include <chemistry/qc/oint3/build.h>
#include <chemistry/qc/intv3/fjt.h>
#include <chemistry/qc/intv2/atoms.h>
#include <chemistry/qc/intv3/types.h>
#include <chemistry/qc/intv2/storage.h>

class Int2eV3: public VRefCount {
  protected:
    BuildIntV3 build;
    RefIntegralStorer storer;

    RefGaussianBasisSet bs1_;
    RefGaussianBasisSet bs2_;
    RefGaussianBasisSet bs3_;
    RefGaussianBasisSet bs4_;

    int bs1_shell_offset_;
    int bs2_shell_offset_;
    int bs3_shell_offset_;
    int bs4_shell_offset_;
    int bs1_func_offset_;
    int bs2_func_offset_;
    int bs3_func_offset_;
    int bs4_func_offset_;
    int bs1_prim_offset_;
    int bs2_prim_offset_;
    int bs3_prim_offset_;
    int bs4_prim_offset_;

    // statics from vrr.cc
  public:
    enum { STORAGE_CHUNK = 4096 };
  protected:
    struct store_list {
        void* data[STORAGE_CHUNK];
        struct store_list* p;
    };
    typedef struct store_list store_list_t;
    int n_store_last;
    store_list_t* store;
    typedef int (BuildIntV3::*intfunc)();
    intfunc build_routine[4][4][4][4][2];
    /* Offset shell numbers. */
    int osh1, osh2, osh3, osh4;
    /* Offset primitive numbers. */
    int opr1, opr2, opr3, opr4;
    /* Boolean array which gives whether or not an array is computed. */
    int_array3_t inthave;
    /* Saved initialization parameters used to free data. */
    int saved_am12,saved_am34,saved_ncon;
    /* Stores the length of the inner loop for integral contraction. */
    int_array3_t contract_length;

    // statics from hrr.cc
  protected:
    /* The general contraction numbers. */
    int g1,g2,g3,g4;
    /* A[] - B[] */
    double AmB[3];
    /* C[] - D[] */
    double CmD[3];
    /* Boolean array which gives whether or not a set of integrals has been
     * computed. */
    int_array4_t shiftinthave;
    int eAB;

    int redundant_;
    int permute_;

  protected:
    RefFJT fjt_;

    int_vector_t int_shell_to_prim;
    double_matrix_t int_shell_r;
    double_matrix_t int_prim_zeta;
    double_matrix_t int_prim_k;
    double_matrix_t int_prim_oo2zeta;
    double_array3_t int_prim_p;

    double *int_buffer;
    double *int_derint_buffer;

    RefGaussianBasisSet int_cs1;
    RefGaussianBasisSet int_cs2;
    RefGaussianBasisSet int_cs3;
    RefGaussianBasisSet int_cs4;

    GaussianShell *int_shell1;
    GaussianShell *int_shell2;
    GaussianShell *int_shell3;
    GaussianShell *int_shell4;

    doublep_array4_t *int_con_ints;
    doublep_array4_t ****int_con_ints_array;  /* The contr. int. inter. */

    int int_expweight1; // For exponent weighted contractions.
    int int_expweight2; // For exponent weighted contractions.
    int int_expweight3; // For exponent weighted contractions.
    int int_expweight4; // For exponent weighted contractions.

    // These are used to compute two and three center electron repulsion
    // integrals.  int_unit2 is 1 if shell 2 is to have value one everywhere
    // and int_unit4 is 1 if shell4 is to be a unit function.  Otherwise,
    // they should be zero.
    //

    int int_unit2;
    int int_unit4;
    GaussianShell* int_unit_shell;

    int int_integral_storage;
    int int_store1;
    int int_store2;
    int int_derivative_bounds;

    // locals from vrr.cc
  protected:
    void add_store(void *p);
    void free_store();
    void _free_store(store_list_t* s, int n);
    void build_not_using_gcs(int nc1, int nc2, int nc3, int nc4,
                             int minam1, int minam3, int maxam12, int maxam34,
                             int dam1, int dam2, int dam3, int dam4, int eAB);
    void build_using_gcs(int nc1, int nc2, int nc3, int nc4,
                         int minam1, int minam3, int maxam12, int maxam34,
                         int dam1, int dam2, int dam3, int dam4, int eAB);
    void gen_prim_intermediates(int pr1, int pr2, int pr3, int pr4, int am);
    void gen_prim_intermediates_with_norm(int pr1, int pr2, int pr3, int pr4,
                                 int am, double norm);
    void gen_shell_intermediates(int sh1, int sh2, int sh3, int sh4);
    double * buildprim(int am12, int am34, int m);
    void buildprim_1(double *I00, int am12, int am34, int m);
    void buildprim_3(double *I00, int am12, int am34, int m);
    void init_inthave(int am12, int am34);
    int choose_center(int am12, int am34, int m);

    // globals from vrr.cc
  protected:
    void int_init_buildgc(int order,
                          int am1, int am2, int am3, int am4,
                          int nc1, int nc2, int nc3, int nc4);
    void int_done_buildgc();
    void int_buildgcam(int minam1, int minam2, int minam3, int minam4,
                       int maxam1, int maxam2, int maxam3, int maxam4,
                       int dam1, int dam2, int dam3, int dam4,
                       int sh1, int sh2, int sh3, int sh4,
                       int eAB);

    // globals from print2e.cc
  protected:
    void int_offset_print(FILE *fp,
                          double *buffer,
                          RefGaussianBasisSet c1, int s1,
                          RefGaussianBasisSet c2, int s2,
                          RefGaussianBasisSet c3, int s3,
                          RefGaussianBasisSet c4, int s4);
    void int_offset_print_n(FILE *fp, double *buffer,
                            int n1, int n2, int n3, int n4,
                            int o1, int o2, int o3, int o4,
                            int e12, int e13e24, int e34);
    void int_print(FILE *fp, double *buffer,
                   RefGaussianBasisSet c1, int s1,
                   RefGaussianBasisSet c2, int s2,
                   RefGaussianBasisSet c3, int s3,
                   RefGaussianBasisSet c4, int s4);
    void int_print_n(FILE *fp, double *buffer,
                     int n1, int n2, int n3, int n4,
                     int e12, int e13e24, int e34);
    void int_print_intermediates(FILE *fp);

    // locals from hrr.cc
  protected:
    void init_shiftinthave(int am1, int am2, int am3, int am);
    double * shiftint(int am1, int am2, int am3, int am4);
    int choose_shift(int am1, int am2, int am3, int am4);
    void shiftam_12(double *I0100, int am1, int am2, int am3, int am4);
    void shiftam_12eAB(double *I0100, int am1, int am2, int am3, int am4);
    void shiftam_34(double *I0001, int am1, int am2, int am3, int am4);
        
    // globals from hrr.cc
  protected:
    void int_init_shiftgc(int order, int am1, int am2, int am3, int am4);
    void int_done_shiftgc();
    void int_shiftgcam(int gc1, int gc2, int gc3, int gc4,
                       int tam1, int tam2, int tam3, int tam4, int peAB);

    // locals from init2e.cc
  protected:
    void alloc_inter(int nprim,int nshell);
    void compute_shell_1(RefGaussianBasisSet cs, int, int);
    void compute_prim_1(RefGaussianBasisSet cs1);
    void compute_shell_2(RefGaussianBasisSet cs1,RefGaussianBasisSet cs2);
    void compute_prim_2(RefGaussianBasisSet cs1,RefGaussianBasisSet cs2);


    // globals from init2e.cc
  protected:
    double *int_initialize_erep(int storage, int order,
                                RefGaussianBasisSet cs1,
                                RefGaussianBasisSet cs2,
                                RefGaussianBasisSet cs3,
                                RefGaussianBasisSet cs4);
    void int_done_erep();

    // locals from comp2e.cc
  protected:
    void compute_erep(int flags, int *psh1, int *psh2, int *psh3, int *psh4,
                      int dam1, int dam2, int dam3, int dam4);
    void compute_erep_1der(int flags, double *buffer,
                           int *psh1, int *psh2, int *psh3, int *psh4,
                           int dercenter);
    void nonredundant_erep(double *buffer, int e12, int e34, int e13e24,
                           int n1, int n2, int n3, int n4,
                           int *red_off, int *nonred_off);
    void compute_erep_bound1der(int flags, double *buffer,
                                int *psh1, int *psh2, int *psh3, int *psh4);

    // globals from comp2e.cc
  protected:
    void int_erep_bound1der(int flags, int bsh1, int bsh2, int *size);


    // global vars from bounds.h
  protected:
    typedef signed char int_bound_t;
    int_bound_t int_Q;    
    int_bound_t int_R;    
    int_bound_t *int_Qvec;
    int_bound_t *int_Rvec;

    // global routines from bounds.cc
  protected:
    void int_init_bounds_nocomp();
    void int_init_bounds_1der_nocomp();
    void int_bounds_comp(int s1, int s2);
    void int_bounds_1der_comp(int s1, int s2);
    int int_erep_2bound(int s1, int s2);
    int int_erep_0bound_1der();
    int int_erep_2bound_1der(int s1, int s2);

    // local routines from bounds.cc
  protected:
    void compute_bounds(int_bound_t *overall, int_bound_t *vec, int flag);
    void compute_bounds_shell(int_bound_t *overall, int_bound_t *vec,
                              int flag, int sh1, int sh2);

    // global routines from storage.cc
  protected:
    int int_have_stored_integral(int sh1,int sh2,int sh3,int sh4,
                                 int p12,int p34,int p13p24);
    void int_store_integral(int sh1,int sh2,int sh3,int sh4,
                            int p12,int p34,int p13p24,
                            int size);

    // from offsets.cc
  protected:
    void int_initialize_offsets2();
    void int_done_offsets2();

    // from comp2e3c.cc
  protected:
    void make_int_unit_shell();
    void delete_int_unit_shell();

  protected:
    int used_storage_;
    int used_storage_build_;

  public:
    Int2eV3(const RefGaussianBasisSet&,
            const RefGaussianBasisSet&,
            const RefGaussianBasisSet&,
            const RefGaussianBasisSet&,
            int order, int storage);
    ~Int2eV3();

    // storage.cc
    void init_storage(int size);
    void done_storage();

    // bounds.cc
    void init_bounds();
    void init_bounds_1der();
    void done_bounds();
    void done_bounds_1der();
    // Covert a bound to/from the log of the bound (returns 2^bound)
    // replace:
    //double int_bound_to_double(int bound);
    //double int_bound_double(int value);
    //int int_bound_log(double value);
    static double logbound_to_bound(int);
    static int bound_to_logbound(double value);

    // If redundant is false the redundant integrals that arise
    // when a shell index is repeated are stored.
    // The default is true.
    int redundant() { return redundant_; }
    void set_redundant(int i) { redundant_ = i; }

    // If permute is true the routines are allowed to permute indices.
    // The default is false.
    int permute() { return permute_; }
    void set_permute(int i) { permute_ = i; }

    int used_storage() const { return used_storage_; }

    // from comp2e.cc
    void erep(int &psh1, int &psh2, int &psh3, int &psh4);
    void erep(int *shells, int  *sizes);
    void erep_all1der(int &psh1, int &psh2, int &psh3, int &psh4,
                      der_centersv3_t *der_centers);
    void erep_all1der(int *shells, int  *sizes,
                      der_centersv3_t *dercenters);

    // from comp2e3c.cc
    void erep_2center(int &psh1, int &psh2);
    void erep_2center(int *shells, int  *sizes);
    void erep_3center(int &psh1, int &psh2, int &psh3);
    void erep_3center(int *shells, int  *sizes);

    // from bounds.cc
    int erep_4bound(int s1, int s2, int s3, int s4);
    int erep_4bound_1der(int s1, int s2, int s3, int s4);

    double *buffer() { return int_buffer; }

    RefGaussianBasisSet basis()
    {
      if (bs1_==bs2_ && bs1_ == bs3_ && bs1_ == bs4_) return bs1_;
      return 0;
    }
    RefGaussianBasisSet basis1() { return bs1_; }
    RefGaussianBasisSet basis2() { return bs2_; }
    RefGaussianBasisSet basis3() { return bs3_; }
    RefGaussianBasisSet basis4() { return bs4_; }

    RefGaussianBasisSet cs1() const { return int_cs1; }
    RefGaussianBasisSet cs2() const { return int_cs2; }
    RefGaussianBasisSet cs3() const { return int_cs3; }
    RefGaussianBasisSet cs4() const { return int_cs4; }
};
REF_dec(Int2eV3);

#endif
