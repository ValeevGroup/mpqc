
#if defined(__cplusplus) && defined(__GNUC__)
#pragma interface
#endif

#ifndef _chemistry_qc_intv2_tranform_h
#define _chemistry_qc_intv2_tranform_h

#include <chemistry/qc/intv2/atoms.h>

#ifdef __cplusplus

class SphericalTransform {
  public:
    class Component {
      private:
        double coef_;
        int a_, b_, c_, cartindex_, pureindex_;
      public:
        int a() const { return a_; }
        int b() const { return b_; }
        int c() const { return c_; }
        int cartindex() const { return cartindex_; }
        int pureindex() const { return pureindex_; }
        double coef() const { return coef_; }
        void init(int a, int b, int c, double coef, int pureindex);
    };
  protected:
    int n_;
    int l_;
    Component *components_;
    SphericalTransform();
  public:
    SphericalTransform(int l);
    ~SphericalTransform();
    int cartindex(int i) const { return components_[i].cartindex(); }
    int pureindex(int i) const { return components_[i].pureindex(); }
    double coef(int i) const { return components_[i].coef(); }
    int a(int i) const { return components_[i].a(); }
    int b(int i) const { return components_[i].b(); }
    int c(int i) const { return components_[i].c(); }
    int l() const { return l_; }
    int n() const { return n_; }

    void add(int a, int b, int c, double coef, int pureindex);
};

// The inverse transforms
class ISphericalTransform: public SphericalTransform {
  public:
    ISphericalTransform(int l);
};

class SphericalTransformIter {
  private:
    const SphericalTransform *transform_;

    int i_;
  public:
    SphericalTransformIter(SphericalTransform*);
    SphericalTransformIter(int l, int inverse = 0);
    void begin() { i_ = 0; }
    void start() { begin(); }
    void next() { i_++; }
    int ready() { return i_ < transform_->n(); }
    operator int() { return ready(); }
    int l() { return transform_->l(); }
    int cartindex() { return transform_->cartindex(i_); }
    int pureindex() { return transform_->pureindex(i_); }
    int bfn() { return pureindex(); }
    double coef() { return transform_->coef(i_); }
    int a() { return transform_->a(i_); }
    int b() { return transform_->b(i_); }
    int c() { return transform_->c(i_); }
    int l(int i) { return i?(i==1?b():c()):a(); }
    int n() { return 2*l() + 1; }
};

#endif

#ifdef __cplusplus
extern "C" {
#endif

/* integrals and target may overlap */
void int_transform_1e(double *integrals, double *target,
                      shell_t *sh1, shell_t *sh2);

/* integrals and target may not overlap */
void int_accum_transform_1e(double *integrals, double *target,
                            shell_t *sh1, shell_t *sh2);

/* integrals and target may overlap */
void int_transform_1e_xyz(double *integrals, double *target,
                          shell_t *sh1, shell_t *sh2);

/* integrals and target may not overlap */
void int_accum_transform_1e_xyz(double *integrals, double *target,
                                shell_t *sh1, shell_t *sh2);

/* integrals and target may overlap */
void int_transform_2e(double *integrals, double *target,
                      shell_t *sh1, shell_t *sh2,
                      shell_t *sh3, shell_t *sh4);

#ifdef __cplusplus
}
#endif

#endif
