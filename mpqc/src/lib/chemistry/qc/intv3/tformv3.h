
#if defined(__GNUC__)
#pragma interface
#endif

#ifndef _chemistry_qc_intv3_tranform_h
#define _chemistry_qc_intv3_tranform_h

#include <chemistry/qc/basis/transform.h>
#include <chemistry/qc/intv3/macros.h>

extern "C" {
#include <chemistry/qc/intv2/atoms.h>
}

/* integrals and target may overlap */
void intv3_transform_1e(double *integrals, double *target,
                      shell_t *sh1, shell_t *sh2);

/* integrals and target may not overlap */
void intv3_accum_transform_1e(double *integrals, double *target,
                            shell_t *sh1, shell_t *sh2);

/* integrals and target may overlap */
void intv3_transform_1e_xyz(double *integrals, double *target,
                          shell_t *sh1, shell_t *sh2);

/* integrals and target may not overlap */
void intv3_accum_transform_1e_xyz(double *integrals, double *target,
                                shell_t *sh1, shell_t *sh2);

/* integrals and target may overlap */
void intv3_transform_2e(double *integrals, double *target,
                      shell_t *sh1, shell_t *sh2,
                      shell_t *sh3, shell_t *sh4);

class SphericalTransformComponentV3 : public SphericalTransformComponent {
  public:
    void init(int a, int b, int c, double coef, int pureindex) {
      a_ = a;
      b_ = b;
      c_ = c;
      coef_ = coef;
      pureindex_ = pureindex;
      cartindex_ = INT_CARTINDEX(a+b+c,a,b);
    }
};

class SphericalTransformV3 : public SphericalTransform {
  public:
    SphericalTransformV3(int l) {
      n_=0;
      l_=l;
      components_=0;
      init();
    }

    SphericalTransformComponent * new_components() {
      return new SphericalTransformComponentV3[n_+1];
    }
};

class ISphericalTransformV3 : public ISphericalTransform {
  public:
    ISphericalTransformV3(int l) {
      n_ = 0;
      l_ = l;
      components_ = 0;
      init();
    }

    SphericalTransformComponent * new_components() {
      return new SphericalTransformComponentV3[n_+1];
    }
};

class SphericalTransformIterV3 : public SphericalTransformIter {
  public:
    SphericalTransformIterV3(int l, int inverse=0);
};

#endif
