/*
 * tformv2.h
 *
 * Copyright (C) 1996 Limit Point Systems, Inc.
 *
 * Author: Curtis Janssen <cljanss@ca.sandia.gov>
 * Maintainer: LPS
 *
 * This file is part of the SC Toolkit.
 *
 * The SC Toolkit is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Library General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * The SC Toolkit is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public License
 * along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
 * the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 * The U.S. Government is granted a limited license as per AL 91-7.
 */

#if defined(__cplusplus) && defined(__GNUC__)
#pragma interface
#endif

#ifndef _chemistry_qc_intv2_tranform_h
#define _chemistry_qc_intv2_tranform_h

#include <chemistry/qc/intv2/atoms.h>

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

#include <chemistry/qc/basis/transform.h>
#include <chemistry/qc/intv2/int_macros.h>

class SphericalTransformComponentV2 : public SphericalTransformComponent {
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

class SphericalTransformV2 : public SphericalTransform {
  public:
    SphericalTransformV2(int l) {
      n_=0;
      l_=l;
      components_=0;
      init();
    }

    SphericalTransformComponent * new_components() {
      return new SphericalTransformComponentV2[n_+1];
    }
};

class ISphericalTransformV2 : public ISphericalTransform {
  public:
    ISphericalTransformV2(int l) {
      n_ = 0;
      l_ = l;
      components_ = 0;
      init();
    }

    SphericalTransformComponent * new_components() {
      return new SphericalTransformComponentV2[n_+1];
    }
};

class SphericalTransformIterV2 : public SphericalTransformIter {
  public:
    SphericalTransformIterV2(int l, int inverse=0);
};

#endif /* __c_plus_plus */

#endif
