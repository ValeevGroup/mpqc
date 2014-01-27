//
// aoints.cc
//
// Copyright (C) 2012 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
// Maintainer: EV
//
// This file is part of the SC Toolkit.
//
// The SC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The SC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#include <iostream>
#include <iomanip>
#include <cassert>
#include <chemistry/molecule/molecule.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/split.h>

#define SKIP_DERIV_INTS 0 // change to 0 to compute derivative integrals also

using namespace sc;

std::string xyzstring(unsigned int xyz) {
  switch(xyz) {
    case 0: return "x";
    case 1: return "y";
    case 2: return "z";
    default: MPQC_ASSERT(false);
  }
  return "";   // unreachable
}

int
main(int argc, char** argv) {

  //
  // construct a Molecule object
  //
  Ref<Molecule> mol = new Molecule;
  mol->add_atom(8, 0.0, 0.0, 0.0);
  mol->add_atom(1, 0.0, 1.0, 1.0);
  mol->add_atom(1, 0.0,-1.0, 1.0);
  const bool use_symmetry = true;
  if (not use_symmetry) {
    Ref<PointGroup> c1_ptgrp = new PointGroup("C1");
    mol->set_point_group(c1_ptgrp);
  }
  else {
    mol->symmetrize(mol->highest_point_group(1e-4));
  }

  ExEnv::out0() << std::endl << indent << "constructed Molecule object:" << std::endl;
  mol->print(ExEnv::out0());
  ExEnv::out0() << std::endl;

  //
  // construct a GaussianBasisSet object
  //
  Ref<AssignedKeyVal> akv = new AssignedKeyVal;
  akv->assign("molecule", mol.pointer());
  akv->assign("name", "cc-pVDZ");
  Ref<GaussianBasisSet> obs = new GaussianBasisSet(Ref<KeyVal>(akv));
  // get rid of general constractions for simplicity
  if (obs->max_ncontraction() > 1) {
    Ref<GaussianBasisSet> split_basis = new SplitBasisSet(obs, obs->name());
    obs = split_basis;
  }
  const int nshell = obs->nshell();

  ExEnv::out0() << std::endl << indent << "constructed GaussianBasisSet object:" << std::endl;
  obs->print(ExEnv::out0());
  ExEnv::out0() << std::endl;

  //
  // construct an Integral object
  // it will produce integral evaluator objects
  //
  Ref<Integral> integral = Integral::initial_integral(argc, argv);
  if (integral) Integral::set_default_integral(integral);
  integral = Integral::get_default_integral()->clone();
  integral->set_basis(obs);


  //
  // compute overlap integrals
  //

  // construct an OneBodyInt object that computes overlap integrals
  Ref<OneBodyInt> s_inteval = integral->overlap();
  const double* buffer = s_inteval->buffer();
  // and compute overlap integrals
  std::cout << "overlap integrals:" << std::endl;
  for(int s1=0; s1<nshell; s1++) {
    const int bf1_offset = obs->shell_to_function(s1);
    const int nbf1 = obs->shell(s1).nfunction();

    for(int s2=0; s2<nshell; s2++) {
      const int bf2_offset = obs->shell_to_function(s2);
      const int nbf2 = obs->shell(s2).nfunction();

      s_inteval->compute_shell(s1, s2);

      int bf12 = 0;
      for(int bf1=0; bf1<nbf1; ++bf1) {
        for(int bf2=0; bf2<nbf2; ++bf2, ++bf12) {
          std::cout << bf1+bf1_offset << " " << bf2+bf2_offset << " "
              << std::setprecision(15) << buffer[bf12] << std::endl;
        }
      }
    }
  }
  s_inteval = 0;


  //
  // compute core Hamiltonian integrals
  //

  // construct an OneBodyInt object that computes overlap integrals
  Ref<OneBodyInt> h_inteval = integral->hcore();
  buffer = h_inteval->buffer();
  // and compute core Hamiltonian integrals
  std::cout << "hcore integrals:" << std::endl;
  for(int s1=0; s1<nshell; s1++) {
    const int bf1_offset = obs->shell_to_function(s1);
    const int nbf1 = obs->shell(s1).nfunction();

    for(int s2=0; s2<nshell; s2++) {
      const int bf2_offset = obs->shell_to_function(s2);
      const int nbf2 = obs->shell(s2).nfunction();

      h_inteval->compute_shell(s1, s2);

      int bf12 = 0;
      for(int bf1=0; bf1<nbf1; ++bf1) {
        for(int bf2=0; bf2<nbf2; ++bf2, ++bf12) {
          std::cout << bf1+bf1_offset << " " << bf2+bf2_offset << " "
              << std::setprecision(15) << buffer[bf12] << std::endl;
        }
      }
    }
  }
  h_inteval = 0;

  //
  // compute 2-e Coulomb integrals
  //

  Ref<TwoBodyInt> twoecoulomb_inteval = integral->electron_repulsion();
  buffer = twoecoulomb_inteval->buffer();
  std::cout << "two-e Coulomb integrals:" << std::endl;
  for(int s1=0; s1<nshell; s1++) {
    const int bf1_offset = obs->shell_to_function(s1);
    const int nbf1 = obs->shell(s1).nfunction();

    for(int s2=0; s2<nshell; s2++) {
      const int bf2_offset = obs->shell_to_function(s2);
      const int nbf2 = obs->shell(s2).nfunction();

      for(int s3=0; s3<nshell; s3++) {
        const int bf3_offset = obs->shell_to_function(s3);
        const int nbf3 = obs->shell(s3).nfunction();

        for(int s4=0; s4<nshell; s4++) {
          const int bf4_offset = obs->shell_to_function(s4);
          const int nbf4 = obs->shell(s4).nfunction();

          twoecoulomb_inteval->compute_shell(s1, s2, s3, s4);

          int bf1234 = 0;
          for(int bf1=0; bf1<nbf1; ++bf1) {
            for(int bf2=0; bf2<nbf2; ++bf2) {
              for(int bf3=0; bf3<nbf3; ++bf3) {
                for(int bf4=0; bf4<nbf4; ++bf4, ++bf1234) {
                  std::cout << bf1+bf1_offset << " " << bf2+bf2_offset << " "
                      << bf3+bf3_offset << " " << bf4+bf4_offset << " "
                      << std::setprecision(15) << buffer[bf1234] << std::endl;
                }
              }
            }
          }

        }
      }
    }
  }
  twoecoulomb_inteval = 0;

#if !SKIP_DERIV_INTS
  //
  // compute geometric derivatives of overlap integrals
  //

  // construct an OneBodyDerivInt object that computes derivative overlap integrals
  Ref<OneBodyDerivInt> s_deriv_inteval = integral->overlap_deriv();
  buffer = s_deriv_inteval->buffer();
  // and compute derivative core Hamiltonian integrals
  std::cout << "overlap deriv integrals:" << std::endl;
  for(int atom=0; atom<mol->natom(); ++atom) {

    for(int s1=0; s1<nshell; s1++) {
      const int bf1_offset = obs->shell_to_function(s1);
      const int nbf1 = obs->shell(s1).nfunction();

      for(int s2=0; s2<nshell; s2++) {
        const int bf2_offset = obs->shell_to_function(s2);
        const int nbf2 = obs->shell(s2).nfunction();

        s_deriv_inteval->compute_shell(s1, s2, atom);

        int bf12 = 0;
        int bf12_deriv = 0;
        for(int bf1=0; bf1<nbf1; ++bf1) {
          for(int bf2=0; bf2<nbf2; ++bf2, ++bf12) {
            for(int xyz=0; xyz<3; ++xyz, ++bf12_deriv) {
              std::cout << bf1+bf1_offset << " " << bf2+bf2_offset
                        << " atom" << atom << "_" << xyzstring(xyz) << " "
                        << std::setprecision(15) << buffer[bf12_deriv] << std::endl;
            }
          }
        }
      }
    }
  }
  s_deriv_inteval = 0;

  //
  // compute geometric derivatives of core Hamiltonian integrals
  //

  // construct an OneBodyInt object that computes overlap integrals
  Ref<OneBodyDerivInt> h_deriv_inteval = integral->hcore_deriv();
  buffer = h_deriv_inteval->buffer();
  // and compute derivative core Hamiltonian integrals
  std::cout << "hcore deriv integrals:" << std::endl;
  for(int atom=0; atom<mol->natom(); ++atom) {

    for(int s1=0; s1<nshell; s1++) {
      const int bf1_offset = obs->shell_to_function(s1);
      const int nbf1 = obs->shell(s1).nfunction();

      for(int s2=0; s2<nshell; s2++) {
        const int bf2_offset = obs->shell_to_function(s2);
        const int nbf2 = obs->shell(s2).nfunction();

        h_deriv_inteval->compute_shell(s1, s2, atom);

        int bf12 = 0;
        int bf12_deriv = 0;
        for(int bf1=0; bf1<nbf1; ++bf1) {
          for(int bf2=0; bf2<nbf2; ++bf2, ++bf12) {
            for(int xyz=0; xyz<3; ++xyz, ++bf12_deriv) {
              std::cout << bf1+bf1_offset << " " << bf2+bf2_offset
                        << " atom" << atom << "_" << xyzstring(xyz) << " "
                        << std::setprecision(15) << buffer[bf12_deriv] << std::endl;
            }
          }
        }
      }
    }
  }
  h_deriv_inteval = 0;

  //
  // compute 2-e Coulomb integral first derivatives
  //

  Ref<TwoBodyDerivInt> twoecoulombderiv_inteval = integral->electron_repulsion_deriv();
  buffer = twoecoulombderiv_inteval->buffer();
  DerivCenters derivcenters;
  std::cout << "two-e deriv Coulomb integrals:" << std::endl;
  for(int s1=0; s1<nshell; s1++) {
    const int bf1_offset = obs->shell_to_function(s1);
    const int nbf1 = obs->shell(s1).nfunction();

    for(int s2=0; s2<nshell; s2++) {
      const int bf2_offset = obs->shell_to_function(s2);
      const int nbf2 = obs->shell(s2).nfunction();

      for(int s3=0; s3<nshell; s3++) {
        const int bf3_offset = obs->shell_to_function(s3);
        const int nbf3 = obs->shell(s3).nfunction();

        for(int s4=0; s4<nshell; s4++) {
          const int bf4_offset = obs->shell_to_function(s4);
          const int nbf4 = obs->shell(s4).nfunction();

          twoecoulombderiv_inteval->compute_shell(s1, s2, s3, s4, derivcenters);
          const int nderiv_centers = derivcenters.n();
          const int nderivs = nderiv_centers * 3;
          // for the sake of efficiency, the two-body derivative interface does not compute derivatives
          // that can be computed by the user using translational invariance ... see below how
          // the translational invariance can be used to recompute the "omitted" derivatives
          const bool has_omitted_center = derivcenters.has_omitted_center();

          int bf1234 = 0;
          int bf1234_deriv = 0;
          for(int bf1=0; bf1<nbf1; ++bf1) {
            for(int bf2=0; bf2<nbf2; ++bf2) {
              for(int bf3=0; bf3<nbf3; ++bf3) {
                for(int bf4=0; bf4<nbf4; ++bf4, ++bf1234, bf1234_deriv+=nderivs) {

                  int cxyz = 0;
                  for(int center=0; center<nderiv_centers; ++center) {
                    for(int xyz=0; xyz<3; ++xyz, ++cxyz) {
                      std::cout << bf1+bf1_offset << " " << bf2+bf2_offset << " "
                                << bf3+bf3_offset << " " << bf4+bf4_offset << " "
                                << " atom" << derivcenters.atom(center) << "_" << xyzstring(xyz) << " "
                                << std::setprecision(15) << buffer[bf1234_deriv + cxyz] << std::endl;
                    }
                  }

                  // if one of the centers is omitted, reconstruct its derivative from the computed ones
                  if (has_omitted_center) {
                    const int atom = derivcenters.omitted_atom();
                    for(int xyz=0; xyz<3; ++xyz, ++cxyz) {

                      // translational invariance == sum of derivatives
                      // with respect to all centers involved in the integral
                      // in the same direction (x,y,z) is zero
                      double value = 0.0;
                      for(int center=0; center<nderiv_centers; ++center)
                        value -= buffer[bf1234_deriv + center*3 + xyz];

                      std::cout << bf1+bf1_offset << " " << bf2+bf2_offset << " "
                                << bf3+bf3_offset << " " << bf4+bf4_offset << " "
                                << " atom" << atom << "_" << xyzstring(xyz) << " "
                                << std::setprecision(15) << value << std::endl;
                    }
                  }

                }
              }
            }
          }

        }
      }
    }
  }
  twoecoulombderiv_inteval = 0;
#endif

  integral = 0;

  return 0;
}
