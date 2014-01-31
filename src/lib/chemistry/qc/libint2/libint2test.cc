//
// libint2test.cc
//
// Copyright (C) 2001 Edward Valeev
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

#include <stdlib.h>
#include <string.h>

#include <util/misc/formio.h>
#include <util/misc/regtime.h>
#include <util/keyval/keyval.h>
#include <util/group/message.h>
#include <util/group/pregtime.h>
#include <math/optimize/gaussianfit.h>
#include <math/optimize/gaussianfit.timpl.h>
#include <chemistry/qc/wfn/wfn.h>
#include <chemistry/qc/basis/integral.h>
#include <chemistry/qc/basis/petite.h>
#include <chemistry/qc/basis/symmint.h>
#include <chemistry/qc/intv3/int1e.h>
#include <chemistry/qc/intv3/int2e.h>
#include <chemistry/qc/intv3/intv3.h>
#include <chemistry/qc/intv3/cartitv3.h>
#include <chemistry/qc/basis/cartit_cca.h>
#include <chemistry/qc/libint2/libint2.h>
#include <chemistry/qc/libint2/int2e.h>

#include <chemistry/qc/scf/linkage.h>
#include <chemistry/qc/libint2/linkage.h>

using namespace std;
using namespace sc;

#define GAMMA12 1.0
#define TEST_ITERATORS 1
#define TEST_1E_INTEGRALS 1
#define TEST_2E_INTEGRALS 1
#define TEST_2E_INTEGRALS_MULTIBASES 1

void
    compare_1e_libint2_vs_v3(Ref<OneBodyInt>& oblibint2, Ref<OneBodyInt>& obv3);
void compare_1e3_libint2_vs_v3(Ref<OneBodyInt>& oblibint2,
                               Ref<OneBodyInt>& obv3);
void compare_2e_libint2_vs_v3(Ref<TwoBodyInt>& tblibint2,
                              Ref<TwoBodyInt>& tbv3, bool print_all);
void compare_2e_puream_libint2_vs_v3(Ref<TwoBodyInt>& tblibint2,
                                     Ref<TwoBodyInt>& tbv3, bool print_all);
void compare_2e_bufsum_libint2_vs_v3(Ref<TwoBodyInt>& tblibint2,
                                     Ref<TwoBodyInt>& tbv3);
void compare_2e_unique_bufsum_libint2_vs_v3(Ref<TwoBodyInt>& tblibint2, Ref<
    TwoBodyInt>& tbv3);
void compare_2e(Ref<TwoBodyInt>& tb1,
                Ref<TwoBodyInt>& tb2,
                TwoBodyOper::type otype,
                bool print_all);
void print_ints(Ref<TwoBodyInt>& tblibint2, unsigned int num_te_types,
                const std::string& prefix);
void print_all_ints(Ref<TwoBodyInt>& tblibint2, unsigned int num_te_types,
                    const std::string& prefix);
void compare_2e_permute(Ref<Integral>& libint2);
void test_int_shell_1e(const Ref<KeyVal>&, const Ref<Int1eV3> &int1ev3,
                       void(Int1eV3::*int_shell_1e)(int, int), int permute);
void test_3_center(const Ref<KeyVal>&, const Ref<Int2eV3> &);
void test_4_center(const Ref<KeyVal>& keyval, const Ref<Int2eV3> &int2ev3);
void test_4der_center(const Ref<KeyVal>&, const Ref<Int2eV3> &int2ev3);
void test_p4(const Ref<KeyVal>& keyval, const Ref<Integral>& ints);

#define maxint 9

void testint(const Ref<OneBodyInt>& in) {
  if (in.null()) {
    cout << "null integral generator" << endl;
    abort();
  }
}

void testint(const Ref<OneBodyDerivInt>& in) {
  if (in.null()) {
    cout << "null integral generator" << endl;
    abort();
  }
}

void testint(const Ref<TwoBodyInt>& in) {
  if (in.null()) {
    cout << "null integral generator" << endl;
    abort();
  }
}

void testint(const Ref<TwoBodyDerivInt>& in) {
  if (in.null()) {
    cout << "null integral generator" << endl;
    abort();
  }
}

int main(int argc, char **argv) {
  int ii, i, j, k, l, m, n;

  Ref<MessageGrp> msg = MessageGrp::initial_messagegrp(argc, argv);
  if (msg.null())
    msg = new ProcMessageGrp();
  MessageGrp::set_default_messagegrp(msg);

  Ref<RegionTimer> tim = new ParallelRegionTimer(msg, "libint2test", 1, 1);

  char *infile = new char[strlen(SRCDIR) + strlen("/libint2test.in") + 1];
  sprintf(infile,SRCDIR "/libint2test.in");
  if (argc == 2) {
    delete[] infile;
    infile = argv[1];
  }

  Ref<KeyVal> pkv(new ParsedKeyVal(infile));
  Ref<KeyVal> tkeyval(new PrefixKeyVal(pkv, ":test"));

  Ref<GaussianBasisSet>
      basis =
          require_dynamic_cast<GaussianBasisSet*> (
                                                   tkeyval->describedclassvalue(
                                                                                "basis").pointer(),
                                                   "main\n");
  Ref<Molecule> mol = basis->molecule();

  int tproc = tkeyval->intvalue("test_processor");
  if (tproc >= msg->n())
    tproc = 0;
  int me = msg->me();

  if (me == tproc)
    cout << "testing on processor " << tproc << endl;

  int storage = tkeyval->intvalue("storage");
  cout << "storage = " << storage << endl;

  tim->enter("Integral");
  Ref<Integral> integral = new IntegralV3(basis);
  Ref<Integral> integrallibint2 = new IntegralLibint2(basis);

  Ref<OneBodyInt> overlapv3 = integral->overlap();
  Ref<OneBodyInt> kineticv3 = integral->kinetic();
  Ref<OneBodyInt> nuclearv3 = integral->nuclear();
  Ref<OneBodyInt> hcorev3 = integral->hcore();
  Ref<OneBodyInt> edipolev3 = integral->dipole(0);

  Ref<OneBodyInt> overlaplibint2 = integrallibint2->overlap();
  testint(overlaplibint2);
  Ref<OneBodyInt> kineticlibint2 = integrallibint2->kinetic();
  testint(kineticlibint2);
  Ref<OneBodyInt> nuclearlibint2 = integrallibint2->nuclear();
  testint(nuclearlibint2);
  Ref<OneBodyInt> hcorelibint2 = integrallibint2->hcore();
  testint(hcorelibint2);
  Ref<OneBodyInt> p4libint2 = integrallibint2->p4();
  testint(p4libint2);
  Ref<OneBodyInt> edipolelibint2 = integrallibint2->dipole(0);
  testint(edipolelibint2);

  Ref<TwoBodyInt> erepv3 = integral->electron_repulsion();

  int storage_needed = integrallibint2->storage_required_eri(basis);
  cout << scprintf("Need %d bytes to create EriLibint2\n", storage_needed);
  Ref<TwoBodyInt> ereplibint2 = integrallibint2->electron_repulsion();
  testint(ereplibint2);

  storage_needed = integrallibint2->storage_required_g12nc(basis);
  cout << scprintf("Need %d bytes to create G12NCLibint2\n",storage_needed);
  Ref<IntParamsG12> g12params;
  Ref<IntParamsG12> g12g12params;
  {
    using namespace sc::math;
    PowerExponential1D* w;
    // fit using weight exp(-0.005*r^6), which is flat to r=1, then falls slowly till r=2, then quickly decays to r=3
    w = new PowerExponential1D(0.005, 6, 0);
    typedef GaussianFit<Slater1D, PowerExponential1D> GTGFit;
    GTGFit fitter(5, *w, 0.0, 10.0, 1001);
    Slater1D stg(GAMMA12, 0, -1.0 / GAMMA12);
    g12params = new IntParamsG12( fitter(stg) );
    g12g12params = new IntParamsG12( g12params->bra(), g12params->bra() );
  }
  Ref<TwoBodyInt> g12nclibint2 = integrallibint2->g12nc<4>(g12params);
  testint(g12nclibint2);

# if LIBINT2_SUPPORT_G12
  storage_needed = integrallibint2->storage_required_g12(basis);
  cout << scprintf("Need %d bytes to create G12Libint2\n",storage_needed);
  Ref<TwoBodyInt> g12libint2 = integrallibint2->g12<4>(g12params);
  testint(g12libint2);
# endif
  tim->exit();

  // Test iterators
#if TEST_ITERATORS
  for (int l = 0; l <= 3; ++l) {
    CartesianIter* citerv[2];
    citerv[0] = integrallibint2->new_cartesian_iter(l);
    citerv[1] = integral->new_cartesian_iter(l);
    for (int i = 0; i < 2; ++i) {
      CartesianIter* citer = citerv[i];
      cout << "Cartesian L=" << l << " shell (" << (i == 0 ? "Libint2" : "V3")
          << "):" << endl;
      for (citer->start(); int(*citer); citer->next())
        cout << "nx = " << citer->a() << " ny = " << citer->b() << " nz = "
            << citer->c() << endl;
    }
    RedundantCartesianIter* rciterv[2];
    rciterv[0] = integrallibint2->new_redundant_cartesian_iter(l);
    rciterv[1] = integral->new_redundant_cartesian_iter(l);
    for (int i = 0; i < 2; ++i) {
      RedundantCartesianIter* rciter = rciterv[i];
      cout << "Redundant Cartesian L=" << l << " shell ("
          << (i == 0 ? "Libint2" : "V3") << "):" << endl;
      for (rciter->start(); int(*rciter); rciter->next())
        cout << "nx = " << rciter->a() << " ny = " << rciter->b() << " nz = "
            << rciter->c() << endl;
    }
  }
#endif

  // Test 1-e integrals
#if TEST_1E_INTEGRALS
  if (basis->has_pure())
    cout << "Cannot test 1e between Libint2 and IntV3 in presence of spherical harmonics yet"
         << endl;
  else {
    cout << "Testing Libint2' overlap integrals against IntV3's" << endl;
    compare_1e_libint2_vs_v3(overlaplibint2, overlapv3);
    cout << "Testing Libint2' kinetic energy integrals against IntV3's" << endl;
    compare_1e_libint2_vs_v3(kineticlibint2, kineticv3);
    cout << "Testing Libint2' nuclear attraction integrals against IntV3's"
        << endl;
    compare_1e_libint2_vs_v3(nuclearlibint2, nuclearv3);
    cout << "Testing Libint2' core hamiltonian integrals against IntV3's"
        << endl;
    compare_1e_libint2_vs_v3(hcorelibint2, hcorev3);
    cout << "Testing Libint2' electric dipole moment integrals against IntV3's"
        << endl;
    compare_1e3_libint2_vs_v3(edipolelibint2, edipolev3);
  }

  cout << "Testing Libint2' p^4 integrals" << endl;
  test_p4(tkeyval, integrallibint2);

#endif

  //  compare_2e_permute(integrallibint2);

#if TEST_2E_INTEGRALS
  bool puream = basis->has_pure();
  cout << "spherical harmonics " << (puream ? "" : "not")
      << " present in the basis" << endl;
  cout << "Testing Libint2' ERIs against IntV3's" << endl;
  const bool print_all = false;
  if (puream)
    compare_2e_puream_libint2_vs_v3(ereplibint2, erepv3, print_all);
  else
    compare_2e_libint2_vs_v3(ereplibint2, erepv3, print_all);

  cout << "Testing Libint2' ERIs (from G12NCLibint2) against IntV3's" << endl;
  if (puream)
    compare_2e_puream_libint2_vs_v3(g12nclibint2, erepv3, print_all);
  else
    compare_2e_libint2_vs_v3(g12nclibint2, erepv3, print_all);

  cout << "Testing Libint2' ERIs (from G12NCLibint2) against Libint2' ERIs (from ERILibint2)" << endl;
  compare_2e(g12nclibint2, ereplibint2, TwoBodyOper::eri, print_all);

  {
    cout << "Testing Libint2' G12: G12NCLibint2 vs R12_0_G12" << endl;
    {
      Ref<TwoBodyInt> tbint = integrallibint2->r12_k_g12<4>(g12params, 0);
      testint(tbint);
      compare_2e(g12nclibint2, tbint, TwoBodyOper::r12_0_g12, print_all);
    }
    cout << "Testing Libint2' G12/R12: G12NCLibint2 vs R12_m1_G12" << endl;
    {
      Ref<TwoBodyInt> tbint = integrallibint2->r12_k_g12<4>(g12params, -1);
      testint(tbint);
      compare_2e(g12nclibint2, tbint, TwoBodyOper::r12_m1_g12, print_all);
    }
    cout << "Testing Libint2' G12^2: G12NCLibint2 vs R12_0_G12" << endl;
    {
      Ref<TwoBodyInt> g12nclibint2 = integrallibint2->g12nc<4>(g12g12params);
      testint(g12nclibint2);

      Ref<TwoBodyInt> tbint = integrallibint2->r12_k_g12<4>(g12g12params, 0);
      testint(tbint);
      compare_2e(g12nclibint2, tbint, TwoBodyOper::r12_0_g12, print_all);
    }
    cout << "Testing Libint2' G12T1G12: G12NCLibint2 vs G12_T1_G12" << endl;
    {
      Ref<TwoBodyInt> g12nclibint2 = integrallibint2->g12nc<4>(g12g12params);
      testint(g12nclibint2);

      Ref<TwoBodyInt> tbint = integrallibint2->g12t1g12<4>(g12g12params);
      testint(tbint);
      compare_2e(g12nclibint2, tbint, TwoBodyOper::g12t1g12, print_all);
    }
  }

#if LIBINT2_SUPPORT_G12
  cout << "Testing Libint2' ERIs (from G12Libint2) against IntV3's" << endl;
  if (puream)
    compare_2e_puream_libint2_vs_v3(g12libint2, erepv3, print_all);
  else
    compare_2e_libint2_vs_v3(g12libint2, erepv3, print_all);
#endif

  cout << "Testing sums of Libint2' ERIs against IntV3's" << endl;
  compare_2e_bufsum_libint2_vs_v3(ereplibint2, erepv3);
# if LIBINT2_SUPPORT_G12
  cout << "Testing sums of Libint2' ERIs (from G12Libint2) against IntV3's" << endl;
  compare_2e_bufsum_libint2_vs_v3(g12libint2,erepv3);
# endif

  cout << "Testing sums of unique Libint2' ERIs against IntV3's" << endl;
  ereplibint2->set_redundant(0);
  erepv3->set_redundant(0);
  compare_2e_unique_bufsum_libint2_vs_v3(ereplibint2, erepv3);
  erepv3->set_redundant(1);
  ereplibint2->set_redundant(1);
  cout << "Printing ERI integrals" << endl;
  print_all_ints(ereplibint2, 1, "eri");
#if LIBINT2_SUPPORT_G12
  cout << "Printing G12 integrals" << endl;
  print_all_ints(g12libint2,6,"g12");
#endif

#if TEST_2E_INTEGRALS_MULTIBASES
  {
    Ref<GaussianBasisSet>
        basis1 =
            require_dynamic_cast<GaussianBasisSet*> (
                                                     tkeyval->describedclassvalue(
                                                                                  "basis1").pointer(),
                                                     "main\n");
    Ref<GaussianBasisSet>
        basis2 =
            require_dynamic_cast<GaussianBasisSet*> (
                                                     tkeyval->describedclassvalue(
                                                                                  "basis2").pointer(),
                                                     "main\n");
    Ref<GaussianBasisSet>
        basis3 =
            require_dynamic_cast<GaussianBasisSet*> (
                                                     tkeyval->describedclassvalue(
                                                                                  "basis3").pointer(),
                                                     "main\n");
    Ref<GaussianBasisSet>
        basis4 =
            require_dynamic_cast<GaussianBasisSet*> (
                                                     tkeyval->describedclassvalue(
                                                                                  "basis4").pointer(),
                                                     "main\n");

    // all bases must have same puream
    const bool can_do_this_test = (basis1->has_pure() == basis2->has_pure())
        && (basis1->has_pure() == basis3->has_pure()) && (basis1->has_pure()
        == basis4->has_pure());
    if (can_do_this_test) {
      const bool puream = basis1->has_pure();
      cout << "Testing Libint2' ERIs against IntV3's with multiple basis sets"
          << endl;
      cout << "spherical harmonics " << (puream ? "" : "not")
          << " present in the basis1" << endl;
      cout << "spherical harmonics " << (puream ? "" : "not")
          << " present in the basis2" << endl;
      cout << "spherical harmonics " << (puream ? "" : "not")
          << " present in the basis3" << endl;
      cout << "spherical harmonics " << (puream ? "" : "not")
          << " present in the basis4" << endl;

      integrallibint2->set_basis(basis1, basis2, basis3, basis4);
      Ref<TwoBodyInt> ereplibint2 = integrallibint2->electron_repulsion();
      integral->set_basis(basis1, basis2, basis3, basis4);
      Ref<TwoBodyInt> erepv3 = integral->electron_repulsion();

      const bool print_all = false;
      if (puream)
        compare_2e_puream_libint2_vs_v3(ereplibint2, erepv3, print_all);
      else {
        compare_2e_libint2_vs_v3(ereplibint2, erepv3, print_all);
      }
    }
  }
#endif

#endif // TEST_2E_INTEGRALS

  //  tim->print();
  return 0;
}

void compare_1e_libint2_vs_v3(Ref<OneBodyInt>& oblibint2, Ref<OneBodyInt>& obv3) {
  Ref<GaussianBasisSet> basis = oblibint2->basis(0);
  for (int sh1 = 4; sh1 < basis->nshell(); sh1++)
    for (int sh2 = 0; sh2 < basis->nshell(); sh2++) {
      int nbf2 = basis->shell(sh2).nfunction();
      obv3->compute_shell(sh1, sh2);
      oblibint2->compute_shell(sh1, sh2);
      const double *bufferlibint2 = oblibint2->buffer();
      const double *bufferv3 = obv3->buffer();

      int bf1_offset = 0;
      for (int gc1 = 0; gc1 < basis->shell(sh1).ncontraction(); gc1++) {
        int am1 = basis->shell(sh1).am(gc1);
        CartesianIter* citer1 = oblibint2->integral()->new_cartesian_iter(am1);
        CartesianIter* iter1 = obv3->integral()->new_cartesian_iter(am1);
        for (citer1->start(); int(*citer1); citer1->next()) {
          int bf1libint2 = bf1_offset + citer1->bfn();
          int bf1v3;
          for (iter1->start(); int(*iter1); iter1->next()) {
            if (iter1->a() == citer1->a() && iter1->b() == citer1->b()
                && iter1->c() == citer1->c()) {
              bf1v3 = bf1_offset + iter1->bfn();
              break;
            }
          }

          int bf2_offset = 0;
          for (int gc2 = 0; gc2 < basis->shell(sh2).ncontraction(); gc2++) {
            int am2 = basis->shell(sh2).am(gc2);
            CartesianIter* citer2 =
                oblibint2->integral()->new_cartesian_iter(am2);
            CartesianIter* iter2 = obv3->integral()->new_cartesian_iter(am2);

            for (citer2->start(); int(*citer2); citer2->next()) {
              int bf2libint2 = bf2_offset + citer2->bfn();
              int bf2v3;
              for (iter2->start(); int(*iter2); iter2->next()) {
                if (iter2->a() == citer2->a() && iter2->b() == citer2->b()
                    && iter2->c() == citer2->c()) {
                  bf2v3 = bf2_offset + iter2->bfn();
                  break;
                }
              }

              double valuelibint2 = bufferlibint2[bf1libint2 * nbf2
                  + bf2libint2];
              double valuev3 = bufferv3[bf1v3 * nbf2 + bf2v3];
              if (fabs(valuelibint2 - valuev3) > 1E-13) {
                cout << scprintf("Discrepancy in OEInt(sh1 = %d, sh2 = %d)\n",
                                 sh1, sh2);
                cout
                    << scprintf(
                                "bf1 = %d   bf2 = %d  OEIntegral(libint2) = %20.15lf\n",
                                bf1libint2, bf2libint2, valuelibint2);
                cout
                    << scprintf(
                                "bf1 = %d   bf2 = %d  OEIntegral(V3)    = %20.15lf\n\n",
                                bf1v3, bf2v3, valuev3);
              }
            }
            bf2_offset += basis->shell(sh2).nfunction(gc2);
          }
        }
        bf1_offset += basis->shell(sh1).nfunction(gc1);
      }
    }
}

void compare_1e3_libint2_vs_v3(Ref<OneBodyInt>& oblibint2,
                               Ref<OneBodyInt>& obv3) {
  Ref<GaussianBasisSet> basis = oblibint2->basis(0);
  for (int sh1 = 4; sh1 < basis->nshell(); sh1++)
    for (int sh2 = 0; sh2 < basis->nshell(); sh2++) {
      int nbf2 = basis->shell(sh2).nfunction();
      obv3->compute_shell(sh1, sh2);
      oblibint2->compute_shell(sh1, sh2);
      const double *bufferlibint2 = oblibint2->buffer();
      const double *bufferv3 = obv3->buffer();

      int bf1_offset = 0;
      for (int gc1 = 0; gc1 < basis->shell(sh1).ncontraction(); gc1++) {
        int am1 = basis->shell(sh1).am(gc1);
        CartesianIter* citer1 = oblibint2->integral()->new_cartesian_iter(am1);
        CartesianIter* iter1 = obv3->integral()->new_cartesian_iter(am1);
        for (citer1->start(); int(*citer1); citer1->next()) {
          int bf1libint2 = bf1_offset + citer1->bfn();
          int bf1v3;
          for (iter1->start(); int(*iter1); iter1->next()) {
            if (iter1->a() == citer1->a() && iter1->b() == citer1->b()
                && iter1->c() == citer1->c()) {
              bf1v3 = bf1_offset + iter1->bfn();
              break;
            }
          }

          int bf2_offset = 0;
          for (int gc2 = 0; gc2 < basis->shell(sh2).ncontraction(); gc2++) {
            int am2 = basis->shell(sh2).am(gc2);
            CartesianIter* citer2 =
                oblibint2->integral()->new_cartesian_iter(am2);
            CartesianIter* iter2 = obv3->integral()->new_cartesian_iter(am2);

            for (citer2->start(); int(*citer2); citer2->next()) {
              int bf2libint2 = bf2_offset + citer2->bfn();
              int bf2v3;
              for (iter2->start(); int(*iter2); iter2->next()) {
                if (iter2->a() == citer2->a() && iter2->b() == citer2->b()
                    && iter2->c() == citer2->c()) {
                  bf2v3 = bf2_offset + iter2->bfn();
                  break;
                }
              }

              for (int xyz = 0; xyz < 3; ++xyz) {
                double valuelibint2 = bufferlibint2[(bf1libint2 * nbf2
                    + bf2libint2) * 3 + xyz];
                double valuev3 = bufferv3[(bf1v3 * nbf2 + bf2v3) * 3
                    + xyz];
                if (fabs(valuelibint2 - valuev3) > 1E-13) {
                  cout
                      << scprintf("Discrepancy in OEInt(sh1 = %d, sh2 = %d)\n",
                                  sh1, sh2);
                  cout
                      << scprintf(
                                  "bf1 = %d   bf2 = %d   xyz = %d   OEIntegral(libint2) = %20.15lf\n",
                                  bf1libint2, bf2libint2, xyz, valuelibint2);
                  cout
                      << scprintf(
                                  "bf1 = %d   bf2 = %d   xyz = %d   OEIntegral(V3)    = %20.15lf\n\n",
                                  bf1v3, bf2v3, xyz, valuev3);
                }
              }
            }

            bf2_offset += basis->shell(sh2).nfunction(gc2);
          }
        }
        bf1_offset += basis->shell(sh1).nfunction(gc1);
      }
    }
}

void compare_2e_libint2_vs_v3(Ref<TwoBodyInt>& tblibint2,
                              Ref<TwoBodyInt>& tbv3, bool print_all) {
#if INTEGRALLIBINT2_NORMCONV != INTEGRALLIBINT2_NORMCONV_CCA
  cout << "In IntegralLibint2 normalization convention != CCA, cannot compare Cartesian integrals vs IntegralV3" << endl;
  return;
#endif
  const double *bufferlibint2 = tblibint2->buffer();
  const double *bufferv3 = tbv3->buffer();

  Ref<GaussianBasisSet> basis1 = tblibint2->basis1();
  Ref<GaussianBasisSet> basis2 = tblibint2->basis2();
  Ref<GaussianBasisSet> basis3 = tblibint2->basis3();
  Ref<GaussianBasisSet> basis4 = tblibint2->basis4();
  for (int sh1 = 0; sh1 < basis1->nshell(); sh1++)
    for (int sh2 = 0; sh2 < basis2->nshell(); sh2++)
      for (int sh3 = 0; sh3 < basis3->nshell(); sh3++)
        for (int sh4 = 0; sh4 < basis4->nshell(); sh4++) {
#if 0
          cout << scprintf("Computing TEInt(sh1 = %d, sh2 = %d, sh3 = %d, sh4 = %d)\n",sh1,sh2,sh3,sh4);
#endif
          tbv3->compute_shell(sh1, sh2, sh3, sh4);
          tblibint2->compute_shell(sh1, sh2, sh3, sh4);

          int nbf2 = basis2->shell(sh2).nfunction();
          int nbf3 = basis3->shell(sh3).nfunction();
          int nbf4 = basis4->shell(sh4).nfunction();

          int bf1_offset = 0;
          for (int gc1 = 0; gc1 < basis1->shell(sh1).ncontraction(); gc1++) {
            int am1 = basis1->shell(sh1).am(gc1);
            CartesianIter* citer1 =
                tblibint2->integral()->new_cartesian_iter(am1);
            CartesianIter* iter1 = tbv3->integral()->new_cartesian_iter(am1);
            for (citer1->start(); int(*citer1); citer1->next()) {
              int bf1libint2 = citer1->bfn();
              int bf1v3;
              for (iter1->start(); int(*iter1); iter1->next()) {
                if (iter1->a() == citer1->a() && iter1->b() == citer1->b()
                    && iter1->c() == citer1->c()) {
                  bf1v3 = iter1->bfn();
                  break;
                }
              }
              bf1libint2 += bf1_offset;
              bf1v3 += bf1_offset;

              int bf2_offset = 0;
              for (int gc2 = 0; gc2 < basis2->shell(sh2).ncontraction(); gc2++) {
                int am2 = basis2->shell(sh2).am(gc2);
                CartesianIter* citer2 =
                    tblibint2->integral()->new_cartesian_iter(am2);
                CartesianIter* iter2 =
                    tbv3->integral()->new_cartesian_iter(am2);
                for (citer2->start(); int(*citer2); citer2->next()) {
                  int bf2libint2 = citer2->bfn();
                  int bf2v3;
                  for (iter2->start(); int(*iter2); iter2->next()) {
                    if (iter2->a() == citer2->a() && iter2->b() == citer2->b()
                        && iter2->c() == citer2->c()) {
                      bf2v3 = iter2->bfn();
                      break;
                    }
                  }
                  bf2libint2 += bf2_offset;
                  bf2v3 += bf2_offset;

                  int bf3_offset = 0;
                  for (int gc3 = 0; gc3 < basis3->shell(sh3).ncontraction(); gc3++) {
                    int am3 = basis3->shell(sh3).am(gc3);
                    CartesianIter* citer3 =
                        tblibint2->integral()->new_cartesian_iter(am3);
                    CartesianIter* iter3 =
                        tbv3->integral()->new_cartesian_iter(am3);
                    for (citer3->start(); int(*citer3); citer3->next()) {
                      int bf3libint2 = citer3->bfn();
                      int bf3v3;
                      for (iter3->start(); int(*iter3); iter3->next()) {
                        if (iter3->a() == citer3->a() && iter3->b()
                            == citer3->b() && iter3->c() == citer3->c()) {
                          bf3v3 = iter3->bfn();
                          break;
                        }
                      }
                      bf3libint2 += bf3_offset;
                      bf3v3 += bf3_offset;

                      int bf4_offset = 0;
                      for (int gc4 = 0; gc4 < basis4->shell(sh4).ncontraction(); gc4++) {
                        int am4 = basis4->shell(sh4).am(gc4);
                        CartesianIter* citer4 =
                            tblibint2->integral()->new_cartesian_iter(am4);
                        CartesianIter* iter4 =
                            tbv3->integral()->new_cartesian_iter(am4);
                        for (citer4->start(); int(*citer4); citer4->next()) {
                          int bf4libint2 = citer4->bfn();
                          int bf4v3;
                          for (iter4->start(); int(*iter4); iter4->next()) {
                            if (iter4->a() == citer4->a() && iter4->b()
                                == citer4->b() && iter4->c() == citer4->c()) {
                              bf4v3 = iter4->bfn();
                              break;
                            }
                          }
                          bf4libint2 += bf4_offset;
                          bf4v3 += bf4_offset;

                          double valuelibint2 = bufferlibint2[((bf1libint2
                              * nbf2 + bf2libint2) * nbf3 + bf3libint2) * nbf4
                              + bf4libint2];
                          double valuev3 = bufferv3[((bf1v3 * nbf2 + bf2v3)
                              * nbf3 + bf3v3) * nbf4 + bf4v3];
                          const bool significant_discrepancy =
                              (fabs(valuelibint2 - valuev3) > 1E-12);
                          if (significant_discrepancy || print_all) {
                            if (significant_discrepancy)
                              cout << scprintf(
                                              "Discrepancy in TEInt(sh1 = %d, sh2 = %d, sh3 = %d, sh4 = %d)\n",
                                              sh1, sh2, sh3, sh4);
                            else
                              cout << scprintf(
                                              "Debugging TEInt(sh1 = %d, sh2 = %d, sh3 = %d, sh4 = %d)\n",
                                              sh1, sh2, sh3, sh4);

                            cout << scprintf(
                                            "bf1 = %d  bf2 = %d  bf3 = %d  bf4 = %d  TEIntegral(libint2) = %20.15lf\n",
                                            bf1libint2, bf2libint2, bf3libint2,
                                            bf4libint2, valuelibint2);
                            cout << scprintf(
                                            "bf1 = %d  bf2 = %d  bf3 = %d  bf4 = %d  TEIntegral(V3)    = %20.15lf\n\n",
                                            bf1v3, bf2v3, bf3v3, bf4v3, valuev3);
                          }
                        }
                        bf4_offset += basis4->shell(sh4).nfunction(gc4);
                      }
                    }
                    bf3_offset += basis3->shell(sh3).nfunction(gc3);
                  }
                }
                bf2_offset += basis2->shell(sh2).nfunction(gc2);
              }
            }
            bf1_offset += basis1->shell(sh1).nfunction(gc1);
          }
          //return;
        }
}

void compare_2e_puream_libint2_vs_v3(Ref<TwoBodyInt>& tblibint2,
                                     Ref<TwoBodyInt>& tbv3, bool print_all) {
  const double *bufferlibint2 = tblibint2->buffer();
  const double *bufferv3 = tbv3->buffer();

  Ref<GaussianBasisSet> basis1 = tblibint2->basis1();
  Ref<GaussianBasisSet> basis2 = tblibint2->basis2();
  Ref<GaussianBasisSet> basis3 = tblibint2->basis3();
  Ref<GaussianBasisSet> basis4 = tblibint2->basis4();
  for (int sh1 = 0; sh1 < basis1->nshell(); sh1++)
    for (int sh2 = 0; sh2 < basis2->nshell(); sh2++)
      for (int sh3 = 0; sh3 < basis3->nshell(); sh3++)
        for (int sh4 = 0; sh4 < basis4->nshell(); sh4++) {
          //	    	    	    sh1 = 0; sh2 = 0; sh3 = 6; sh4 = 13;
          /*	    if ( !((basis->shell(sh1).has_pure() || basis->shell(sh1).max_am()==0) &&
           (basis->shell(sh2).has_pure() || basis->shell(sh2).max_am()==0) &&
           (basis->shell(sh3).has_pure() || basis->shell(sh3).max_am()==0) &&
           (basis->shell(sh4).has_pure() || basis->shell(sh4).max_am()==0)) )
           continue;*/
          tblibint2->compute_shell(sh1, sh2, sh3, sh4);
          tbv3->compute_shell(sh1, sh2, sh3, sh4);

          int nbf1 = basis1->shell(sh1).nfunction();
          int nbf2 = basis2->shell(sh2).nfunction();
          int nbf3 = basis3->shell(sh3).nfunction();
          int nbf4 = basis4->shell(sh4).nfunction();

#if 0
          cout << scprintf(
                           "Computing TEInt(sh1 = %d, sh2 = %d, sh3 = %d, sh4 = %d)\n",
                           sh1, sh2, sh3, sh4);
#endif
          int i1234 = 0;
          for(int bf1 = 0; bf1 < nbf1; ++bf1) {
            for(int bf2 = 0; bf2 < nbf2; ++bf2) {
              for(int bf3 = 0; bf3 < nbf3; ++bf3) {
                for(int bf4 = 0; bf4 < nbf4; ++bf4, ++i1234) {

                  double valuelibint2 = bufferlibint2[i1234];
                  double valuev3 = bufferv3[i1234];
                  const bool significant_discrepancy = (fabs(valuelibint2 - valuev3) > 1E-12);
                  if (significant_discrepancy || print_all) {
                    if (significant_discrepancy)
                      cout << scprintf(
                                       "Discrepancy in TEInt(sh1 = %d, sh2 = %d, sh3 = %d, sh4 = %d)\n",
                                       sh1, sh2, sh3, sh4);
                    else
                      cout << scprintf(
                                       "Debugging TEInt(sh1 = %d, sh2 = %d, sh3 = %d, sh4 = %d)\n",
                                       sh1, sh2, sh3, sh4);

                    cout << scprintf(
                                    "bf1 = %d  bf2 = %d  bf3 = %d  bf4 = %d  TEIntegral(libint2) = %20.15lf\n",
                                    bf1, bf2, bf3, bf4, valuelibint2);
                    cout << scprintf(
                                    "bf1 = %d  bf2 = %d  bf3 = %d  bf4 = %d  TEIntegral(V3)    = %20.15lf\n\n",
                                    bf1, bf2, bf3, bf4, valuev3);
                  }
                }
              }
            }
          }
          //  	    	    return;
        }
}

void compare_2e(Ref<TwoBodyInt>& tb1,
                Ref<TwoBodyInt>& tb2,
                TwoBodyOper::type otype,
                bool print_all) {
  const double *buf1 = tb1->buffer(otype); MPQC_ASSERT(buf1 != 0);
  const double *buf2 = tb2->buffer(otype); MPQC_ASSERT(buf2 != 0);

  Ref<GaussianBasisSet> basis1 = tb1->basis1(); MPQC_ASSERT(basis1 == tb2->basis1());
  Ref<GaussianBasisSet> basis2 = tb1->basis2(); MPQC_ASSERT(basis2 == tb2->basis2());
  Ref<GaussianBasisSet> basis3 = tb1->basis3(); MPQC_ASSERT(basis3 == tb2->basis3());
  Ref<GaussianBasisSet> basis4 = tb1->basis4(); MPQC_ASSERT(basis4 == tb2->basis4());
  for (int sh1 = 0; sh1 < basis1->nshell(); sh1++)
    for (int sh2 = 0; sh2 < basis2->nshell(); sh2++)
      for (int sh3 = 0; sh3 < basis3->nshell(); sh3++)
        for (int sh4 = 0; sh4 < basis4->nshell(); sh4++) {
          //                        sh1 = 0; sh2 = 0; sh3 = 6; sh4 = 13;
          /*        if ( !((basis->shell(sh1).has_pure() || basis->shell(sh1).max_am()==0) &&
           (basis->shell(sh2).has_pure() || basis->shell(sh2).max_am()==0) &&
           (basis->shell(sh3).has_pure() || basis->shell(sh3).max_am()==0) &&
           (basis->shell(sh4).has_pure() || basis->shell(sh4).max_am()==0)) )
           continue;*/
          tb1->compute_shell(sh1, sh2, sh3, sh4);
          tb2->compute_shell(sh1, sh2, sh3, sh4);

          int nbf1 = basis1->shell(sh1).nfunction();
          int nbf2 = basis2->shell(sh2).nfunction();
          int nbf3 = basis3->shell(sh3).nfunction();
          int nbf4 = basis4->shell(sh4).nfunction();

#if 0
          cout << scprintf(
                           "Computing TEInt(sh1 = %d, sh2 = %d, sh3 = %d, sh4 = %d)\n",
                           sh1, sh2, sh3, sh4);
#endif
          int i1234 = 0;
          for(int bf1 = 0; bf1 < nbf1; ++bf1) {
            for(int bf2 = 0; bf2 < nbf2; ++bf2) {
              for(int bf3 = 0; bf3 < nbf3; ++bf3) {
                for(int bf4 = 0; bf4 < nbf4; ++bf4, ++i1234) {

                  double value1 = buf1[i1234];
                  double value2 = buf2[i1234];
                  const bool significant_discrepancy = (fabs(value1 - value2) > 1E-12);
                  if (significant_discrepancy || print_all) {
                    if (significant_discrepancy)
                      cout << scprintf(
                                       "Discrepancy in TEInt(sh1 = %d, sh2 = %d, sh3 = %d, sh4 = %d)\n",
                                       sh1, sh2, sh3, sh4);
                    else
                      cout << scprintf(
                                       "Debugging TEInt(sh1 = %d, sh2 = %d, sh3 = %d, sh4 = %d)\n",
                                       sh1, sh2, sh3, sh4);

                    cout << scprintf(
                                    "bf1 = %d  bf2 = %d  bf3 = %d  bf4 = %d  TEIntegral(1) = %20.15lf\n",
                                    bf1, bf2, bf3, bf4, value1);
                    cout << scprintf(
                                    "bf1 = %d  bf2 = %d  bf3 = %d  bf4 = %d  TEIntegral(2)    = %20.15lf\n\n",
                                    bf1, bf2, bf3, bf4, value2);
                  }
                }
              }
            }
          }
          //                return;
        }
}

void compare_2e_bufsum_libint2_vs_v3(Ref<TwoBodyInt>& tblibint2,
                                     Ref<TwoBodyInt>& tbv3) {
  Ref<GaussianBasisSet> basis = tblibint2->basis();
  const double *bufferlibint2 = tblibint2->buffer();
  const double *bufferv3 = tbv3->buffer();

  for (int sh1 = 0; sh1 < basis->nshell(); sh1++)
    for (int sh2 = 0; sh2 < basis->nshell(); sh2++)
      for (int sh3 = 0; sh3 < basis->nshell(); sh3++)
        for (int sh4 = 0; sh4 < basis->nshell(); sh4++) {
          //	    sh1=12;sh2=12;sh3=12;sh4=12;
          tblibint2->compute_shell(sh1, sh2, sh3, sh4);
          tbv3->compute_shell(sh1, sh2, sh3, sh4);

          int nbf1 = basis->shell(sh1).nfunction();
          int nbf2 = basis->shell(sh2).nfunction();
          int nbf3 = basis->shell(sh3).nfunction();
          int nbf4 = basis->shell(sh4).nfunction();

          double sum_libint2 = 0.0;
          double sum_v3 = 0.0;

          int index = 0;
          for (int i = 0; i < nbf1; i++) {
            for (int j = 0; j < nbf2; j++) {
              for (int k = 0; k < nbf3; k++) {
                for (int l = 0; l < nbf4; l++) {
                  sum_libint2 += bufferlibint2[index];
                  sum_v3 += bufferv3[index];
                  /*		    cout << scprintf("index = %d TEIntegral(libint2) = %20.15lf\n",
                   index,bufferlibint2[index]);
                   cout << scprintf("index = %d TEIntegral(V3) = %20.15lf\n\n",
                   index,bufferv3[index]);
                   */
                  index++;
                }
              }
            }
          }

          if (fabs(sum_libint2 - sum_v3) > 1E-10) {
            cout
                << scprintf(
                            "Discrepancy in TEInt(sh1 = %d, sh2 = %d, sh3 = %d, sh4 = %d)\n",
                            sh1, sh2, sh3, sh4);
            cout
                << scprintf("TEIntegralSum(libint2) = %20.15lf\n", sum_libint2);
            cout << scprintf("TEIntegralSum(V3)    = %20.15lf\n\n", sum_v3);
          }
          //return;
        }
}

void compare_2e_unique_bufsum_libint2_vs_v3(Ref<TwoBodyInt>& tblibint2, Ref<
    TwoBodyInt>& tbv3) {
  Ref<GaussianBasisSet> basis = tblibint2->basis();
  const double *bufferlibint2 = tblibint2->buffer();
  const double *bufferv3 = tbv3->buffer();

  for (int sh1 = 0; sh1 < basis->nshell(); sh1++)
    for (int sh2 = 0; sh2 <= sh1; sh2++)
      for (int sh3 = 0; sh3 <= sh1; sh3++)
        for (int sh4 = 0; sh4 <= ((sh1 == sh3) ? sh2 : sh3); sh4++) {
          tblibint2->compute_shell(sh1, sh2, sh3, sh4);
          tbv3->compute_shell(sh1, sh2, sh3, sh4);

          int nbf1 = basis->shell(sh1).nfunction();
          int nbf2 = basis->shell(sh2).nfunction();
          int nbf3 = basis->shell(sh3).nfunction();
          int nbf4 = basis->shell(sh4).nfunction();

          double sum_libint2 = 0.0;
          double sum_v3 = 0.0;

          int e12 = (sh1 == sh2) ? 1 : 0;
          int e34 = (sh3 == sh4) ? 1 : 0;
          int e13e24 = (((sh1 == sh3) && (sh2 == sh4)) || ((sh1 == sh4) && (sh2
              == sh3))) ? 1 : 0;

          int index = 0;
          for (int i = 0; i < nbf1; i++) {
            int jmax = e12 ? i : nbf2 - 1;
            for (int j = 0; j <= jmax; j++) {
              int kmax = e13e24 ? i : nbf3 - 1;
              for (int k = 0; k <= kmax; k++) {
                int lmax = e34 ? ((e13e24 && (i == k)) ? j : k) : ((e13e24
                    && (i == k)) ? j : nbf4 - 1);
                for (int l = 0; l <= lmax; l++) {
                  sum_libint2 += bufferlibint2[index];
                  sum_v3 += bufferv3[index];
                  /*		    cout << scprintf("index = %d TEIntegral(libint2) = %20.15lf\n",
                   index,bufferlibint2[index]);
                   cout << scprintf("index = %d TEIntegral(V3) = %20.15lf\n\n",
                   index,bufferv3[index]);
                   */
                  index++;
                }
              }
            }
          }

          if (fabs(sum_libint2 - sum_v3) > 1E-10) {
            cout
                << scprintf(
                            "Discrepancy in TEInt(sh1 = %d, sh2 = %d, sh3 = %d, sh4 = %d)\n",
                            sh1, sh2, sh3, sh4);
            cout
                << scprintf("TEIntegralSum(libint2) = %20.15lf\n", sum_libint2);
            cout << scprintf("TEIntegralSum(V3)    = %20.15lf\n\n", sum_v3);
          }
        }
}

void print_ints(Ref<TwoBodyInt>& tblibint2, unsigned int num_te_types,
                const std::string& prefix) {
  Ref<GaussianBasisSet> basis = tblibint2->basis();
  const double *buffer[num_te_types];
  for (int te_type = 0; te_type < num_te_types; te_type++) {
    buffer[te_type]
        = tblibint2->buffer(static_cast<TwoBodyOper::type> (te_type));
  }
  char teout_filename[] = "teout0.dat";
  FILE *teout[num_te_types];

  for (int te_type = 0; te_type < num_te_types; te_type++) {
    teout_filename[5] = te_type + '0';
    std::string fname = prefix + teout_filename;
    teout[te_type] = fopen(fname.c_str(), "w");
  }

  for (int ush1 = 0; ush1 < basis->nshell(); ush1++)
    for (int ush2 = 0; ush2 <= ush1; ush2++)
      for (int ush3 = 0; ush3 <= ush2; ush3++)
        for (int ush4 = 0; ush4 <= ush3; ush4++) {
          int S1[3], S2[3], S3[3], S4[4];
          int num = 1;
          S1[0] = ush1;
          S2[0] = ush2;
          S3[0] = ush3;
          S4[0] = ush4;
          if (ush1 == ush2 && ush1 == ush3 || ush2 == ush3 && ush2 == ush4)
            num = 1;
          else if (ush1 == ush3 || ush2 == ush4) {
            num = 2;
            S1[1] = ush1;
            S2[1] = ush3;
            S3[1] = ush2;
            S4[1] = ush4;
          } else if (ush2 == ush3) {
            num = 2;
            S1[1] = ush1;
            S2[1] = ush4;
            S3[1] = ush2;
            S4[1] = ush3;
          } else if (ush1 == ush2 || ush3 == ush4) {
            num = 2;
            S1[1] = ush1;
            S2[1] = ush3;
            S3[1] = ush2;
            S4[1] = ush4;
          } else {
            num = 3;
            S1[1] = ush1;
            S2[1] = ush3;
            S3[1] = ush2;
            S4[1] = ush4;
            S1[2] = ush1;
            S2[2] = ush4;
            S3[2] = ush2;
            S4[2] = ush3;
          }

          for (int uq = 0; uq < num; uq++) {
            int sh1 = S1[uq];
            int sh2 = S2[uq];
            int sh3 = S3[uq];
            int sh4 = S4[uq];

            tblibint2->compute_shell(sh1, sh2, sh3, sh4);

            int nbf1 = basis->shell(sh1).nfunction();
            int nbf2 = basis->shell(sh2).nfunction();
            int nbf3 = basis->shell(sh3).nfunction();
            int nbf4 = basis->shell(sh4).nfunction();

            int e12 = (sh1 == sh2) ? 1 : 0;
            int e34 = (sh3 == sh4) ? 1 : 0;
            int e13e24 = (((sh1 == sh3) && (sh2 == sh4)) || ((sh1 == sh4)
                && (sh2 == sh3))) ? 1 : 0;

            for (int i = 0; i < nbf1; i++) {
              int jmax = e12 ? i : nbf2 - 1;
              for (int j = 0; j <= jmax; j++) {
                int kmax = e13e24 ? i : nbf3 - 1;
                for (int k = 0; k <= kmax; k++) {
                  int lmax = e34 ? ((e13e24 && (i == k)) ? j : k) : ((e13e24
                      && (i == k)) ? j : nbf4 - 1);
                  for (int l = 0; l <= lmax; l++) {

                    int index = (((i * nbf2 + j) * nbf3 + k) * nbf4 + l);

                    for (int te_type = 0; te_type < num_te_types; te_type++) {
                      double integral = buffer[te_type][index];
                      if (fabs(integral) > 1E-15) {
                        fprintf(teout[te_type], "%5d%5d%5d%5d%20.10lf\n",
                                basis->shell_to_function(sh1) + i + 1,
                                basis->shell_to_function(sh2) + j + 1,
                                basis->shell_to_function(sh3) + k + 1,
                                basis->shell_to_function(sh4) + l + 1, integral);
                      }
                    }
                  }
                }
              }
            }
          }
        }
  for (int te_type = 0; te_type < num_te_types; te_type++)
    fclose(teout[te_type]);
}

void print_all_ints(Ref<TwoBodyInt>& tblibint2, unsigned int num_te_types,
                    const std::string& prefix) {
  Ref<GaussianBasisSet> basis = tblibint2->basis();
  const double *buffer[num_te_types];
  for (int te_type = 0; te_type < num_te_types; te_type++) {
    buffer[te_type]
        = tblibint2->buffer(static_cast<TwoBodyOper::type> (te_type));
  }
  char teout_filename[] = "teout0.dat";
  FILE *teout[num_te_types];

  for (int te_type = 0; te_type < num_te_types; te_type++) {
    teout_filename[5] = te_type + '0';
    std::string fname = prefix + teout_filename;
    teout[te_type] = fopen(fname.c_str(), "w");
  }

  for (int sh1 = 0; sh1 < basis->nshell(); sh1++)
    for (int sh2 = 0; sh2 < basis->nshell(); sh2++)
      for (int sh3 = 0; sh3 < basis->nshell(); sh3++)
        for (int sh4 = 0; sh4 < basis->nshell(); sh4++) {

          tblibint2->compute_shell(sh1, sh2, sh3, sh4);

          int nbf1 = basis->shell(sh1).nfunction();
          int nbf2 = basis->shell(sh2).nfunction();
          int nbf3 = basis->shell(sh3).nfunction();
          int nbf4 = basis->shell(sh4).nfunction();

          for (int i = 0; i < nbf1; i++) {
            int jmax = nbf2 - 1;
            for (int j = 0; j <= jmax; j++) {
              int kmax = nbf3 - 1;
              for (int k = 0; k <= kmax; k++) {
                int lmax = nbf4 - 1;
                for (int l = 0; l <= lmax; l++) {

                  int index = (((i * nbf2 + j) * nbf3 + k) * nbf4 + l);

                  for (int te_type = 0; te_type < num_te_types; te_type++) {
                    double integral = buffer[te_type][index];
                    if (fabs(integral) > 1E-15) {
                      fprintf(teout[te_type], "%5d%5d%5d%5d%20.10lf\n",
                              basis->shell_to_function(sh1) + i + 1,
                              basis->shell_to_function(sh2) + j + 1,
                              basis->shell_to_function(sh3) + k + 1,
                              basis->shell_to_function(sh4) + l + 1, integral);
                    }
                  }
                }
              }
            }
          }
        }

  for (int te_type = 0; te_type < num_te_types; te_type++)
    fclose(teout[te_type]);
}

void compare_2e_permute(Ref<Integral>& libint2) {
  Ref<TwoBodyInt> tb1 = libint2->electron_repulsion();
  Ref<TwoBodyInt> tb2 = libint2->electron_repulsion();
  Ref<GaussianBasisSet> basis = tb1->basis();
  const double *buffer1 = tb1->buffer();
  const double *buffer2 = tb2->buffer();

  int sh1 = 0;
  int sh2 = 0;
  int sh3 = 4;
  int sh4 = 0;

  tb1->compute_shell(sh1, sh2, sh3, sh4);
  tb2->compute_shell(sh1, sh2, sh4, sh3);

  int nbf1 = basis->shell(sh1).nfunction();
  int nbf2 = basis->shell(sh2).nfunction();
  int nbf3 = basis->shell(sh3).nfunction();
  int nbf4 = basis->shell(sh4).nfunction();

  for (int index = 0; index < nbf1 * nbf2 * nbf3 * nbf4; index++)
    if (fabs(buffer1[index] - buffer2[index]) > 1E-13) {
      cout
          << scprintf(
                      "Discrepancy in TEInt(sh1 = %d, sh2 = %d, sh3 = %d, sh4 = %d)\n",
                      sh1, sh2, sh3, sh4);
      cout << scprintf("TEIntegral(libint21)    = %20.15lf\n", buffer1[index]);
      cout
          << scprintf("TEIntegral(libint22)    = %20.15lf\n\n", buffer2[index]);
    }
}

void do_shell_test_1e(const Ref<Int1eV3> &int1ev3,
                      void(Int1eV3::*int_shell_1e)(int, int), int permute,
                      int i, int j, int na, int nb, double *buf, double *pbuf) {
  int ii = 0;
  int a;
  double *buffer = int1ev3->buffer();
  (int1ev3 ->* int_shell_1e)(i, j);
  for (a = 0; a < na * nb; a++) {
    buf[a] = buffer[a];
  }
  (int1ev3 ->* int_shell_1e)(j, i);
  for (a = 0; a < na * nb; a++) {
    pbuf[a] = buffer[a];
  }
  for (a = 0; a < na; a++) {
    for (int b = 0; b < nb; b++) {
      if (fabs(buf[ii] - pbuf[a + na * b]) > 1.0e-13) {
        cout << scprintf("----- 1e perm failed:"
          "<%d %d|%d %d>:"
          " %18.14f != %18.14f "
          "<%d %d|%d %d>\n", i, a, j, b, buf[ii], pbuf[a + na * b], j, b, i, a);
      }
      if (fabs(buf[ii]) > 1.0e-15) {
        cout << scprintf(" <(%d %d)|(%d %d)> = %15.11f\n", i, a, j, b, buf[ii]);
      }
      ii++;
    }
  }
}

void test_int_shell_1e(const Ref<KeyVal>& keyval, const Ref<Int1eV3> &int1ev3,
                       void(Int1eV3::*int_shell_1e)(int, int), int permute) {
  int flags = 0;
  Ref<GaussianBasisSet> basis = int1ev3->basis();
  int maxfunc = basis->max_nfunction_in_shell();
  int size = maxfunc * maxfunc;
  double *buf = new double[size];
  double *pbuf = new double[size];
  int nshell = int1ev3->basis()->nshell();

  for (int i = 0; i < nshell; i++) {
    int na = basis->shell(i).nfunction();
    for (int j = 0; j < nshell; j++) {
      int nb = basis->shell(j).nfunction();
      do_shell_test_1e(int1ev3, int_shell_1e, permute, i, j, na, nb, buf, pbuf);

    }
  }

  delete[] buf;
  delete[] pbuf;
}

void test_3_center(const Ref<KeyVal>& keyval, const Ref<Int2eV3> &int2ev3) {
  int ii, i, j, k, l, m, n;

  int2ev3->set_redundant(1);
  int2ev3->set_permute(0);
  double *buffer = int2ev3->buffer();
  int nshell = int2ev3->basis()->nshell();

  for (i = 0; i < nshell; i++) {
    for (j = 0; j < nshell; j++) {
      int sh[2], sizes[2];
      sh[0] = i;
      sh[1] = j;
      int2ev3->erep_2center(sh, sizes);
      ii = 0;
      for (k = 0; k < sizes[0]; k++) {
        for (l = 0; l < sizes[1]; l++) {
          if (fabs(buffer[ii]) > 1.0e-15)
            cout << scprintf(" ((%d %d)|(%d %d)) = %15.11f\n", sh[0], k, sh[1],
                             l, buffer[ii]);
          ii++;
        }
      }
    }
  }

  for (i = 0; i < nshell; i++) {
    for (j = 0; j < nshell; j++) {
      for (m = 0; m < nshell; m++) {
        int sh[3], sizes[3];
        sh[0] = i;
        sh[1] = j;
        sh[2] = m;
        int2ev3->erep_3center(sh, sizes);
        ii = 0;
        for (k = 0; k < sizes[0]; k++) {
          for (l = 0; l < sizes[1]; l++) {
            for (n = 0; n < sizes[2]; n++) {
              if (fabs(buffer[ii]) > 1.0e-15)
                cout << scprintf(" ((%d %d)|(%d %d)(%d %d)) = %15.11f\n",
                                 sh[0], k, sh[1], l, sh[2], n, buffer[ii]);
              ii++;
            }
          }
        }
      }
    }
  }

}

void init_shell_perm(const Ref<Int2eV3> &int2ev3, double *integrals,
                     double buff[maxint][maxint][maxint][maxint], int sh[4],
                     int sizes[4]) {
  int i, j, k, l;
  int oldp = int2ev3->permute();
  int2ev3->set_permute(0);
  int2ev3->erep(sh, sizes);
  int2ev3->set_permute(oldp);
  for (i = 0; i < sizes[0]; i++) {
    for (j = 0; j < sizes[1]; j++) {
      for (k = 0; k < sizes[2]; k++) {
        for (l = 0; l < sizes[3]; l++) {
          buff[i][j][k][l] = *integrals++;
        }
      }
    }
  }
}

void check_shell_perm(const Ref<Int2eV3> &int2ev3, double *integrals,
                      double buff[maxint][maxint][maxint][maxint], int sh[4],
                      int sizes[4], int p0, int p1, int p2, int p3) {
  int ip[4], p[4];
  int psizes[4];
  int psh[4];
  int index = 0;
  int i[4];
  p[0] = p0;
  p[1] = p1;
  p[2] = p2;
  p[3] = p3;
  ip[p0] = 0;
  ip[p1] = 1;
  ip[p2] = 2;
  ip[p3] = 3;
  psh[0] = sh[p0];
  psh[1] = sh[p1];
  psh[2] = sh[p2];
  psh[3] = sh[p3];
  int oldp = int2ev3->permute();
  int2ev3->set_permute(0);
  int2ev3->erep(psh, psizes);
  int2ev3->set_permute(oldp);
  for (i[0] = 0; i[0] < psizes[0]; i[0]++) {
    for (i[1] = 0; i[1] < psizes[1]; i[1]++) {
      for (i[2] = 0; i[2] < psizes[2]; i[2]++) {
        for (i[3] = 0; i[3] < psizes[3]; i[3]++) {
          if (fabs(buff[i[ip[0]]][i[ip[1]]][i[ip[2]]][i[ip[3]]]
              - integrals[index]) > 1.0e-13) {
            cout << scprintf("perm %d %d %d %d failed:"
              "((%d %d)(%d %d)|(%d %d)(%d %d)):"
              " %18.14f != %18.14f "
              "((%d %d)(%d %d)|(%d %d)(%d %d))\n", p0, p1, p2, p3, sh[0], i[0],
                             sh[1], i[1], sh[2], i[2], sh[3], i[3],
                             buff[i[ip[0]]][i[ip[1]]][i[ip[2]]][i[ip[3]]],
                             integrals[index], psh[0], i[p[0]], psh[1],
                             i[p[1]], psh[2], i[p[2]], psh[3], i[p[3]]);
          }
          index++;
        }
      }
    }
  }
}

void do_shell_quartet_test(const Ref<Int2eV3> &int2ev3, int print,
                           int printbounds, int bounds, int permute, const Ref<
                               KeyVal>& keyval, int i, int j, int k, int l) {
  int sh[4], sizes[4];
  int ibuf;
  int ii, jj, kk, ll;
  sh[0] = i;
  sh[1] = j;
  sh[2] = k;
  sh[3] = l;
  double maxintegral, integralbound;
  int boundijkl;
  if (bounds) {
    integralbound = int2ev3->logbound_to_bound((boundijkl
        = int2ev3->erep_4bound(i, j, k, l)));
  }
  double *buffer = int2ev3->buffer();
  int2ev3->erep(sh, sizes);
  ibuf = 0;
  maxintegral = 0.0;
  for (ii = 0; ii < sizes[0]; ii++) {
    for (jj = 0; jj < sizes[1]; jj++) {
      for (kk = 0; kk < sizes[2]; kk++) {
        for (ll = 0; ll < sizes[3]; ll++) {
          double absint = fabs(buffer[ibuf]);
          if (absint > maxintegral) {
            maxintegral = absint;
          }
          if (bounds && absint > integralbound) {
            cout << scprintf("((%d %d)(%d %d)|(%d %d)(%d %d)) = %15.11f, "
              "bound = %15.11f\n", sh[0], ii, sh[1], jj, sh[2], kk, sh[3], ll,
                             buffer[ibuf], integralbound);
            abort();
          }
          if (print && (absint > 1.0e-9 || (bounds && integralbound > 1.0e-9))) {
            cout << scprintf(" ((%d %d)(%d %d)|(%d %d)(%d %d))"
              " = %15.11f", sh[0], ii, sh[1], jj, sh[2], kk, sh[3], ll,
                             buffer[ibuf]);
            if (bounds) {
              cout << scprintf(" (%2d%% of bound)", (int) (100 * (absint
                  / integralbound)));
            }
            cout << scprintf("\n");
          }
          ibuf++;
        }
      }
    }
  }

  if (permute) {
    double buff1[maxint][maxint][maxint][maxint];
    sh[0] = i;
    sh[1] = j;
    sh[2] = k;
    sh[3] = l;
    init_shell_perm(int2ev3, buffer, buff1, sh, sizes);
    check_shell_perm(int2ev3, buffer, buff1, sh, sizes, 0, 1, 2, 3);
    check_shell_perm(int2ev3, buffer, buff1, sh, sizes, 1, 0, 2, 3);
    check_shell_perm(int2ev3, buffer, buff1, sh, sizes, 0, 1, 3, 2);
    check_shell_perm(int2ev3, buffer, buff1, sh, sizes, 1, 0, 3, 2);
    check_shell_perm(int2ev3, buffer, buff1, sh, sizes, 2, 3, 0, 1);
    check_shell_perm(int2ev3, buffer, buff1, sh, sizes, 2, 3, 1, 0);
    check_shell_perm(int2ev3, buffer, buff1, sh, sizes, 3, 2, 0, 1);
    check_shell_perm(int2ev3, buffer, buff1, sh, sizes, 3, 2, 1, 0);
  }

  if (bounds) {
    int boundij = int2ev3->erep_4bound(i, j, -1, -1);
    int boundkl = int2ev3->erep_4bound(-1, -1, k, l);
    int badbound = 0;
    if (boundij < boundijkl || boundkl < boundijkl) {
      badbound = 1;
    }
    if (badbound || printbounds) {
      cout << scprintf("max(%d,%d,%d,%d)=%7.4f, bnd=%7.4f, "
        "bnd(%d,%d,*,*)=%7.4f, bnd(*,*,%d,%d)=%7.4f\n", i, j, k, l,
                       maxintegral, integralbound, i, j,
                       int2ev3->logbound_to_bound(boundij), k, l,
                       int2ev3->logbound_to_bound(boundkl));
    }
    if (badbound) {
      cout << scprintf("ERROR: bad bound\n");
      abort();
    }
  }
}

void do_4_center_test(const Ref<Int2eV3> &int2ev3, int print, int printbounds,
                      int bounds, int permute, const Ref<KeyVal>& keyval) {
  int ii, jj, kk, ll, i, j, k, l, ibuf;
  int nshell = int2ev3->basis()->nshell();
  int unique = keyval->booleanvalue("unique");
  int timestats = keyval->booleanvalue("timestats");
  Ref<RegionTimer> timer = new RegionTimer();

  if (!timestats) {
    for (i = 0; i < nshell; i++) {
      int jmax = nshell - 1;
      if (unique)
        jmax = i;
      for (j = 0; j <= jmax; j++) {
        int kmax = nshell - 1;
        if (unique)
          kmax = i;
        for (k = 0; k <= kmax; k++) {
          int lmax = nshell - 1;
          if (unique) {
            if (k == i)
              lmax = j;
            else
              lmax = k;
          }
          for (l = 0; l <= lmax; l++) {
            do_shell_quartet_test(int2ev3, print, printbounds, bounds, permute,
                                  keyval, i, j, k, l);
          }
        }
      }
    }
  }
  if (timestats && nshell) {
    unsigned short seed = 1234;
    seed48(&seed);
    const int nsample = 5000;
    const int ntrials = 50;
    double times[ntrials];
    for (i = 0; i < ntrials; i++) {
      double t1 = timer->get_cpu_time();
      for (j = 0; j < nsample; j++) {
        // pick an integral at random
        int ish = int(drand48() * nshell);
        int jsh = int(drand48() * ish);
        int ksh = int(drand48() * ish);
        int lsh;
        if (ish == ksh)
          lsh = int(drand48() * jsh);
        else
          lsh = int(drand48() * ksh);
        int sh[4], sizes[4];
        if (ish >= nshell)
          ish = nshell - 1;
        if (jsh >= nshell)
          jsh = nshell - 1;
        if (ksh >= nshell)
          ksh = nshell - 1;
        if (lsh >= nshell)
          lsh = nshell - 1;
        sh[0] = ish;
        sh[1] = jsh;
        sh[2] = ksh;
        sh[3] = lsh;
        int2ev3->erep(sh, sizes);
      }
      double t2 = timer->get_cpu_time();
      times[i] = t2 - t1;
    }
    double ave = 0.0;
    for (i = 0; i < ntrials; i++) {
      ave += times[i];
    }
    ave /= ntrials;
    double sigma2 = 0.0;
    for (i = 0; i < ntrials; i++) {
      double diff = times[i] - ave;
      sigma2 += diff * diff;
    }
    double sigma = sqrt(sigma2 / ntrials);
    // adjust sigma and ave from the trial average results to
    // the integral results
    ave /= nsample;
    sigma /= sqrt(double(nsample));
    cout << scprintf(" ave = %10.8f sigma = %10.8f (microsecs)\n"
      " sigma/ave = %10.4f", ave * 1e6, sigma * 1e6, sigma / ave) << endl;
  }
}

void test_4_center(const Ref<KeyVal>& keyval, const Ref<Int2eV3> &int2ev3) {
  int i;

  cout << scprintf("4 center test:\n");
  cout << scprintf("  on entry int2ev3 used %d bytes\n",
                   int2ev3->used_storage());

  int2ev3->set_permute(0);
  int2ev3->set_redundant(1);

  int storage = keyval->intvalue("storage") - int2ev3->used_storage();
  if (storage < 0)
    storage = 0;
  if (keyval->booleanvalue("store_integrals"))
    storage = 0;
  int niter = keyval->intvalue("niter");
  int print = keyval->booleanvalue("print");
  int bounds = keyval->booleanvalue("bounds");
  int permute = keyval->booleanvalue("permute");
  int printbounds = keyval->booleanvalue("printbounds");

  cout << scprintf("  storage   = %d\n", storage);
  cout << scprintf("  niter     = %d\n", niter);
  cout << scprintf("  print     = %d\n", print);
  cout << scprintf("  bounds    = %d\n", bounds);
  cout << scprintf("  permute   = %d\n", permute);
  cout << scprintf("printbounds = %d\n", printbounds);

  if (bounds)
    int2ev3->init_bounds();

  int2ev3->init_storage(storage);

  for (i = 0; i < niter; i++) {
    do_4_center_test(int2ev3, print, printbounds, bounds, permute, keyval);
  }

  if (keyval->count("quartet") == 4) {
    do_shell_quartet_test(int2ev3, print, printbounds, bounds, permute, keyval,
                          keyval->intvalue("quartet", 0),
                          keyval->intvalue("quartet", 1),
                          keyval->intvalue("quartet", 2),
                          keyval->intvalue("quartet", 3));
  }

  int2ev3->done_storage();
  int2ev3->done_bounds();
}

void do_shell_quartet_der_test(const Ref<Int2eV3> &int2ev3, double* buffer,
                               int print, int printbounds, int bounds,
                               int permute, const Ref<KeyVal>& keyval, int i,
                               int j, int k, int l) {
  int ii, jj, kk, ll, ibuf, ider, xyz;
  der_centersv3_t dercenters;

  int sh[4], sizes[4];
  sh[0] = i;
  sh[1] = j;
  sh[2] = k;
  sh[3] = l;
  double maxintegral = 0.0, integralbound;
  int boundijkl;
  if (bounds) {
    integralbound = int2ev3->logbound_to_bound((boundijkl
        = int2ev3->erep_4bound_1der(i, j, k, l)));
  }
  int2ev3->erep_all1der(sh, sizes, &dercenters);
  ibuf = 0;
  for (ider = 0; ider < dercenters.n; ider++) {
    for (xyz = 0; xyz < 3; xyz++) {
      for (ii = 0; ii < sizes[0]; ii++) {
        for (jj = 0; jj < sizes[1]; jj++) {
          for (kk = 0; kk < sizes[2]; kk++) {
            for (ll = 0; ll < sizes[3]; ll++) {
              double absint = fabs(buffer[ibuf]);
              if (absint > maxintegral) {
                maxintegral = absint;
              }
              if (bounds && absint > integralbound) {
                cout << scprintf("((%d %d)(%d %d)|(%d %d)(%d %d))"
                  " = %15.11f, bound = %15.11f\n", sh[0], ii, sh[1], jj, sh[2],
                                 kk, sh[3], ll, buffer[ibuf], integralbound);
                abort();
              }
              if (print && absint > 1.0e-15) {
                cout << scprintf(" ((%d %d)(%d %d)"
                  "|(%d %d)(%d %d))(%d %d)"
                  " = %15.11f\n", sh[0], ii, sh[1], jj, sh[2], kk, sh[3], ll,
                                 dercenters.num[ider], xyz, buffer[ibuf]);
              }
              ibuf++;
            }
          }
        }
      }
    }
  }

  if (bounds) {
    int boundij = int2ev3->erep_4bound_1der(i, j, -1, -1);
    int boundkl = int2ev3->erep_4bound_1der(-1, -1, k, l);
    int badbound = 0;
    if (boundij < boundijkl || boundkl < boundijkl) {
      badbound = 1;
    }
    if (badbound || printbounds) {
      cout << scprintf("max(%d,%d,%d,%d)=%7.4f, bnd=%7.4f, "
        "bnd(%d,%d,*,*)=%8.4f, bnd(*,*,%d,%d)=%8.4f\n", i, j, k, l,
                       maxintegral, integralbound, i, j,
                       int2ev3->logbound_to_bound(boundij), k, l,
                       int2ev3->logbound_to_bound(boundkl));
    }
    if (badbound) {
      cout << scprintf("ERROR: bad bound\n");
      abort();
    }
  }
}

void do_test_4der_center(const Ref<Int2eV3> &int2ev3, double* buffer,
                         int print, int printbounds, int bounds, int permute,
                         const Ref<KeyVal>& keyval) {
  int i, j, k, l;
  int nshell = int2ev3->basis()->nshell();
  for (i = 0; i < nshell; i++) {
    for (j = 0; j < nshell; j++) {
      for (k = 0; k < nshell; k++) {
        for (l = 0; l < nshell; l++) {
          do_shell_quartet_der_test(int2ev3, buffer, print, printbounds,
                                    bounds, permute, keyval, i, j, k, l);
        }
      }
    }
  }
}

void test_4der_center(const Ref<KeyVal>& keyval, const Ref<Int2eV3> &int2ev3) {
  int i;

  int2ev3->set_permute(0);
  int2ev3->set_redundant(1);
  double *buffer = int2ev3->buffer();

  int niter = keyval->intvalue("niter");
  int print = keyval->booleanvalue("print");
  int bounds = keyval->booleanvalue("bounds");
  int printbounds = keyval->booleanvalue("printbounds");
  int permute = keyval->booleanvalue("permute");

  cout << scprintf("4 center derivative test:\n");
  cout << scprintf("  niter      = %d\n", niter);
  cout << scprintf("  print      = %d\n", print);
  cout << scprintf("  bounds     = %d\n", bounds);
  cout << scprintf("printbounds  = %d\n", printbounds);
  cout << scprintf("  permute    = %d\n", permute);

  if (bounds)
    int2ev3->init_bounds_1der();

  for (i = 0; i < niter; i++) {
    do_test_4der_center(int2ev3, buffer, print, printbounds, bounds, permute,
                        keyval);
  }

  if (keyval->count("quartet") == 4) {
    do_shell_quartet_der_test(int2ev3, buffer, print, printbounds, bounds,
                              permute, keyval, keyval->intvalue("quartet", 0),
                              keyval->intvalue("quartet", 1),
                              keyval->intvalue("quartet", 2),
                              keyval->intvalue("quartet", 3));
  }

  if (bounds)
    int2ev3->done_bounds_1der();
}

void test_p4(const Ref<KeyVal>& tkeyval, const Ref<Integral>& ints) {

  Ref<Wavefunction>
      wfn =
          require_dynamic_cast<Wavefunction*> (
                                               tkeyval->describedclassvalue(
                                                                            "wfn").pointer(),
                                               "test_p4\n");

  const RefSymmSCMatrix P = wfn->density();
  const Ref<GaussianBasisSet> bs = wfn->basis();
  ints->set_basis(bs, bs);
  Ref<OneBodyInt> obint = ints->p4();
  Ref<PetiteList> pl = ints->petite_list();
  RefSCDimension so_dim = pl->SO_basisdim();
  RefSCDimension ao_dim = pl->AO_basisdim();
  Ref<SCMatrixKit> so_kit = bs->so_matrixkit();
  Ref<SCMatrixKit> kit = bs->matrixkit();

  // compute p^4 in SO basis
  RefSymmSCMatrix p4_skel(ao_dim, kit);
  p4_skel.assign(0.0);
  Ref<SCElementOp> p4intop =
      new OneBodyIntOp(new SymmOneBodyIntIter(obint, pl));
  p4_skel.element_op(p4intop);
  p4intop = 0;
  RefSymmSCMatrix p4(so_dim, so_kit);
  pl->symmetrize(p4_skel, p4);
  p4_skel = 0;

  const double c = 137.0359895;
  const double mv = -(p4 * P).trace() / (8.0 * c * c);
  std::cout << "mass-velocity term = " << mv << std::endl;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
