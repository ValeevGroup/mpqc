//
// inttest.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
// Maintainer: LPS
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
#include <chemistry/qc/intv3/int1e.h>
#include <chemistry/qc/intv3/int2e.h>
#include <chemistry/qc/intv3/intv3.h>

void test_int_shell_1e(const RefKeyVal&, const RefInt1eV3 &int1ev3,
                       void (Int1eV3::*int_shell_1e)(int,int),
                       int permute);
void test_3_center(const RefKeyVal&, const RefInt2eV3 &);
void test_4_center(const RefKeyVal& keyval, const RefInt2eV3 &int2ev3);
void test_4der_center(const RefKeyVal&, const RefInt2eV3 &int2ev3);

#define maxint 9

void
testint(const RefOneBodyInt& in)
{
  if (in.null()) {
      cout << "null integral generator" << endl;
      abort();
    }
}

void
testint(const RefOneBodyDerivInt& in)
{
  if (in.null()) {
      cout << "null integral generator" << endl;
      abort();
    }
}

void
testint(const RefTwoBodyInt& in)
{
  if (in.null()) {
      cout << "null integral generator" << endl;
      abort();
    }
}

void
testint(const RefTwoBodyDerivInt& in)
{
  if (in.null()) {
      cout << "null integral generator" << endl;
      abort();
    }
}

main()
{
  int ii, i,j,k,l,m,n;

  RefRegionTimer tim = new RegionTimer("inttest", 1, 1);

  char *infile = new char[strlen(SRCDIR)+strlen("/inttest.in")+1];
  sprintf(infile,SRCDIR "/inttest.in");

  RefKeyVal pkv(new ParsedKeyVal(infile));
  RefKeyVal tkeyval(new PrefixKeyVal(":test", pkv));

  RefGaussianBasisSet basis = tkeyval->describedclassvalue("basis");
  RefMolecule mol = basis->molecule();

  int storage = tkeyval->intvalue("storage");
  RefInt1eV3 int1ev3 = new Int1eV3(basis,basis,1);
  RefInt2eV3 int2ev3 = new Int2eV3(basis,basis,basis,basis,
                                   1, storage);

  int permute = tkeyval->booleanvalue("permute");
  tim->enter("overlap");
  if (tkeyval->booleanvalue("overlap")) {
      cout << scprintf("testing overlap:\n");
      test_int_shell_1e(tkeyval, int1ev3, Int1eV3::overlap, permute);
    }
  tim->change("kinetic");
  if (tkeyval->booleanvalue("kinetic")) {
      cout << scprintf("testing kinetic:\n");
      test_int_shell_1e(tkeyval, int1ev3, Int1eV3::kinetic, permute);
    }
  tim->change("hcore");
  if (tkeyval->booleanvalue("hcore")) {
      cout << scprintf("testing hcore:\n");
      test_int_shell_1e(tkeyval, int1ev3, Int1eV3::hcore, permute);
    }
  tim->change("nuclear");
  if (tkeyval->booleanvalue("nuclear")) {
      cout << scprintf("testing nuclear:\n");
      test_int_shell_1e(tkeyval, int1ev3, Int1eV3::nuclear, permute);
    }
  tim->change("3 center");
  if (tkeyval->booleanvalue("3")) {
      test_3_center(tkeyval, int2ev3);
    }
  tim->change("4 center");
  if (tkeyval->booleanvalue("4")) {
      test_4_center(tkeyval, int2ev3);
    }
  tim->change("4 center der");
  if (tkeyval->booleanvalue("4der")) {
      test_4der_center(tkeyval, int2ev3);
    }

  tim->change("IntegralV3");
  RefIntegral integral = new IntegralV3(basis);
  RefOneBodyInt overlap = integral->overlap();
  testint(overlap);
  RefOneBodyInt kinetic = integral->kinetic();
  testint(kinetic);
  RefOneBodyDerivInt kinetic_der = integral->kinetic_deriv();
  testint(kinetic_der);
  RefOneBodyDerivInt overlap_der = integral->overlap_deriv();
  testint(overlap_der);
  RefOneBodyDerivInt nuclear_der = integral->nuclear_deriv();
  testint(nuclear_der);
  RefTwoBodyInt erep = integral->electron_repulsion();
  testint(erep);
  RefTwoBodyDerivInt erep_der = integral->electron_repulsion_deriv();
  testint(erep_der);
  tim->exit();

  tim->print();
  return 0;
}

void
do_shell_test_1e(const RefInt1eV3 &int1ev3,
                 void (Int1eV3::*int_shell_1e)(int,int),
                 int permute, int i, int j, int na, int nb,
                 double *buf, double *pbuf)
{
  int ii = 0;
  int a;
  double *buffer = int1ev3->buffer();
  (int1ev3->*int_shell_1e)(i, j);
  for (a=0; a<na*nb; a++) {
      buf[a] = buffer[a];
    }
  (int1ev3->*int_shell_1e)(j, i);
  for (a=0; a<na*nb; a++) {
      pbuf[a] = buffer[a];
    }
  for (a=0; a<na; a++) {
      for (int b=0; b<nb; b++) {
          if (fabs(buf[ii] - pbuf[a + na*b]) > 1.0e-13) {
              cout << scprintf("----- 1e perm failed:"
                               "<%d %d|%d %d>:"
                               " %18.14f != %18.14f "
                               "<%d %d|%d %d>\n",
                               i, a, j, b,
                               buf[ii],
                               pbuf[a + na*b],
                               j, b, i, a);
            }
          if (fabs(buf[ii]) > 1.0e-15) {
              cout << scprintf(" <(%d %d)|(%d %d)> = %15.11f\n",
                               i,a,j,b, buf[ii]);
            }
          ii++;
        }
    }
}

void
test_int_shell_1e(const RefKeyVal& keyval, const RefInt1eV3 &int1ev3,
                  void (Int1eV3::*int_shell_1e)(int,int),
                  int permute)
{
  int flags = 0;
  RefGaussianBasisSet basis = int1ev3->basis();
  int maxfunc = basis->max_nfunction_in_shell();
  int size = maxfunc * maxfunc;
  double *buf = new double[size];
  double *pbuf = new double[size];
  int nshell = int1ev3->basis()->nshell();

  for (int i=0; i<nshell; i++) {
      int na = basis->shell(i).nfunction();
      for (int j=0; j<nshell; j++) {
          int nb = basis->shell(j).nfunction();
          do_shell_test_1e(int1ev3, int_shell_1e, permute,
                           i, j, na, nb, buf, pbuf);

        }
    }

  delete[] buf;
  delete[] pbuf;
}

void
test_3_center(const RefKeyVal& keyval, const RefInt2eV3 &int2ev3)
{
  int ii, i,j,k,l,m,n;

  int2ev3->set_redundant(1);
  int2ev3->set_permute(0);
  double *buffer = int2ev3->buffer();
  int nshell = int2ev3->basis()->nshell();

  for (i=0; i<nshell; i++) {
      for (j=0; j<nshell; j++) {
          int sh[2], sizes[2];
          sh[0] = i;
          sh[1] = j;
          int2ev3->erep_2center(sh,sizes);
          ii = 0;
          for (k=0; k<sizes[0]; k++) {
              for (l=0; l<sizes[1]; l++) {
                  if (fabs(buffer[ii])>1.0e-15)
                      cout << scprintf(" ((%d %d)|(%d %d)) = %15.11f\n",
                                       sh[0],k,sh[1],l, buffer[ii]);
                  ii++;
                }
            }
        }
    }

  for (i=0; i<nshell; i++) {
      for (j=0; j<nshell; j++) {
          for (m=0; m<nshell; m++) {
              int sh[3], sizes[3];
              sh[0] = i;
              sh[1] = j;
              sh[2] = m;
              int2ev3->erep_3center(sh,sizes);
              ii = 0;
              for (k=0; k<sizes[0]; k++) {
                  for (l=0; l<sizes[1]; l++) {
                      for (n=0; n<sizes[2]; n++) {
                          if (fabs(buffer[ii])>1.0e-15)
                              cout << scprintf(
                                  " ((%d %d)|(%d %d)(%d %d)) = %15.11f\n",
                                  sh[0],k,sh[1],l,sh[2],n, buffer[ii]);
                          ii++;
                        }
                    }
                }
            }
        }
    }

}

void
init_shell_perm(const RefInt2eV3 &int2ev3, double *integrals,
                double buff[maxint][maxint][maxint][maxint],
                int sh[4], int sizes[4])
{
  int i, j, k, l;
  int oldp = int2ev3->permute();
  int2ev3->set_permute(0);
  int2ev3->erep(sh, sizes);
  int2ev3->set_permute(oldp);
  for (i=0; i<sizes[0]; i++) {
      for (j=0; j<sizes[1]; j++) {
          for (k=0; k<sizes[2]; k++) {
              for (l=0; l<sizes[3]; l++) {
                  buff[i][j][k][l] = *integrals++;
                }
            }
        }
    }
}

void
check_shell_perm(const RefInt2eV3 &int2ev3, double *integrals,
                 double buff[maxint][maxint][maxint][maxint],
                 int sh[4], int sizes[4], int p0, int p1, int p2, int p3)
{
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
  for (i[0]=0; i[0]<psizes[0]; i[0]++) {
      for (i[1]=0; i[1]<psizes[1]; i[1]++) {
          for (i[2]=0; i[2]<psizes[2]; i[2]++) {
              for (i[3]=0; i[3]<psizes[3]; i[3]++) {
                  if (fabs(buff[i[ip[0]]][i[ip[1]]][i[ip[2]]][i[ip[3]]]
                           - integrals[index]) > 1.0e-13) {
                      cout << scprintf("perm %d %d %d %d failed:"
                             "((%d %d)(%d %d)|(%d %d)(%d %d)):"
                             " %18.14f != %18.14f "
                             "((%d %d)(%d %d)|(%d %d)(%d %d))\n",
                             p0, p1, p2, p3,
                             sh[0],i[0], sh[1],i[1], sh[2],i[2], sh[3],i[3],
                             buff[i[ip[0]]][i[ip[1]]][i[ip[2]]][i[ip[3]]],
                             integrals[index],
                             psh[0],i[p[0]], psh[1],i[p[1]],
                             psh[2],i[p[2]], psh[3],i[p[3]]);
                    }
                  index++;
                }
            }
        }
    }
}

void
do_shell_quartet_test(const RefInt2eV3 &int2ev3,
                      int print, int printbounds, int bounds, int permute,
                      const RefKeyVal& keyval,
                      int i, int j, int k, int l)
{
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
      integralbound
          = int2ev3->logbound_to_bound(
              (boundijkl = int2ev3->erep_4bound(i,j,k,l))
              );
    }
  double *buffer = int2ev3->buffer();
  int2ev3->erep(sh,sizes);
  ibuf = 0;
  maxintegral = 0.0;
  for (ii=0; ii<sizes[0]; ii++) {
      for (jj=0; jj<sizes[1]; jj++) {
          for (kk=0; kk<sizes[2]; kk++) {
              for (ll=0; ll<sizes[3]; ll++) {
                  double absint = fabs(buffer[ibuf]);
                  if (absint > maxintegral) {
                      maxintegral = absint;
                    }
                  if (bounds &&  absint > integralbound) {
                      cout << scprintf("((%d %d)(%d %d)|(%d %d)(%d %d)) = %15.11f, "
                             "bound = %15.11f\n",
                             sh[0], ii, sh[1], jj, sh[2], kk, sh[3], ll,
                             buffer[ibuf], integralbound);
                      abort();
                    }
                  if (print && (absint > 1.0e-9
                                || (bounds && integralbound > 1.0e-9))) {
                      cout << scprintf(" ((%d %d)(%d %d)|(%d %d)(%d %d))"
                             " = %15.11f",
                             sh[0],ii,
                             sh[1],jj,
                             sh[2],kk,
                             sh[3],ll,
                             buffer[ibuf]);
                      if (bounds) {
                          cout << scprintf(" (%2d%% of bound)",
                                 (int)(100*(absint/integralbound)));
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
      int boundij = int2ev3->erep_4bound(i,j,-1,-1);
      int boundkl = int2ev3->erep_4bound(-1,-1,k,l);
      int badbound = 0;
      if (boundij < boundijkl || boundkl < boundijkl) {
          badbound = 1;
        }
      if (badbound || printbounds) {
          cout << scprintf("max(%d,%d,%d,%d)=%7.4f, bnd=%7.4f, "
                 "bnd(%d,%d,*,*)=%7.4f, bnd(*,*,%d,%d)=%7.4f\n",
                 i, j, k, l, maxintegral, integralbound,
                 i,j, int2ev3->logbound_to_bound(boundij),
                 k,l, int2ev3->logbound_to_bound(boundkl));
        }
      if (badbound) {
          cout << scprintf("ERROR: bad bound\n");
          abort();
        }
    }
}

void
do_4_center_test(const RefInt2eV3 &int2ev3, int print, int printbounds,
                 int bounds, int permute,
                 const RefKeyVal& keyval)
{
  int ii,jj,kk,ll, i,j,k,l, ibuf;
  int nshell = int2ev3->basis()->nshell();

  for (i=0; i<nshell; i++) {
      for (j=0; j<nshell; j++) {
          for (k=0; k<nshell; k++) {
              for (l=0; l<nshell; l++) {
                  do_shell_quartet_test(int2ev3, print, printbounds,
                                        bounds, permute,
                                        keyval, i, j, k, l);
                }
            }
        }
    }
}

void
test_4_center(const RefKeyVal& keyval, const RefInt2eV3 &int2ev3)
{
  int i;

  cout << scprintf("4 center test:\n");
  cout << scprintf("  on entry int2ev3 used %d bytes\n", int2ev3->used_storage());

  int2ev3->set_permute(0);
  int2ev3->set_redundant(1);

  int storage = keyval->intvalue("storage") - int2ev3->used_storage();
  if (storage < 0) storage = 0;
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

  if (bounds) int2ev3->init_bounds();

  int2ev3->init_storage(storage);

  for (i=0; i<niter; i++) {
      do_4_center_test(int2ev3, print, printbounds, bounds, permute, keyval);
    }

  if (keyval->count("quartet") == 4) {
      do_shell_quartet_test(int2ev3, print, printbounds, bounds, permute,
                            keyval,
                            keyval->intvalue("quartet", 0),
                            keyval->intvalue("quartet", 1),
                            keyval->intvalue("quartet", 2),
                            keyval->intvalue("quartet", 3));
    }

  int2ev3->done_storage();
  int2ev3->done_bounds();
}

void
do_shell_quartet_der_test(const RefInt2eV3 &int2ev3,
                          double* buffer, int print, int printbounds,
                          int bounds, int permute,
                          const RefKeyVal& keyval,
                          int i, int j, int k, int l)
{
  int ii,jj,kk,ll, ibuf, ider, xyz;
  der_centersv3_t dercenters;

  int sh[4], sizes[4];
  sh[0] = i;
  sh[1] = j;
  sh[2] = k;
  sh[3] = l;
  double maxintegral = 0.0, integralbound;
  int boundijkl;
  if (bounds) {
      integralbound
          = int2ev3->logbound_to_bound(
              (boundijkl = int2ev3->erep_4bound_1der(i,j,k,l))
              );
    }
  int2ev3->erep_all1der(sh,sizes,&dercenters);
  ibuf = 0;
  for (ider=0; ider<dercenters.n; ider++) {
      for (xyz=0; xyz<3; xyz++) {
          for (ii=0; ii<sizes[0]; ii++) {
              for (jj=0; jj<sizes[1]; jj++) {
                  for (kk=0; kk<sizes[2]; kk++) {
                      for (ll=0; ll<sizes[3]; ll++) {
                          double absint = fabs(buffer[ibuf]);
                          if (absint > maxintegral) {
                              maxintegral = absint;
                            }
                          if (bounds &&  absint > integralbound) {
                              cout << scprintf("((%d %d)(%d %d)|(%d %d)(%d %d))"
                                     " = %15.11f, bound = %15.11f\n",
                                     sh[0], ii, sh[1], jj,
                                     sh[2], kk, sh[3], ll,
                                     buffer[ibuf], integralbound);
                              abort();
                            }
                          if (print && absint > 1.0e-15) {
                              cout << scprintf(" ((%d %d)(%d %d)"
                                     "|(%d %d)(%d %d))(%d %d)"
                                     " = %15.11f\n",
                                     sh[0],ii,
                                     sh[1],jj,
                                     sh[2],kk,
                                     sh[3],ll,
                                     dercenters.num[ider], xyz,
                                     buffer[ibuf]
                                  );
                            }
                          ibuf++;
                        }
                    }
                }
            }
        }
    }

  if (bounds) {
      int boundij = int2ev3->erep_4bound_1der(i,j,-1,-1);
      int boundkl = int2ev3->erep_4bound_1der(-1,-1,k,l);
      int badbound = 0;
      if (boundij < boundijkl || boundkl < boundijkl) {
          badbound = 1;
        }
      if (badbound || printbounds) {
          cout << scprintf("max(%d,%d,%d,%d)=%7.4f, bnd=%7.4f, "
                 "bnd(%d,%d,*,*)=%8.4f, bnd(*,*,%d,%d)=%8.4f\n",
                 i, j, k, l, maxintegral, integralbound,
                 i,j, int2ev3->logbound_to_bound(boundij),
                 k,l, int2ev3->logbound_to_bound(boundkl));
        }
      if (badbound) {
          cout << scprintf("ERROR: bad bound\n");
          abort();
        }
    }
}

void
do_test_4der_center(const RefInt2eV3 &int2ev3,
                    double* buffer, int print, int printbounds,
                    int bounds, int permute,
                    const RefKeyVal& keyval)
{
  int i,j,k,l;
  int nshell = int2ev3->basis()->nshell();
  for (i=0; i<nshell; i++) {
      for (j=0; j<nshell; j++) {
          for (k=0; k<nshell; k++) {
              for (l=0; l<nshell; l++) {
                  do_shell_quartet_der_test(int2ev3, buffer,
                                            print, printbounds,
                                            bounds, permute,
                                            keyval,
                                            i, j, k, l);
                }
            }
        }
    }
}

void
test_4der_center(const RefKeyVal& keyval, const RefInt2eV3 &int2ev3)
{
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

  if (bounds) int2ev3->init_bounds_1der();

  for (i=0; i<niter; i++) {
      do_test_4der_center(int2ev3, buffer,
                          print, printbounds, bounds, permute, keyval);
    }

  if (keyval->count("quartet") == 4) {
      do_shell_quartet_der_test(int2ev3, buffer, print, printbounds,
                                bounds, permute,
                                keyval,
                                keyval->intvalue("quartet", 0),
                                keyval->intvalue("quartet", 1),
                                keyval->intvalue("quartet", 2),
                                keyval->intvalue("quartet", 3));
    }

  if (bounds) int2ev3->done_bounds_1der();
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
// End:
