//
// print2e.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
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

#include <util/misc/formio.h>
#include <chemistry/qc/intv3/macros.h>
#include <chemistry/qc/intv3/int2e.h>

using namespace std;
using namespace sc;

/* Prints out an integral buffer given
 * fp = where to print
 * buffer = the integrals (>>> nonredundant <<<)
 * c1 = centers structure for center 1
 * s1 = shell number on center 1
 * ...
 * This prints out integrals using the offset arrays in the
 * centers structure.  Only nonzero integrals are printed.
 */
void
Int2eV3::int_offset_print(ostream &o,
                          double *buffer,
                          Ref<GaussianBasisSet> c1, int s1,
                          Ref<GaussianBasisSet> c2, int s2,
                          Ref<GaussianBasisSet> c3, int s3,
                          Ref<GaussianBasisSet> c4, int s4)
{
  int nfunc1,nfunc2,nfunc3,nfunc4;

  nfunc1 = c1->shell(s1).nfunction();
  nfunc2 = c1->shell(s2).nfunction();
  nfunc3 = c1->shell(s3).nfunction();
  nfunc4 = c1->shell(s4).nfunction();

  int_offset_print_n(o,buffer,nfunc1,nfunc2,nfunc3,nfunc4
    ,bs1_func_offset_ + c1->shell_to_function(s1)
    ,bs2_func_offset_ + c2->shell_to_function(s2)
    ,bs3_func_offset_ + c3->shell_to_function(s3)
    ,bs4_func_offset_ + c4->shell_to_function(s4)
    ,(c2==c1)&&(s2==s1)
    ,(c3==c1)&&(s3==s1) && (c4==c2)&&(s4==s2)
    ,(c4==c3)&&(s4==s3)
    );
  }

/* Prints out an integrals buffer given the number of functions
 * on each center and shell equivalency information.
 * fp = where to print
 * buffer = the integrals (>>> nonredundant <<<)
 * n1 = number of functions in shell 1
 * ...
 * o1 = the basis function offset for shell1
 * ...
 * e12 = shell 1 == shell 2
 * e13e24 = (shell 1 == shell 3) && (shell 2 == shell 4)
 * e34 = shell 3 == shell 4
 */
void
Int2eV3::int_offset_print_n(ostream &o, double *buffer,
                            int n1, int n2, int n3, int n4,
                            int o1, int o2, int o3, int o4,
                            int e12, int e13e24, int e34)
{
  int i,j,k,l;
  int index;

  index = 0;
  for (i=0; i<=INT_MAX1(n1); i++) {
    for (j=0; j<=INT_MAX2(e12,i,n2); j++) {
      for (k=0; k<=INT_MAX3(e13e24,i,n3); k++) {
        for (l=0; l<=INT_MAX4(e13e24,e34,i,j,k,n4); l++) {
          if (INT_NONZERO(buffer[index]))
            o << scprintf(" (%2d %2d|%2d %2d) = %11.7f",
                          o1+i,o2+j,o3+k,o4+l,buffer[index])
              << endl;
          index++;
          }
        }
      }
    }
  }

/* Prints out an integral buffer given
 * fp = where to print
 * buffer = the integrals (>>> nonredundant <<<)
 * c1 = centers structure for center 1
 * s1 = shell number on center 1
 * ...
 */
void
Int2eV3::int_print(ostream &o, double *buffer,
                   Ref<GaussianBasisSet> c1, int s1,
                   Ref<GaussianBasisSet> c2, int s2,
                   Ref<GaussianBasisSet> c3, int s3,
                   Ref<GaussianBasisSet> c4, int s4)
{
  int nfunc1,nfunc2,nfunc3,nfunc4;

  nfunc1 = c1->shell(s1).nfunction();
  nfunc2 = c1->shell(s2).nfunction();
  nfunc3 = c1->shell(s3).nfunction();
  nfunc4 = c1->shell(s4).nfunction();

  int_print_n(o,buffer,nfunc1,nfunc2,nfunc3,nfunc4
    ,(c2==c1)&&(s2==s1)
    ,(c3==c1)&&(s3==s1) && (c4==c2)&&(s4==s2)
    ,(c4==c3)&&(s4==s3)
    );
  }

/* Prints out an integrals buffer given the number of functions
 * on each center and shell equivalency information.
 * fp = where to print
 * buffer = the integrals (>>> nonredundant <<<)
 * n1 = number of functions in shell 1
 * ...
 * e12 = shell 1 == shell 2
 * e13e24  = (shell 1 == shell 3) && (shell 2 == shell 4)
 * e34 = shell 3 == shell 4
 */
void
Int2eV3::int_print_n(ostream &o, double *buffer,
                     int n1, int n2, int n3, int n4,
                     int e12, int e13e24, int e34)
{
  int i,j,k,l;
  int index;

  index = 0;
  for (i=0; i<=INT_MAX1(n1); i++) {
    for (j=0; j<=INT_MAX2(e12,i,n2); j++) {
      for (k=0; k<=INT_MAX3(e13e24,i,n3); k++) {
        for (l=0; l<=INT_MAX4(e13e24,e34,i,j,k,n4); l++) {
          if (INT_NONZERO(buffer[index]))
            o << scprintf(" (%2d %2d|%2d %2d) = (%4d) = %11.7f",
                          i,j,k,l,index,buffer[index])
              << endl;
          index++;
          }
        }
      }
    }
  }

void
Int2eV3::int_print_intermediates(ostream &o)
{
  o << "The integral intermediates:" << endl;

  o << "  int_prim_zeta:" << endl;
  int_prim_zeta.print(o);

  o << "  int_prim_k:" << endl;
  int_prim_k.print(o);

  o << "  int_prim_oo2zeta:" << endl;
  int_prim_oo2zeta.print(o);

  o << "  int_prim_p:" << endl;
  int_prim_p.print(o);
  }

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
