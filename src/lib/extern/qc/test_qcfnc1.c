//
// test_qcfnc1.c
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

#include <stdlib.h>
#include <stdio.h>
#include <extern/qc/qc_fnc.h>

int
main(int argc, char** argv) {

  int natoms = 3;
  double *Z = (double*)malloc(natoms*sizeof(double));
  double *geom = (double*)malloc(natoms*3*sizeof(double));
  int use_symmetry = 1;
  int nshell;
  int s1, s2, s3, s4;
  int bf1_offset, bf2_offset, bf3_offset, bf4_offset;
  int bf1, bf2, bf3, bf4, bf12, bf1234;
  int nbf1, nbf2, nbf3, nbf4;
  const double* buffer;
  int* shell2nfunction;
  int* shell2function;

  Z[0] = 8.0; geom[0] = 0.0; geom[1] = 0.0; geom[2] = 0.0;
  Z[1] = 1.0; geom[3] = 0.0; geom[4] = 1.0; geom[5] = 1.0;
  Z[2] = 1.0; geom[6] = 0.0; geom[7] =-1.0; geom[8] = 1.0;

  init_molecule_(natoms, Z, geom, use_symmetry);

  init_basis_set_("cc-pVDZ",7);   /* 7 is the number of characters in string "cc-pVDZ".
                                     this argument is required to allow calls from Fortran programs.
                                     in practice you want to call this function as such:

                                     const char basis_set_name[] = "cc-pVDZ";
                                     init_basis_set_(basis_set_name, strlen(basis_set_name));

                                   */
  nshell = basis_set_nshell_();
  shell2nfunction = basis_set_shell_to_nfunction_();
  shell2function = basis_set_shell_to_function_();

  init_integrals_();

  init_overlap_integrals_();
  buffer = overlap_integrals_buffer_();
  fprintf(stdout, "overlap integrals:\n");
  for(s1=0; s1<nshell; s1++) {
    bf1_offset = shell2function[s1];
    nbf1 = shell2nfunction[s1];

    for(s2=0; s2<nshell; s2++) {
      bf2_offset = shell2function[s2];
      nbf2 = shell2nfunction[s2];

      compute_overlap_shell_(s1, s2);

      bf12 = 0;
      for(bf1=0; bf1<nbf1; ++bf1) {
        for(bf2=0; bf2<nbf2; ++bf2) {
          fprintf(stdout, "%d %d %15.10lf\n", bf1+bf1_offset, bf2+bf2_offset, buffer[bf12]);
          ++bf12;
        }
      }
    }
  }
  done_overlap_integrals_();

  init_hcore_integrals_();
  buffer = hcore_integrals_buffer_();
  fprintf(stdout, "hcore integrals:\n");
  for(s1=0; s1<nshell; s1++) {
    bf1_offset = shell2function[s1];
    nbf1 = shell2nfunction[s1];

    for(s2=0; s2<nshell; s2++) {
      bf2_offset = shell2function[s2];
      nbf2 = shell2nfunction[s2];

      compute_hcore_shell_(s1, s2);

      bf12 = 0;
      for(bf1=0; bf1<nbf1; ++bf1) {
        for(bf2=0; bf2<nbf2; ++bf2) {
          fprintf(stdout, "%d %d %15.10lf\n", bf1+bf1_offset, bf2+bf2_offset, buffer[bf12]);
          ++bf12;
        }
      }
    }
  }
  done_hcore_integrals_();

  init_twoecoulomb_integrals_();
  buffer = twoecoulomb_integrals_buffer_();
  fprintf(stdout, "two-e Coulomb integrals:\n");
  for(s1=0; s1<nshell; s1++) {
    bf1_offset = shell2function[s1];
    nbf1 = shell2nfunction[s1];

    for(s2=0; s2<nshell; s2++) {
      bf2_offset = shell2function[s2];
      nbf2 = shell2nfunction[s2];

      for(s3=0; s3<nshell; s3++) {
        bf3_offset = shell2function[s3];
        nbf3 = shell2nfunction[s3];

        for(s4=0; s4<nshell; s4++) {
          bf4_offset = shell2function[s4];
          nbf4 = shell2nfunction[s4];

          compute_twoecoulomb_shell_(s1, s2, s3, s4);

          bf1234 = 0;
          for(bf1=0; bf1<nbf1; ++bf1) {
            for(bf2=0; bf2<nbf2; ++bf2) {
              for(bf3=0; bf3<nbf3; ++bf3) {
                for(bf4=0; bf4<nbf4; ++bf4) {
                  fprintf(stdout, "%d %d %d %d %15.10lf\n",
                          bf1+bf1_offset, bf2+bf2_offset,
                          bf3+bf3_offset, bf4+bf4_offset,
                          buffer[bf1234]);
                  ++bf1234;
                }
              }
            }
          }

        }
      }
    }
  }
  done_twoecoulomb_integrals_();

  free(geom);
  free(Z);
  done_integrals_();
  return 0;
}
