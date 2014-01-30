
/*
 * Copyright 2009 Sandia Corporation. Under the terms of Contract
 * DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government
 * retains certain rights in this software.
 *
 * This file is a part of the MPQC LMP2 library.
 *
 * The MPQC LMP2 library is free software: you can redistribute it
 * and/or modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation, either
 * version 3 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this program.  If not, see
 * <http://www.gnu.org/licenses/>.
 *
 */

#include <iostream>
#include <map>
#include <vector>
#include <stdexcept>

#include <math/scmat/matrix.h>
#include <math/scmat/local.h>

#include <chemistry/qc/lmp2/sma.h>

#undef r2 // r2 is defined in gcc 3.4

namespace sc {

template <int N>
void
initialize_array_values(sma2::Array<N> &array, double offset = 0.0)
{
  // initialize the data in array with some values
  array.allocate_blocks();
  for (typename sma2::Array<N>::blockmap_t::const_iterator
           iter = array.blockmap().begin();
       iter != array.blockmap().end();
       iter++) {
      const sma2::BlockInfo<N> &bi = iter->first;
      double *data = iter->second;
      int sz = array.block_size(bi);
      int blockval = 0;
      int blockmul = 1;
      for (int i=N-1; i>=0; i--) {
          blockval += bi.block(i) * blockmul;
          blockmul *= 10;
        }
      for (int i=0; i<sz; i++) {
          *data++ = 10000*i + blockval + offset;
        }
    }
}

void
initialize_sparse_array(sma2::Array<4>&fullarray,
                        std::vector<std::pair<int,int> > &pairs)
{
  const sma2::Range &r1 = fullarray.index(0);
  const sma2::Range &r2 = fullarray.index(2);
  if (r1 != fullarray.index(1)
      || r2 != fullarray.index(3)) {
      throw std::runtime_error("bad ranges to initialize_sparse_array");
    }

  // compute the domains of those r2 interacting with each r1
  std::vector<std::vector<int> > domains(r1.nblock());
  for (int i=0; i<domains.size(); i++) {
      for (int j=(i?i-1:0); j<=(i==r2.nblock()-1?i:i+1); j++) {
          domains[i].push_back(j);
        }
    }

  // compute the pair domains and fill fullarray with blocks
  for (int i=0; i<r1.nblock(); i++) {
      sma2::BlockInfo<4> bi;
      bi.block(0) = i;
      for (int j=(i?i-1:0); j<=(i==r1.nblock()-1?i:i+1); j++) {
          bi.block(1) = j;
          pairs.push_back(std::make_pair(i,j));
          for (int k=0; k<domains[i].size(); k++) {
              bi.block(2) = k;
              for (int l=0; l<domains[j].size(); l++) {
                  bi.block(3) = l;
                  fullarray.add_unallocated_block(bi);
                }
            }
        }
    }

  initialize_array_values(fullarray);
}

template <int N>
void
initialize_dense_array(sma2::Array<N> &array, double offset = 0.0)
{
  array.add_all_unallocated_blocks();
  initialize_array_values(array, offset);
}

void
test_fixed_values()
{
  sma2::Range r1(4, 1);
  sma2::Range r2(7, 2);

  sma2::Array<4> fullarray(r1, r1, r2, r2);
  sma2::Array<4> fullarray2(r1, r1, r2, r2);

  std::vector<std::pair<int,int> > pairs;
  initialize_sparse_array(fullarray,pairs);

  std::cout << "fullarray:" << std::endl;
  fullarray.print();

  std::map<std::pair<int,int>,sma2::Array<2>*> subarraymap;

  for (int pairiter=0; pairiter<pairs.size(); pairiter++) {
      int i = pairs[pairiter].first;
      int j = pairs[pairiter].second;
      sma2::Index I(i), J(j), P("p"), Q("q");
      subarraymap[pairs[pairiter]] = new sma2::Array<2>(r2,r2);
      sma2::Array<2> &subarray = *subarraymap[pairs[pairiter]];
//        sma2::BlockInfo<2> subbi;
//        for (int k=0; k<domains[i].size(); k++) {
//            subbi.block(0) = k;
//            for (int l=0; l<domains[j].size(); l++) {
//                subbi.block(1) = l;
//                subarray.add_unallocated_block(subbi);
//              }
//          }
      subarray(P,Q) |= fullarray(I,J,P,Q);
      subarray(P,Q) = fullarray(I,J,P,Q);
      std::cout << "subarray(" << i << "," << j << "):" << std::endl;
      subarray.print();
    }

  fullarray2.init_blocks(fullarray);
  fullarray2.allocate_blocks();
  fullarray2.zero();

  for (int pairiter=0; pairiter<pairs.size(); pairiter++) {
      int i = pairs[pairiter].first;
      int j = pairs[pairiter].second;
      sma2::Index I(i), J(j), P("p"), Q("q");
      sma2::Array<2> &subarray = *subarraymap[pairs[pairiter]];
      fullarray2(I,J,P,Q) += subarray(P,Q);
    }

  std::cout << "fullarray2:" << std::endl;
  fullarray2.print();
}

void
test_remapping()
{
  sma2::Range r1(7, 1);
  sma2::Range r2(4, 2);

  sma2::Array<4> T(r1, r2, r1, r2);
  sma2::Array<4> D(r2, r2, r1, r1);

  std::vector<std::pair<int,int> > pairs;
  initialize_sparse_array(D,pairs);

  std::cout << "D:" << std::endl;
  D.print();

  T("p","i","q","j") |= D("i","j","p","q");
  T("p","i","q","j") = D("i","j","p","q");

  std::cout << "T:" << std::endl;
  T.print();

}

void
test_assign()
{
  std::cout << "entered test_assign" << std::endl;

  sma2::Range o(1,1);
  sma2::Range v(4,2);

  sma2::Array<4> C(o,o,v,v), A(v,v,o,o);

  initialize_dense_array(A);
  std::cout << "initial A:" << std::endl;
  A.print();

  sma2::Index I("i"), J("j"), R("r"), S("s");

  C(I,J,R,S) |= A(R,S,I,J);

  std::cout << "C after |=:" << std::endl;
  C.print();

  C(I,J,R,S) = A(R,S,I,J);

  std::cout << "C after =:" << std::endl;
  C.print();
}

void
test_sum()
{
  std::cout << "entered test_sum" << std::endl;

  sma2::Range o(1,1);
  sma2::Range v(3,2);

  sma2::Array<2> C(o,v), A(v,o);

  initialize_dense_array(A);
  std::cout << "initial A:" << std::endl;
  A.print();

  initialize_dense_array(C);
  std::cout << "initial C:" << std::endl;
  C.print();

  sma2::Index I("i"), J("j");

  C(I,J) += A(J,I);

  std::cout << "C after +=:" << std::endl;
  C.print();
}

void
test_div()
{
  std::cout << "entered test_div" << std::endl;

  sma2::Range o(1,1);
  sma2::Range v(3,2);

  sma2::Array<2> C(o,v), A(v,o);

  initialize_dense_array(A,1.0);
  std::cout << "initial A:" << std::endl;
  A.print();

  initialize_dense_array(C);
  std::cout << "initial C:" << std::endl;
  C.print();

  sma2::Index I("i"), J("j");

  C(I,J) /= A(J,I);

  std::cout << "C after /=:" << std::endl;
  C.print();
}

void
test_contract_2x2x2()
{
  std::cout << "entered test_contract_2x2x2" << std::endl;

  sma2::Range o(1,1);
  sma2::Range r(2,1);

  sma2::Array<2> A(r,r),B(r,o),C(r,o);

  sma2::Index I("i"),J("j"),K("k");

  initialize_dense_array(A);
  initialize_dense_array(B);

  C(I,J)|= A(I,K) * B(K,J);
  C(I,J) = A(I,K) * B(K,J);

  std::cout << "C:" << std::endl;
  C.print();

}

void
test_contract_2x2x2_reordered()
{
  std::cout << "entered test_contract_2x2x2_reordered" << std::endl;

  sma2::Range r(3,3);

  sma2::Array<2> A(r,r),B(r,r),C(r,r),Cref(r,r);
  sma2::Array<2> At(r,r),Bt(r,r);

  sma2::Index I("i"),J("j"),K("k");

  initialize_dense_array(A);
  initialize_dense_array(B);

  At(I,J) |= A(J,I);
  At(I,J) = A(J,I);

  Bt(I,J) |= B(J,I);
  Bt(I,J) = B(J,I);

  Cref(I,J)|= A(I,K) * B(J,K);
  Cref(I,J) = A(I,K) * B(J,K);
  std::cout << "Cref(I,J) = A(I,K) * B(J,K):" << std::endl;
  Cref.print();

  C.clear();
  C(I,J)|= At(K,I) * B(J,K);
  C(I,J) = At(K,I) * B(J,K);
  C(I,J) -= Cref(I,J);
  std::cout << "At(K,I) * B(J,K) - Cref(I,J):" << std::endl;
  C.print();

  C.clear();
  C(I,J)|= A(I,K) * Bt(K,J);
  C(I,J) = A(I,K) * Bt(K,J);
  C(I,J) -= Cref(J,I);
  std::cout << "A(I,K) * Bt(K,J) - Cref(I,J):" << std::endl;
  C.print();

  C.clear();
  C(I,J)|= At(K,I) * Bt(K,J);
  C(I,J) = At(K,I) * Bt(K,J);
  std::cout << "At(K,I) * Bt(K,J) - Cref(I,J):" << std::endl;
  C(I,J) -= Cref(I,J);
  C.print();

  // tests for cases where contract will transpose C:

  C.clear();
  C(J,I)|= A(I,K) * B(J,K);
  C(J,I) = A(I,K) * B(J,K);
  C(J,I) -= Cref(I,J);
  std::cout << "C(J,I) = A(I,K) * B(J,K) - Cref(I,J):" << std::endl;
  C.print();

  C.clear();
  C(J,I)|= At(K,I) * B(J,K);
  C(J,I) = At(K,I) * B(J,K);
  C(J,I) -= Cref(I,J);
  std::cout << "C(J,I) = At(K,I) * B(J,K) - Cref(I,J):" << std::endl;
  C.print();

  C.clear();
  C(J,I)|= A(I,K) * Bt(K,J);
  C(J,I) = A(I,K) * Bt(K,J);
  C(J,I) -= Cref(I,J);
  std::cout << "C(J,I) = A(I,K) * Bt(K,J) - Cref(I,J):" << std::endl;
  C.print();

  C.clear();
  C(J,I)|= At(K,I) * Bt(K,J);
  C(J,I) = At(K,I) * Bt(K,J);
  std::cout << "C(J,I) = At(K,I) * Bt(K,J) - Cref(I,J):" << std::endl;
  C(J,I) -= Cref(I,J);
  C.print();
}

void
test_contract_empty()
{
  std::cout << "entered test_contract_empty" << std::endl;

  sma2::Range o(1,1);
  sma2::Range r(2,1);

  sma2::Array<2> A(r,r),B(r,o),C(r,o),D(r,o);
  sma2::Array<4> E(o,o,r,o);

  sma2::Index I("i"),J("j"),K("k"), I0(0);

  // In this case C is not empty, and A and B are empty.
  initialize_dense_array(C);
  std::cout << "C initial:" << std::endl;
  C.print();
  C(I,J) += A(I,K) * B(K,J);
  std::cout << "C final (should be the same as C initial):" << std::endl;
  C.print();

  // In this case D, A, and B are empty.
  D(I,J)|= A(I,K) * B(K,J);
  D(I,J) = A(I,K) * B(K,J);

}

void
test_contract_4x6x6()
{
  std::cout << "entered test_contract_4x6x6" << std::endl;

  sma2::Range o(1,1);
  sma2::Range r(2,1);

  sma2::Array<6> A(o,o,o,r,r,r),B(o,o,o,r,r,r);
  sma2::Array<4> C(o,o,r,r);

  sma2::Index I("i"),J("j"),K("k"),L("l");
  sma2::Index P("p"),Q("q"),R("r"),S("s");

  initialize_dense_array(A);
  initialize_dense_array(B);

  C(I,J,P,Q)|= A(I,K,L,P,R,S) * B(J,K,L,Q,R,S);
  C(I,J,P,Q) = A(I,K,L,P,R,S) * B(J,K,L,Q,R,S);

  std::cout << "C:" << std::endl;
  C.print();

}

void
test_fixed_contract_1x2x2()
{
  std::cout << "entered test_fixed_contract_1x2x2" << std::endl;

  sma2::Range o(1,1);
  sma2::Range r(2,1);

  sma2::Array<2> A(r,r),B(r,o);
  sma2::Array<1> C(o);

  sma2::Index I(0),J("j"),K("k");

  initialize_dense_array(A);
  initialize_dense_array(B);

  C(J)|= A(I,K) * B(K,J);
  C(J) = A(I,K) * B(K,J);

  std::cout << "C:" << std::endl;
  C.print();

}

void
test_fixed_contract_3x2x2()
{
  std::cout << "entered test_fixed_contract_3x2x2" << std::endl;

  sma2::Range occ_act(2,1);
  sma2::Range vir(4,2);

  sma2::Index iJ(0), iL("L"), iT("T"), iS("S");

  sma2::Array<2> t(occ_act, vir);

  initialize_dense_array(t);

  sma2::Array<3> tt_tmp(occ_act,vir,vir);

  tt_tmp(iL,iT,iS)|= t(iL,iT) * t(iJ,iS);
  tt_tmp(iL,iT,iS) = t(iL,iT) * t(iJ,iS);

  std::cout << "t" << std::endl;
  t.print();
  std::cout << "tt_tmp" << std::endl;
  tt_tmp.print();
}

void
test_pack_vector()
{
  std::cout << "entered test_pack_vector" << std::endl;

  sma2::Range r(5,2);

  sma2::Array<1> A(r);

  sc::Ref<sc::SCMatrixKit> kit = sc::SCMatrixKit::default_matrixkit();
  sc::RefSCDimension dim = new sc::SCDimension(r.nindex());
  sc::RefSCVector v(dim, kit);

  initialize_dense_array(A);

  pack_array_into_vector(A, v);

  sma2::Array<1> B(r);
  pack_vector_into_empty_array(v, B, 1e-10);

  v.print("v");
  std::cout << "A:" << std::endl;
  A.print();
  std::cout << "B:" << std::endl;
  B.print();

}

// this is not implemented
// void
// test_contract_4x2x2()
// {
//   std::cout << "entered test_contract_4x2x2" << std::endl;

//   sma2::Range o(1,1);
//   sma2::Range r(2,1);

//   sma2::Array<2> A(r,o),B(r,o),C(r,r,o,o),D(r,r,o,o);

//   sma2::Index I("i"),J("j"),K("k"),L("l");

//   initialize_dense_array(A);
//   initialize_dense_array(B);
//   initialize_dense_array(D);

//   C(I,J,K,L)|= A(I,J) * B(K,L);
//   C(I,J,K,L) = D(I,J,K,L) + A(I,K) * B(J,L);

//   std::cout << "C:" << std::endl;
//   C.print();

// }

void
test_fixed_contract_Bfixed()
{
  std::cout << "entered test_fixed_contract_Bfixed" << std::endl;

  sma2::Range r(2, 1);

  sma2::Array<2> C(r, r);
  sma2::Array<2> A(r, r);
  sma2::Array<4> B(r, r, r, r);

  initialize_dense_array(A);
  initialize_dense_array(B);

  std::cout << "A:" << std::endl;
  A.print();

  std::cout << "B:" << std::endl;
  B.print();

  for (int x = 0; x<r.nindex(); x++) {
      for (int y = 0; y<r.nindex(); y++) {
          sma2::Index X(x), Y(y), I("i"), J("j"), K("k"), L("l");
          C.clear();
          C(I,K)|= A(K,J) * B(X, Y, I, J);
          C(I,K) = A(K,J) * B(X, Y, I, J);
          std::cout << "x = " << x << " y = " << y << ":" << std::endl;
          std::cout << "C:" << std::endl;
          C.print();
        }
    }
}

void
test_fixed_contract_Afixed()
{
  std::cout << "entered test_fixed_contract_Afixed" << std::endl;

  sma2::Range r(2, 1);

  sma2::Array<2> C(r, r);
  sma2::Array<4> A(r, r, r, r);
  sma2::Array<2> B(r, r);

  initialize_dense_array(A);
  initialize_dense_array(B);

  std::cout << "A:" << std::endl;
  A.print();

  std::cout << "B:" << std::endl;
  B.print();

  for (int x = 0; x<r.nindex(); x++) {
      for (int y = 0; y<r.nindex(); y++) {
          sma2::Index X(x), Y(y), I("i"), J("j"), K("k"), L("l");
          C.clear();
          C(I,K)|= A(X, Y, K, J) * B(I, J);
          C(I,K) = A(X, Y, K, J) * B(I, J);
          std::cout << "x = " << x << " y = " << y << ":" << std::endl;
          std::cout << "C:" << std::endl;
          C.print();
        }
    }
}

void
test_fixed_contract_Cfixed()
{
  std::cout << "entered test_fixed_contract_Cfixed" << std::endl;

  sma2::Range r(2, 1);

  sma2::Array<4> C(r, r, r, r);
  sma2::Array<2> A(r, r);
  sma2::Array<2> B(r, r);

  initialize_dense_array(A);
  initialize_dense_array(B);
  initialize_dense_array(C);

  std::cout << "A:" << std::endl;
  A.print();

  std::cout << "B:" << std::endl;
  B.print();

  std::cout << "C before:" << std::endl;
  C.print();

  for (int x = 0; x<r.nindex(); x++) {
      for (int y = 0; y<r.nindex(); y++) {
          sma2::Index X(x), Y(y), I("i"), J("j"), K("k"), L("l");
          C(X,Y,I,K) += A(K, J) * B(I, J);
        }
    }

  std::cout << "C after:" << std::endl;
  C.print();
}

void
test_fixed_assign()
{
  std::cout << "entered test_fixed_assign" << std::endl;

  sma2::Range r(2, 1);

  sma2::Array<2> A(r, r);
  sma2::Array<2> B(r, r);
  sma2::Array<3> B3(r, r, r);

  initialize_dense_array(A);
  initialize_dense_array(B);
  initialize_dense_array(B3);

  sma2::Index X(0), Y(1), I("i"), J("j");

  std::cout << "B:" << std::endl;
  B.print();

  std::cout << "B3:" << std::endl;
  B3.print();

  try {
      std::cout << "A(Y,I) = B(Y,I)" << std::endl;
      A.clear();
      A(Y,I) |= B(Y,I);
      A(Y,I) = B(Y,I);
      std::cout << "A:" << std::endl;
      A.print();
    }
  catch (std::exception &e) {
      std::cout << "exception: " << e.what() << std::endl;
    }

  try {
      std::cout << "A(X,I) = B(Y,I)" << std::endl;
      A.clear();
      A(X,I) |= B(Y,I);
      A(X,I) = B(Y,I);
      std::cout << "A:" << std::endl;
      A.print();
    }
  catch (std::exception &e) {
      std::cout << "exception: " << e.what() << std::endl;
    }

  try {
      std::cout << "A(I,J) = B3(Y,I,J)" << std::endl;
      A.clear();
      A(I,J) |= B3(Y,I,J);
      A(I,J) = B3(Y,I,J);
      std::cout << "A:" << std::endl;
      A.print();
    }
  catch (std::exception &e) {
      std::cout << "exception: " << e.what() << std::endl;
    }
}

void
test_fixed_assign_2()
{
  std::cout << "entered test_fixed_assign_2" << std::endl;

  sma2::Index iR("r"), iS("s");
  sma2::Index iK("k"), iL("l");

  sma2::Range r(3, 2);
  sma2::Range o(2, 1);

  sma2::Array<4> C_ref(o,o,r,r);
  sma2::Array<4> C(o,o,r,r);
  sma2::Array<4> A(o,o,r,r);

  initialize_dense_array(A);

  C_ref(iK,iL,iR,iS) |= A(iK,iL,iR,iS);
  C_ref(iK,iL,iR,iS)  =  2.0 * A(iK,iL,iR,iS);
  C_ref(iK,iL,iR,iS) += -1.0 * A(iK,iL,iS,iR);

  C(iK,iL,iR,iS) |= A(iK,iL,iR,iS);
  C.zero();
  for (int i=0; i<o.nindex(); i++) {
      sma2::Index iI(i);
      C(iI,iL,iR,iS) +=  2.0 * A(iI,iL,iR,iS);
      C(iI,iL,iR,iS) += -1.0 * A(iI,iL,iS,iR);
    }

  std::cout << "C:" << std::endl;
  C.print();

  C(iK,iL,iR,iS) -= C_ref(iK,iL,iR,iS);

  std::cout << "C - Cref:" << std::endl;
  C.print();
}

void
test_denom_2()
{
  std::cout << "entered test_denom_2" << std::endl;

  sma2::Range r(3, 2);

  sma2::Array<2> A(r, r);

  initialize_dense_array(A);

  std::cout << "A:" << std::endl;
  A.print();

  std::vector<std::vector<double> *> denoms;
  std::vector<double> denoms0, denoms1;
  denoms0.resize(r.nindex());
  denoms1.resize(r.nindex());
  for (int i=0; i<r.nindex(); i++) {
      denoms0[i] = 1.0 + i;
      denoms1[i] = 1.0 + i;
    }
  denoms.push_back(&denoms0);
  denoms.push_back(&denoms1);

  apply_denominator(A,1.0,denoms);

  std::cout << "A/denom:" << std::endl;
  A.print();
}

void
test_all_fixed()
{
  std::cout << "entered test_all_fixed" << std::endl;

  sma2::Range r(3, 1);

  sma2::Array<2> A(r, r), B(r, r);

  initialize_dense_array(A);

  sma2::Index I("i"), J("j");
  B(I,J) |= A(I,J);

  B.zero();

  for (int i=0; i<r.nindex(); i++) {
      for (int j=0; j<r.nindex(); j++) {
          sma2::Index ifix(i), jfix(j);
          B(ifix,jfix) += A(ifix,jfix);
        }
    }

  std::cout << "A:" << std::endl;
  A.print();

  std::cout << "B (should equal A):" << std::endl;
  B.print();
}

void
test_all_bigblock_fixed()
{
  std::cout << "entered test_all_bigblock_fixed" << std::endl;

  sma2::Range r(4, 2);

  sma2::Array<2> A(r, r), B(r, r);

  initialize_dense_array(A);

  sma2::Index I("i"), J("j");
  B(I,J) |= A(I,J);

  B.zero();

  for (int i=0; i<r.nblock(); i++) {
      for (int j=0; j<r.nblock(); j++) {
          sma2::Index ifix("if",i), jfix("jf",j);
          B(ifix,jfix) += A(ifix,jfix);
          if (i==0 && j==0) {
              std::cout << "B should be incomplete (only 0 0 is assigned)"
                        << std::endl;
              B.print();
            }
        }
    }

  std::cout << "A:" << std::endl;
  A.print();

  std::cout << "B (should equal A):" << std::endl;
  B.print();
}

void
test_all_bigblock_fixed_reversed()
{
  std::cout << "entered test_all_bigblock_fixed_reversed" << std::endl;

  sma2::Range r(4, 2);

  sma2::Array<2> A(r, r), B(r, r), C(r, r);

  initialize_dense_array(A);

  sma2::Index I("i"), J("j");
  B(J,I) |= A(I,J);

  B.zero();

  for (int i=0; i<r.nblock(); i++) {
      for (int j=0; j<r.nblock(); j++) {
          sma2::Index ifix("if",i), jfix("jf",j);
          B(jfix,ifix) += A(ifix,jfix);
          if (i==0 && j==0) {
              std::cout << "B should be incomplete (only 0 0 is assigned)"
                        << std::endl;
              B.print();
            }
        }
    }

  std::cout << "B:" << std::endl;
  B.print();

  C(J,I) |= A(I,J);
  C(J,I) = A(I,J);

  std::cout << "C (should equal B):" << std::endl;
  C.print();
}

void
test_contract_bigblock_fixed1()
{
  std::cout << "entered test_contract_bigblock_fixed1" << std::endl;

  sma2::Range r(4, 2);

  sma2::Array<4> C(r, r, r, r);
  sma2::Array<4> C_correct(r, r, r, r);
  sma2::Array<4> A(r, r, r, r);
  sma2::Array<2> B(r, r);

  initialize_dense_array(A);
  initialize_dense_array(B);
  C.add_all_unallocated_blocks();
  C.zero();
  C_correct.add_all_unallocated_blocks();
  C_correct.zero();

  for (int mu = 0; mu<r.nindex(); mu++) {
      for (int rho = 0; rho<r.nindex(); rho++) {
          sma2::Index iMu("mu",mu), iRho("rho",rho);
          sma2::Index iNu("nu"), iJ("j"), iL("l");
          C(iMu,iRho,iJ,iL) += A(iMu,iRho,iNu,iJ) * B(iNu,iL);
          if (mu == 0 && rho == 0) {
              std::cout << "C after mu = 0 rho = 0 set" << std::endl;
              C.print();
            }
        }
    }

  sma2::Index iMu("mu"), iRho("rho");
  sma2::Index iNu("nu"), iJ("j"), iL("l");
  C_correct(iMu,iRho,iJ,iL) += A(iMu,iRho,iNu,iJ) * B(iNu,iL);
  C_correct(iMu,iRho,iJ,iL) -= C(iMu,iRho,iJ,iL);
  std::cout << "C_correct - C:" << std::endl;
  C_correct.print();
}

void
test_contract_bigblock_fixed2()
{
  std::cout << "entered test_contract_bigblock_fixed2" << std::endl;

  sma2::Range r(4, 2);

  sma2::Array<4> C(r, r, r, r);
  sma2::Array<4> C_correct(r, r, r, r);
  sma2::Array<3> A(r, r, r);
  sma2::Array<3> B(r, r, r);

  initialize_dense_array(A);
  initialize_dense_array(B);
  C.add_all_unallocated_blocks();
  C.zero();
  C_correct.add_all_unallocated_blocks();
  C_correct.zero();

  for (int x = 0; x<r.nindex(); x++) {
      for (int y = 0; y<r.nindex(); y++) {
          sma2::Index X("x",x), Y("y",y), I("i"), J("j"), K("k"), L("l");
          C(X,Y,I,K) += A(X, K, J) * B(Y, I, J);
        }
    }

  sma2::Index X("x"), Y("y"), I("i"), J("j"), K("k"), L("l");
  C_correct(X,Y,I,K) += A(X, K, J) * B(Y, I, J);
  C_correct(X,Y,I,K) -= C(X,Y,I,K);
  std::cout << "C_correct - C:" << std::endl;
  C_correct.print();
}

void
test_scalar()
{
  std::cout << "entered test_scalar" << std::endl;

  sma2::Range r(4, 2);

  sma2::Array<2> A(r, r);
  sma2::Array<2> B(r, r);
  sma2::Array<0> C;

  initialize_dense_array(A);
  initialize_dense_array(B);

  C() = A("p","q") * B("p","q");

  double C_scalar_contract = A("p","q") * B("p","q");

  std::cout << "C_scalar_project = " << C_scalar_contract << std::endl;
  std::cout << "C                = " << C.value() << std::endl;
}

void
try_main(int argc, char *argv[])
{
  test_fixed_values();
  test_remapping();
  test_assign();
  test_sum();
  test_div();
  test_contract_2x2x2();
  test_contract_2x2x2_reordered();
  test_contract_empty();
  test_contract_4x6x6();
  test_fixed_contract_Cfixed();
  test_fixed_contract_Bfixed();
  test_fixed_contract_Afixed();
  test_fixed_assign();
  test_fixed_assign_2();
  //test_contract_4x2x2();
  test_fixed_contract_1x2x2();
  test_fixed_contract_3x2x2();
  test_denom_2();
  test_pack_vector();
  test_all_fixed();
  test_all_bigblock_fixed();
  test_all_bigblock_fixed_reversed();
  test_contract_bigblock_fixed1();
  test_contract_bigblock_fixed2();
  test_scalar();

  std::cout << "---- all tests ran to completion ----" << std::endl;
}

}

int
main(int argc, char *argv[])
{
  try {
      sc::try_main(argc, argv);
  }
  catch (std::bad_alloc &e) {
      std::cout << argv[0] << ": ERROR: MEMORY ALLOCATION FAILED:" << std::endl
                << e.what()
                << std::endl;
      throw;
  }
  catch (std::exception &e) {
      std::cout << argv[0] << ": ERROR: EXCEPTION RAISED:" << std::endl
                << e.what()
                << std::endl;
      throw;
  }
  catch (...) {
      std::cout << argv[0] << ": ERROR: UNKNOWN EXCEPTION RAISED" << std::endl;
      throw;
  }
  return 0;
}
