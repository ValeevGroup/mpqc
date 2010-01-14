
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

#include <util/misc/regtime.h>

#include <chemistry/qc/lmp2/sma.h>

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

template <class T>
bool
ready(std::vector<T> &vec,
     const std::vector<T> &start, const std::vector<T> &fence)
{
  if (vec.size() == 0) return false;
  for (int i=0; i<vec.size(); i++) {
      if (vec[i] < start[i]) return false;
      if (vec[i] >= fence[i]) return false;
    }
  return true;
}


template <class T>
void
incr(std::vector<T> &vec,
     const std::vector<T> &start, const std::vector<T> &fence)
{
  if (vec.size() == 0) return;
  for (int i=vec.size()-1; i>0; i--) {
      if (vec[i] >= fence[i]-1) vec[i] = start[i];
      else {
          vec[i]++;
          return;
        }
    }
  vec[0]++;
}

template <int N>
void
initialize_dense_array(sma2::Array<N> &array, double offset = 0.0)
{
  std::vector<sma2::bi_t> b_indices(N);
  std::vector<sma2::bi_t> b_start(N);
  std::vector<sma2::bi_t> b_fence(N);

  std::fill(b_start.begin(), b_start.end(), 0);
  for (int i=0; i<N; i++) b_fence[i] = array.index(i).nblock();

  array.clear();

  for (b_indices=b_start;
       ready(b_indices, b_start, b_fence);
       incr(b_indices, b_start, b_fence)) {
      sma2::BlockInfo<N> bi(b_indices);
      array.add_unallocated_block(bi);
    }

  initialize_array_values(array, offset);
}

void
time_contract_444_plain(sma2::Range&o, sma2::Range&v,
                        sma2::Array<4> &A,
                        sma2::Array<4> &B,
                        const sc::Ref<sc::RegionTimer> &timer)
{

  sma2::Array<4> C(B.index(0),B.index(1),B.index(2),B.index(3));

  timer->enter("blocklist");
  C("i","j","c","d") |= A("c","d","a","b") * B("i","j","a","b");
  timer->change("contract");
  C("i","j","c","d") = (A("c","d","a","b") * B("i","j","a","b")).timer(timer);
  timer->exit();

}

void
time_contract_444_bswap(sma2::Range&o, sma2::Range&v,
                        sma2::Array<4> &A,
                        sma2::Array<4> &B,
                        const sc::Ref<sc::RegionTimer> &timer)
{

  sma2::Array<4> C(B.index(0),B.index(1),B.index(2),B.index(3));
  sma2::Array<4> Bnew(B.index(2),B.index(3),B.index(0),B.index(1));

  timer->enter("reorder B");
  Bnew("i","j","k","l") |= B("k","l","i","j");
  Bnew("i","j","k","l") = B("k","l","i","j");
  timer->change("blocklist");
  C("i","j","c","d") |= A("c","d","a","b") * Bnew("a","b","i","j");
  timer->change("contract");
  C("i","j","c","d") = (A("c","d","a","b") * Bnew("a","b","i","j")).timer(timer);
  timer->exit();

}

void
time_contract_444_swap2(sma2::Range&o, sma2::Range&v,
                        sma2::Array<4> &A,
                        sma2::Array<4> &B,
                        const sc::Ref<sc::RegionTimer> &timer)
{
  // modeled after:
  // Res(iI,iJ,iR,iS) += C(iI,iJ,iT,iU) * K_vir(iR,iT,iU,iS);

  sma2::Array<4> C(B.index(0),B.index(1),B.index(2),B.index(3));

  timer->enter("blocklist");
  C("i","j","r","s") |= B("i","j","t","u") * A("r","t","u","s");
  timer->change("contract");
  C("i","j","r","s") = (B("i","j","t","u") * A("r","t","u","s")).timer(timer);
  timer->exit();

}

void
time_contract_444_swap3(sma2::Range&o, sma2::Range&v,
                        sma2::Array<4> &A,
                        sma2::Array<4> &B,
                        const sc::Ref<sc::RegionTimer> &timer)
{
  // modeled after:
  // Res(iI,iJ,iR,iS) += C(iI,iJ,iT,iU) * K_vir(iR,iT,iU,iS);

  sma2::Array<4> C(B.index(0),B.index(1),B.index(2),B.index(3));

  timer->enter("blocklist");
  C("i","j","r","s") |= A("r","t","u","s") * B("i","j","t","u");
  timer->change("contract");
  C("i","j","r","s") = (A("r","t","u","s") * B("i","j","t","u")).timer(timer);
  timer->exit();

}

void
time_contract_444_fixed_1(sma2::Range&o, sma2::Range&v,
                          sma2::Array<4> &A,
                          sma2::Array<4> &B,
                          const sc::Ref<sc::RegionTimer> &timer)
{
  // arguments are
  // A(v,v,v,v)
  // B(o,o,v,v)

  // C(v,v)
  sma2::Array<2> C(B.index(2),B.index(3));

  sma2::Index iA("a"), iB("b"), iC("c"), iD("d");

  timer->enter("loops");
  for (int k=0; k<B.index(0).nindex(); k++) {
      sma2::Index iK(k);
      for (int l=0; l<B.index(1).nindex(); l++) {
          sma2::Index iL(l);
          C.clear();
          timer->enter("blocklist");
          C(iA,iB) |= A(iA,iB,iC,iD) * B(iK,iL,iC,iD);
          timer->change("contract");
          C(iA,iB) = (A(iA,iB,iC,iD) * B(iK,iL,iC,iD)).timer(timer);
          timer->exit();
        }
    }
  timer->exit();
}

void
time_contract_444_fixed_2(sma2::Range&o, sma2::Range&v,
                          sma2::Array<4> &A,
                          sma2::Array<4> &B,
                          const sc::Ref<sc::RegionTimer> &timer)
{
  // arguments are
  // A(v,v,v,v)
  // B(o,o,v,v)

  // C(v,v)
  sma2::Array<2> C(B.index(2),B.index(3));

  sma2::Index iA("a"), iB("b"), iC("c"), iD("d");

  timer->enter("loops");
  for (int k=0; k<B.index(0).nindex(); k++) {
      sma2::Index iK(k);
      for (int l=0; l<B.index(1).nindex(); l++) {
          sma2::Index iL(l);
          timer->enter("blocklist");
          C.clear();
          C(iA,iB) |= A(iC,iD,iA,iB) * B(iK,iL,iC,iD);
          timer->change("contract");
          C(iA,iB) = (A(iC,iD,iA,iB) * B(iK,iL,iC,iD)).timer(timer);
          timer->exit();
        }
    }
  timer->exit();
}

void
time_contract_444_fixed_3(sma2::Range&o, sma2::Range&v,
                          sma2::Array<4> &A,
                          sma2::Array<4> &B,
                          const sc::Ref<sc::RegionTimer> &timer)
{
  // arguments are
  // A(v,v,v,v)
  // B(o,o,v,v)

  // C(v,v)
  sma2::Array<2> C(B.index(2),B.index(3));

  sma2::Index iA("a"), iB("b"), iC("c"), iD("d");

  timer->enter("loops");
  for (int k=0; k<B.index(0).nindex(); k++) {
      sma2::Index iK(k);
      for (int l=0; l<B.index(1).nindex(); l++) {
          sma2::Index iL(l);
          C.clear();
          timer->enter("blocklist");
          C(iA,iB) |= B(iK,iL,iC,iD) * A(iC,iD,iA,iB);
          timer->change("contract");
          C(iA,iB) = (B(iK,iL,iC,iD) * A(iC,iD,iA,iB)).timer(timer);
          timer->exit();
        }
    }
  timer->exit();
}

void
time_contract_444()
{
  std::cout << "entered time_contract_444" << std::endl;

  sc::Ref<sc::RegionTimer> timer = new sc::RegionTimer("total", 1, 1);

  // set up for sto-3g uracil: 4 h and 8 c/n/o
//   sma2::Range o(29,1);
//   std::vector<int> blocks;
//   for (int i=0; i<8; i++) blocks.push_back(5);
//   for (int i=0; i<4; i++) blocks.push_back(1);
//   sma2::Range v(blocks);

  // set up for 3-21g uracil: 4 h and 8 c/n/o
//   sma2::Range o(29,1);
//   std::vector<int> blocks;
//   for (int i=0; i<8; i++) blocks.push_back(9);
//   for (int i=0; i<4; i++) blocks.push_back(2);
//   sma2::Range v(blocks);

  // set up for CH3F 6-31G**
  sma2::Range o(9,1);
  std::vector<int> blocks;
  // C, F
  for (int i=0; i<2; i++) blocks.push_back(15);
  // 3 H
  for (int i=0; i<3; i++) blocks.push_back(5);
  sma2::Range v(blocks);

  sma2::Array<4> A(v,v,v,v),B(o,o,v,v);

  initialize_dense_array(A);
  initialize_dense_array(B);
  
  std::cout << "  plain" << std::endl;
  timer->enter("plain");
  time_contract_444_plain(o,v,A,B,timer);
  timer->exit();

  std::cout << "  bswap" << std::endl;
  timer->enter("bswap");
  time_contract_444_bswap(o,v,A,B,timer);
  timer->exit();

  std::cout << "  swap2" << std::endl;
  timer->enter("swap2");
  time_contract_444_swap2(o,v,A,B,timer);
  timer->exit();

  std::cout << "  swap3" << std::endl;
  timer->enter("swap3");
  time_contract_444_swap3(o,v,A,B,timer);
  timer->exit();
  
  std::cout << "  fixed 1" << std::endl;
  timer->enter("fixed 1");
  time_contract_444_fixed_1(o,v,A,B,timer);
  timer->exit();
  
  std::cout << "  fixed 2" << std::endl;
  timer->enter("fixed 2");
  time_contract_444_fixed_2(o,v,A,B,timer);
  timer->exit();
  
  std::cout << "  fixed 3" << std::endl;
  timer->enter("fixed 3");
  time_contract_444_fixed_3(o,v,A,B,timer);
  timer->exit();

//   /////////////////////////////////////////////////

//   sma2::Range o2(o.nindex(),2);
//   sma2::Range v2(v.nindex(),8);

//   sma2::Array<4> A2(v2,v2,v2,v2),B2(v2,v2,o2,o2);

//   initialize_dense_array(A2);
//   initialize_dense_array(B2);

//   std::cout << "  plain2" << std::endl;
//   timer->enter("plain2");
//   time_contract_444_plain(o2,v2,A2,B2,timer);
//   timer->exit();

//   /////////////////////////////////////////////////

//   sma2::Range os(o.nindex(),o.nindex());
//   sma2::Range vs(v.nindex(),v.nindex());

//   sma2::Array<4> As(vs,vs,vs,vs),Bs(vs,vs,os,os);

//   initialize_dense_array(As);
//   initialize_dense_array(Bs);

//   std::cout << "  plain single block" << std::endl;
//   timer->enter("plain single block");
//   time_contract_444_plain(os,vs,As,Bs,timer);
//   timer->exit();

  timer->print();

}

void
try_main(int argc, char *argv[])
{
  time_contract_444();
}

int
main(int argc, char *argv[])
{
  try {
      try_main(argc, argv);
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

}
