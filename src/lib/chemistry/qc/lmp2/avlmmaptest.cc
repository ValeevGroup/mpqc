
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

#include <algorithm>
#include <vector>

#include <util/misc/regtime.h>

#include <chemistry/qc/lmp2/avlmmap.h>

namespace sc {

template <class K, class T>
int
count(const AVLMMap<K,T> &m, K n)
{
  typename AVLMMap<K,T>::const_iterator begin = m.lower_bound(n);
  typename AVLMMap<K,T>::const_iterator end = m.upper_bound(n);
  int nmatch = 0;
  bool bad = false;
  for (typename AVLMMap<K,T>::const_iterator i = begin; i != end; i++) {
      if (i == m.end()) {
          std::cout << "count: i is m.end()" << std::endl;
          std::cout << "searching for " << n << std::endl;
          bad = true;
          break;
        }
//       std::cout << " element key: " << i->first
//                 << " search key: " << n
//                 << std::endl;
      if (n != i->first) {
          if (!bad) std::cout << "count: mismatched data" << std::endl;
          bad = true;
        }
      nmatch++;
    }
  if (bad) {
      if (begin != m.end()) std::cout << "lb = " << begin->first << std::endl;
      if (end != m.end()) std::cout << "ub = " << end->first << std::endl;
      m.print();
      abort();
    }
  return nmatch;
}

void
const_tests(const AVLMMap<int,int> &m)
{
  AVLMMap<int,int>::const_iterator iter;

  iter = m.find(4);
  std::cout << " 4: " << iter->first << " " << iter->second << std::endl;

}

void
nonconst_tests(AVLMMap<int,int> &m)
{
  AVLMMap<int,int>::iterator iter;

  iter = m.find(4);
  std::cout << " 4: " << iter->first << " " << iter->second << std::endl;
}

void
print_iter_info(AVLMMap<int,double>::const_iterator lb, const char *name,
                AVLMMap<int,double>::const_iterator end)
{
  std::cout << name << ": ";
  if (lb != end) {
      std::cout << lb->first;
    }
  else {
      std::cout << "end";
    }
  std::cout << std::endl;
}

void
random_test_1(const std::vector<int> &v, int n)
{
  AVLMMap<int,double> m;
  int n0 = 0;
  int n99 = 0;
  AVLMMap<int,double>::iterator iter = m.end();
  for (int i=0; i<n; i++) {
      int number_added = 1;
      // try different insert methods
      if (i%5 == 0) {
          if (m.find(v[i]) != m.end()) {
              number_added = 0;
              iter = m.insert_unique(std::make_pair(v[i],i));
            }
          else if (i%2 == 0) {
              iter = m.insert_unique(std::make_pair(v[i],i));
            }
          else {
              iter = m.insert_new(std::make_pair(v[i],i));
            }
        }
      else if (i%2) iter = m.insert(std::make_pair(v[i],i));
      else iter = m.insert(iter, std::make_pair(v[i],i));
      if (v[i] == 0) n0 += number_added;
      else if (v[i] == 99) n99 += number_added;
    }
  if (n0 != count(m,0)) {
      std::cout << "miscounted n0" << std::endl;
    }
  if (n99 != count(m,99)) {
      std::cout << "miscounted n99" << std::endl;
    }
  for (int i=-1; i<=n; i++) {
      ////////////////
      // find with hint tests
      for (AVLMMap<int,double>::const_iterator iter = m.begin();
           iter != m.end();
           iter++) {
          AVLMMap<int,double>::const_iterator f1 = m.find(iter, i);
          AVLMMap<int,double>::const_iterator f2 = m.find(i);
          if (f1 != f2
              && (f1 == m.end()
                  || f2 == m.end()
                  || m.key_comp()(f1->first, f2->first)
                  || m.key_comp()(f2->first, f1->first))
              ) {
              std::cout << "bad hint find" << std::endl;
              abort();
            }
        }
      // lower_bound/upper_bound tests
      AVLMMap<int,double>::const_iterator lb_check = m.end();
      AVLMMap<int,double>::const_iterator lb = m.lower_bound(i);
      AVLMMap<int,double>::const_iterator ub_check = m.end();
      AVLMMap<int,double>::const_iterator ub = m.upper_bound(i);
      std::pair<
          AVLMMap<int,double>::const_iterator,
          AVLMMap<int,double>::const_iterator
          > er = m.equal_range(i);
      AVLMMap<int,double>::const_iterator er_lb = er.first;
      AVLMMap<int,double>::const_iterator er_ub = er.second;
      //std::cout << "---" << std::endl;
      for (AVLMMap<int,double>::const_iterator iter = m.begin();
           iter != m.end();
           iter++) {
          //std::cout << "processing " << iter->first << std::endl;
          if (m.key_comp()(i, iter->first)) {
              // iter->first > i
              if (lb_check == m.end()) lb_check = iter;
              if (ub_check == m.end()) ub_check = iter;
              break;
            }
          else if (!m.key_comp()(iter->first,i)) {
              // iter->first == i
              if (lb_check == m.end()) lb_check = iter;
            }
          else {
              // iter->first < i
            }
        }
      bool bad = false;
      if (lb != lb_check) {
          std::cout << "lower_bound problem" << std::endl;
          bad = true;
        }
      if (ub != ub_check) {
          std::cout << "upper_bound problem" << std::endl;
          bad = true;
        }
      if (lb != er_lb) {
          std::cout << "equal_range lower_bound problem" << std::endl;
          bad = true;
        }
      if (ub != er_ub) {
          std::cout << "equal_range upper_bound problem" << std::endl;
          bad = true;
        }
      if (bad) {
          std::cout << "search key: " << i << std::endl;
          print_iter_info(lb,"lb",m.end());
          print_iter_info(lb_check,"lb_check",m.end());
          print_iter_info(ub,"ub",m.end());
          print_iter_info(ub_check,"ub_check",m.end());
          m.print();
          abort();
        }
    }
}

void
random_tests()
{
  int nmax = 50;
  std::vector<int> v(nmax);
  for (int i=0; i<nmax; i++) {
      v[i] = i;
    }
  // make some entries duplicated
  for (int i=0; i<10; i++) {
      v[i] = 0;
      v[i+10] = 99;
    }
  for (int i=0; i<=nmax; i++) {
      for (int j=0; j<=i*3; j++) {
          std::random_shuffle(v.begin(), v.end());
          random_test_1(v, i);
        }
    }
}

void
mmap_timings()
{
  const int nrepeat = 3;
  const int nelement = 10000000;
//   const int nrepeat = 2;
//   const int nelement = 600;

  std::cout << "node size = " << sizeof(AVLMMapNode<int,int>) << std::endl;

  sc::Ref<sc::RegionTimer> timer = new sc::RegionTimer;
  AVLMMap<int,int> m;
  timer->enter("insert_equal nohint");
  for (int i=0; i<nrepeat; i++) {
      for (int j=0; j<nelement; j++) {
          m.insert_equal(std::pair<int,int>(j,0));
        }
      m.clear();
    }
  timer->exit();
  timer->enter("insert_unique nohint");
  for (int i=0; i<nrepeat; i++) {
      for (int j=0; j<nelement; j++) {
          m.insert_unique(std::pair<int,int>(j,0));
        }
      m.clear();
    }
  timer->exit();
  timer->enter("insert_new nohint");
  for (int i=0; i<nrepeat; i++) {
      for (int j=0; j<nelement; j++) {
          m.insert_new(std::pair<int,int>(j,0));
        }
      m.clear();
    }
  timer->exit();
  timer->enter("insert_new hint");
  for (int i=0; i<nrepeat; i++) {
      AVLMMap<int,int>::iterator hint = m.end();
      for (int j=0; j<nelement; j++) {
          hint = m.insert_new(hint, std::pair<int,int>(j,0));
        }
      m.clear();
    }
  timer->exit();
  timer->enter("insert_new bad hint");
  for (int i=0; i<nrepeat; i++) {
      for (int j=0; j<nelement; j++) {
          m.insert_new(m.begin(), std::pair<int,int>(j,0));
        }
      m.clear();
    }
  timer->exit();
  timer->print();
}

}

using namespace sc;

int
main(int argc, char* argv[])
{
  mmap_timings();
  
  AVLMMap<int,int> m;
  for (int i=0; i<10; i++) {
      m.insert(std::make_pair(i,i));
    }
  for (int i=0; i<5; i++) {
      m.insert(std::make_pair(5,6+i));
    }
  for (AVLMMap<int,int>::iterator i = m.begin();
       i != m.end();
       i++) {
      std::cout << " " << i->first << " " << i->second
                << " n=" << count(m,i->first)
                << std::endl;
    }
  const_tests(m);
  nonconst_tests(m);
  random_tests();
  return 0;
}
