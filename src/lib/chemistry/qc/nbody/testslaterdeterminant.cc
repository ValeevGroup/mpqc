//
// testslaterdeterminant.cc
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

#ifdef __GNUG__
#pragma implementation
#endif

#include <map>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <time.h>
#include <cxxabi.h>
#include "string.h"
#include <cassert>

using namespace sc;

/**
 * Implements a "dense" set of Slater determinants
 * @tparam stringset_type container that provides random-access iterators
 */
template <typename stringset_type>
class SlaterDeterminantDenseSet {

  public:
    typedef typename stringset_type::const_iterator string_iterator;
    typedef typename string_iterator::value_type String;

    struct SlaterDeterminant {
      public:
        SlaterDeterminant(const string_iterator& as,
                          const string_iterator& bs) : astr(as), bstr(bs) {}

        bool operator<(const SlaterDeterminant& other) const {
          if (astr == other.astr)
            return bstr < other.bstr;
          else
            return astr < other.astr;
        }

        /// alpha string
        string_iterator astr;
        /// beta string
        string_iterator bstr;
    };

    SlaterDeterminantDenseSet(string_iterator abegin,
                         string_iterator aend,
                         string_iterator bbegin,
                         string_iterator bend,
                         bool full_ci = false) :
                         abegin_(abegin),
                         aend_(aend),
                         bbegin_(bbegin),
                         bend_(bend),
                         num_bstr_(bend-bbegin),
                         full_ci_(full_ci)
                         {
                         }

    size_t add(string_iterator astr,
               string_iterator bstr) {
      const size_t num_sd = sdset_.size();
      sdset_[SlaterDeterminant(astr,bstr)] = num_sd;
    }
    size_t safe_add(string_iterator astr,
                    string_iterator bstr) {
      SlaterDeterminant sd(astr,bstr);
      if (sdset_.find(sd) == sdset_.end()) {
        const size_t num_sd = sdset_.size();
        sdset_[SlaterDeterminant(astr,bstr)] = num_sd;
      }
    }

    long long find(string_iterator astr,
                   string_iterator bstr) const {
      if (full_ci_) {
        return (astr - abegin_) * num_bstr_ + (bstr - bbegin_);
      }
      else {
        typename SDContainer::const_iterator result_iter = sdset_.find(SlaterDeterminant(astr,bstr));
        return (*result_iter).second;
      }
    }

    long long safe_find(string_iterator astr,
                        string_iterator bstr) const {
      if (full_ci_) {
        return (astr - abegin_) * num_bstr_ + (bstr - bbegin_);
      }
      else {
        typename SDContainer::const_iterator result_iter = sdset_.find(SlaterDeterminant(astr,bstr));
        if (result_iter != sdset_.end())
          return (*result_iter).second;
        else
          return -1;
      }
    }

    SlaterDeterminant operator[](size_t i) const {
      if (full_ci_) {
        const size_t a = i / num_bstr_;
        const size_t b = i % num_bstr_;
        return SlaterDeterminant(abegin_ + a, bbegin_ + b);
      }
      else
        throw std::logic_error("SlaterDeterminantSet::operator[] only usable in full CI");
    }

  private:
    string_iterator abegin_;
    string_iterator aend_;
    string_iterator bbegin_;
    string_iterator bend_;
    size_t num_bstr_;
    bool full_ci_;
    typedef std::map<SlaterDeterminant, size_t> SDContainer;
    SDContainer sdset_;
};

/**
 * Implements "sparse" set of Slater determinants.
 * @tparam stringset_type container that provides forward iterators
 */
template <typename stringset_type>
class SlaterDeterminantSparseSet {

  public:
    typedef typename stringset_type::const_iterator string_iterator;
    typedef typename stringset_type::value_type String;

    struct SlaterDeterminant {
      public:
        SlaterDeterminant(string_iterator as,
                          string_iterator bs) : str(std::make_pair(as,bs)) {}

        bool operator==(const SlaterDeterminant& other) const {
          return str == other.str;
        }

        size_t hash_value() const {
          boost::hash< std::pair<const String*,const String*> > h;
          return h(std::make_pair(&(*str.first),&(*str.second)));
        }

        struct hash {
            size_t operator()(const SlaterDeterminant& d) const {
              return d.hash_value();
            }
        };

        /// {alpha,beta} string
        std::pair<string_iterator,string_iterator> str;
    };

    SlaterDeterminantSparseSet(const stringset_type& astrings,
                               const stringset_type& bstrings) :
                                 astrings_(astrings),
                                 bstrings_(bstrings)
                         {
                         }

    size_t insert(string_iterator astr,
                  string_iterator bstr) {
      const size_t num_sd = sdset_.size();
      sdset_[SlaterDeterminant(astr,bstr)] = num_sd;
    }
    size_t safe_insert(string_iterator astr,
                      string_iterator bstr) {
      SlaterDeterminant sd(astr,bstr);
      if (sdset_.find(sd) == sdset_.end()) {
        const size_t num_sd = sdset_.size();
        sdset_[SlaterDeterminant(astr,bstr)] = num_sd;
      }
    }

    long long find(string_iterator astr,
                   string_iterator bstr) const {
      auto result_iter = sdset_.find(SlaterDeterminant(astr,bstr));
      return (*result_iter).second;
    }

    long long safe_find(string_iterator astr,
                        string_iterator bstr) const {
      auto result_iter = sdset_.find(SlaterDeterminant(astr,bstr));
      if (result_iter != sdset_.end())
        return (*result_iter).second;
      else
        return -1;
    }

  private:
    const stringset_type& astrings_;
    const stringset_type& bstrings_;
    typedef std::unordered_map<SlaterDeterminant, size_t, typename SlaterDeterminant::hash> SDContainer;
    SDContainer sdset_;
};

int try_main(int argc, char **argv);
template <typename FString>
  int test1(size_t nstates, size_t nparticles);

int main(int argc, char **argv) {

  try {
    try_main(argc, argv);
  }
  catch(...) {
    std::cerr << "caught an exception, bye-bye" << std::endl;
    return 1;
  }

  return 0;
}

int try_main(int argc, char **argv) {

  size_t nstates, nparticles;

  for(int occtype=0; occtype<2; ++occtype) {
    switch (occtype) {
      case 0:
        // low occupancy
        nstates = 50;
        nparticles = 5;
        break;

      case 1:
        // half-occupancy
        nstates = 24;
        nparticles = 12;
        break;

      default:
        MPQC_ASSERT(false);
        break;
    }

    test1<FermionOccupationNBitString<64> >(nstates, nparticles);
    test1<FermionOccupationNBitString<128> >(nstates, nparticles);
    test1<FermionOccupationDBitString>(nstates, nparticles);
    test1<FermionOccupationBlockString>(nstates, nparticles);
  }

  return 0;
}

template <typename FString>
int test1(size_t nstates, size_t nparticles) {

  typedef FermionStringSparseSet<FString> stringset_type;
  typedef typename stringset_type::const_iterator stringsetiter_type;

  ///////// generate a full set of strings //////////
  clock_t start0 = clock();
  int status;
  char* realname = abi::__cxa_demangle(typeid(FString).name(), 0, 0, &status);
  std::cout << "generating a (m=" << nstates << ",n=" << nparticles << ") string set (type=" << realname << ")" << std::endl;
  stringset_type sset;
  FullFermionStringSetBuild<stringset_type> builder(nstates, nparticles);
  builder(sset);
  std::cout << "# of strings = " << sset.size() << std::endl;
  clock_t stop0 = clock();
  std::cout << "elapsed time = " << (stop0-start0) * 1.0/ CLOCKS_PER_SEC << " sec"<< std::endl;

  ///////// generate replacement lists for the first few //////////
  clock_t start1 = clock();
  const size_t rank = 1;
  const size_t nstrings = 1;
  auto current_string = sset.begin();
  for(size_t s=0; s<nstrings; ++s, ++current_string) {
    std::cout << "iterating over rank-" << rank << " replacement list for string " << *current_string << std::endl;
    StringReplacementListIterator<FString, rank, true> repl_list_iter(*current_string);
  }
  clock_t stop1 = clock();
  std::cout << "elapsed time = " << (stop0-start0) * 1.0/ CLOCKS_PER_SEC << " sec"<< std::endl;

}

void test0() {

  //typedef FermionOccupationDBitString FString;
  typedef FermionOccupationBlockString FString;
#if 0
  if (use_sparse_sdset) {
    typedef FermionStringSparseSet<FString> stringset_type;
    typedef stringset_type::const_iterator stringsetiter_type;

  stringset_type strings;
  std::vector<FString::state_index_type> sv1;
  sv1.push_back(0);
  sv1.push_back(1);

  strings.insert( FString(4, sv1) );
  const FString a_ref(4, sv1);
  const FString b_ref(4, sv1);

  sv1[1] = 2;
  strings.insert( FString(4, sv1) );
  sv1[1] = 3;
  strings.insert( FString(4, sv1) );
  sv1[0] = 2;
  strings.insert( FString(4, sv1) );

  SlaterDeterminantSparseSet<stringset_type> sdset(strings,
                                                   strings);

  for(auto a=strings.begin(); a!=strings.end(); ++a) {
    for(auto b=strings.begin(); b!=strings.end(); ++b) {
      sdset.insert(a,b);
    }
  }

#else
  typedef FermionStringDenseSet<FString> stringset_type;
  typedef stringset_type::const_iterator stringsetiter_type;

  stringset_type strings;
  std::vector<FString::state_index_type> sv1;
  sv1.push_back(0);
  sv1.push_back(1);

  strings.insert(FString(4, sv1));
  const FString a_ref(4, sv1);
  const FString b_ref(4, sv1);

  sv1[1] = 2;
  strings.insert(FString(4, sv1));
  sv1[1] = 3;
  strings.insert(FString(4, sv1));
  sv1[0] = 2;
  strings.insert(FString(4, sv1));

  SlaterDeterminantDenseSet<stringset_type> sdset(strings.begin(),
                                                  strings.end(),
                                                  strings.begin(),
                                                  strings.end());

#endif

  for(auto ai=strings.begin(); ai!=strings.end(); ++ai) {
    const auto astr_ex_lvl = ((*ai) ^ a_ref ).count() / 2;
    for(auto bi=strings.begin(); bi!=strings.end(); ++bi) {
      const auto bstr_ex_lvl = ((*bi) ^ b_ref).count() / 2;
      auto sdi = sdset.find(ai, bi);
      std::cout << "astr = " << *ai << " bstr = " << *bi << " sd index = " << sdi << " ex_lvl = " << astr_ex_lvl+bstr_ex_lvl << std::endl;
    }
  }

}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
