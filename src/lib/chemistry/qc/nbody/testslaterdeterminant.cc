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
#include "string.h"
#include <iostream>

using namespace sc;

template <typename RandomAccessStringIterator>
class SlaterDeterminantSet {

  public:

    typedef typename RandomAccessStringIterator::value_type String;

    struct SlaterDeterminant {
      public:
        SlaterDeterminant(const RandomAccessStringIterator& as,
                          const RandomAccessStringIterator& bs) : astr(as), bstr(bs) {}

        bool operator<(const SlaterDeterminant& other) const {
          if (astr == other.astr)
            return bstr < other.bstr;
          else
            return astr < other.astr;
        }

        /// alpha string
        RandomAccessStringIterator astr;
        /// beta string
        RandomAccessStringIterator bstr;
    };

    SlaterDeterminantSet(RandomAccessStringIterator abegin,
                         RandomAccessStringIterator aend,
                         RandomAccessStringIterator bbegin,
                         RandomAccessStringIterator bend,
                         bool full_ci = false) :
                         abegin_(abegin),
                         aend_(aend),
                         bbegin_(bbegin),
                         bend_(bend),
                         num_bstr_(bend-bbegin),
                         full_ci_(full_ci)
                         {
                         }

    size_t add(RandomAccessStringIterator astr,
               RandomAccessStringIterator bstr) {
      const size_t num_sd = sdset_.size();
      sdset_[SlaterDeterminant(astr,bstr)] = num_sd;
    }
    size_t safe_add(RandomAccessStringIterator astr,
                    RandomAccessStringIterator bstr) {
      SlaterDeterminant sd(astr,bstr);
      if (sdset_.find(sd) == sdset_.end()) {
        const size_t num_sd = sdset_.size();
        sdset_[SlaterDeterminant(astr,bstr)] = num_sd;
      }
    }

    long long find(RandomAccessStringIterator astr,
                   RandomAccessStringIterator bstr) const {
      if (full_ci_) {
        return (astr - abegin_) * num_bstr_ + (bstr - bbegin_);
      }
      else {
        typename SDContainer::const_iterator result_iter = sdset_.find(SlaterDeterminant(astr,bstr));
        return (*result_iter).second;
      }
    }

    long long safe_find(RandomAccessStringIterator astr,
                        RandomAccessStringIterator bstr) const {
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
    RandomAccessStringIterator abegin_;
    RandomAccessStringIterator aend_;
    RandomAccessStringIterator bbegin_;
    RandomAccessStringIterator bend_;
    size_t num_bstr_;
    bool full_ci_;
    typedef std::map<SlaterDeterminant, size_t> SDContainer;
    SDContainer sdset_;
};

int main(int argc, char **argv) {

  typedef FermionOccupationDBitString FString;
  //typedef FermionOccupationBlockString FString;
  typedef std::vector<FString> stringset_t;
  typedef stringset_t::const_iterator stringsetiter_t;

  stringset_t strings;
  {
    std::vector<FString::state_index_t> sv1;
    sv1.push_back(0);
    sv1.push_back(1);

    strings.push_back( FString(4, sv1) );

    sv1[1] = 2;
    strings.push_back( FString(4, sv1) );
    sv1[1] = 3;
    strings.push_back( FString(4, sv1) );
    sv1[0] = 2;
    strings.push_back( FString(4, sv1) );

  }

  SlaterDeterminantSet<stringsetiter_t> sdset(strings.begin(), strings.end(), strings.begin(), strings.end(), true);
  for(stringsetiter_t ai=strings.begin(); ai!=strings.end(); ++ai) {
    const auto astr_ex_lvl = ((*ai) ^ strings.front()).count() / 2;
    for(stringsetiter_t bi=strings.begin(); bi!=strings.end(); ++bi) {
      const auto bstr_ex_lvl = ((*bi) ^ strings.front()).count() / 2;
      auto sdi = sdset.find(ai, bi);
      std::cout << "astr = " << *ai << " bstr = " << *bi << " sd index = " << sdi << " ex_lvl = " << astr_ex_lvl+bstr_ex_lvl << std::endl;
    }
  }

  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
