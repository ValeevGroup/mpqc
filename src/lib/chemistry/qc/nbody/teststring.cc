//
// teststring.cc
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

#include "string.h"
#include <iostream>

using namespace sc;

int main(int argc, char **argv) {

  ////////////////////////////////////////////
  // test N-bit strings
  ////////////////////////////////////////////
  {
    //
    // test empty string
    //
    FermionOccupationNBitString<128> os0;
    std::cout << "os0 = " << os0 << std::endl;
    std::cout << "os0 empty? " << (os0.empty() ? "true" : "false")
        << " (should be true)" << std::endl;

    //
    // test non-empty string
    //
    std::vector<FermionOccupationNBitString<128>::state_index_t> sv1;
    sv1.push_back(0);
    sv1.push_back(5);
    sv1.push_back(73);

    FermionOccupationNBitString<128> os1(sv1);
    std::cout << "os1 = " << os1 << std::endl;
    std::cout << "os1 empty? " << (os1.empty() ? "true" : "false")
        << " (should be false)" << std::endl;
    std::cout << "# of occupied states in os1 = " << os1.count()
        << " (should be 3)" << std::endl;

    //
    // apply operators to a non-empty string
    //
    FermionBasicNCOper<1, FermionOccupationNBitString<128> > op1(72, 0);
    std::pair<bool, FermionOccupationNBitString<128> > res1 = op1(os1);
    std::cout << op1 << "(os1) = " << res1.second << " sign="
        << (res1.first ? "-" : "+") << " (should be -)" << std::endl;
  }

  ////////////////////////////////////////////
  // test dynamic bitstrings
  ////////////////////////////////////////////
  {
    //
    // test empty string
    //
    FermionOccupationDBitString os0(127);
    std::cout << "os0 = " << os0 << std::endl;
    std::cout << "os0 empty? " << (os0.empty() ? "true" : "false")
        << " (should be true)" << std::endl;

    //
    // test non-empty string
    //
    std::vector<FermionOccupationDBitString::state_index_t> sv1;
    sv1.push_back(0);
    sv1.push_back(5);
    sv1.push_back(73);

    FermionOccupationDBitString os1(127, sv1);
    std::cout << "os1 = " << os1 << std::endl;
    std::cout << "os1 empty? " << (os1.empty() ? "true" : "false")
        << " (should be false)" << std::endl;
    std::cout << "# of occupied states in os1 = " << os1.count()
        << " (should be 3)" << std::endl;
  }

  ////////////////////////////////////////////
  // test block strings
  ////////////////////////////////////////////
  {
    //
    // test empty string
    //
    FermionOccupationBlockString os0(127);
    std::cout << "os0 = " << os0 << std::endl;
    std::cout << "os0 empty? " << (os0.empty() ? "true" : "false")
        << " (should be true)" << std::endl;

    //
    // test non-empty string
    //
    std::vector<FermionOccupationBlockString::state_index_t> sv1;
    sv1.push_back(0);
    sv1.push_back(1);
    sv1.push_back(2);
    sv1.push_back(5);
    sv1.push_back(7);
    sv1.push_back(6);
    sv1.push_back(73);

    FermionOccupationBlockString os1(127, sv1);
    std::cout << "os1 = " << os1 << std::endl;
    std::cout << "os1 empty? " << (os1.empty() ? "true" : "false")
        << " (should be false)" << std::endl;
    std::cout << "# of occupied states in os1 = " << os1.count()
        << " (should be 7)" << std::endl;
    //
    // apply operators to a non-empty string
    //
    FermionBasicNCOper<1, FermionOccupationBlockString > op1(72, 0);
    std::pair<bool, FermionOccupationBlockString > res1 = op1(os1);
    std::cout << op1 << "(os1) = " << res1.second << " sign="
        << (res1.first ? "-" : "+") << " (should be -)" << std::endl;
  }

  return 0;
}


/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ-CONDENSED"
// End:
