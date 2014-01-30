//
// tawfn_test.cpp
//
// Copyright (C) 2013 Drew Lewis
//
// Authors: Drew Lewis
// Maintainer: Drew Lewis and Edward Valeev
//
// This file is part of the MPQC Toolkit.
//
// The MPQC Toolkit is free software; you can redistribute it and/or modify
// it under the terms of the GNU Library General Public License as published by
// the Free Software Foundation; either version 2, or (at your option)
// any later version.
//
// The MPQC Toolkit is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Library General Public License for more details.
//
// You should have received a copy of the GNU Library General Public License
// along with the MPQC Toolkit; see the file COPYING.LIB.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#include <mpqc/tawfn/mock_wfn.hpp>
#include <iostream>
#define BOOST_TEST_MODULE test_mock_wfn
#include <boost/test/included/unit_test.hpp>

using namespace boost::unit_test;
using namespace mpqc;
BOOST_AUTO_TEST_CASE( test_mock_wfn ){
    BOOST_MESSAGE("Checking value in vector");
    std::vector<double> a(10,5);
    for(auto &elem : a){
        BOOST_CHECK_EQUAL(5, elem);
    }
}
