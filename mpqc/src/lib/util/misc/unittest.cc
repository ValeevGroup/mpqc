//
// unittest.cc
//
// Copyright (C) 1998 Sandia National Laboratories 
//
// Author: Ida Nielsen <ibniels@idap2.ca.sandia.gov>
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

#include <util/misc/units.h>
#include <util/misc/formio.h>

//////////////////////////////////////////////////////////////////////
// Unit conversion test program 

int main(int argc, char **argv) {
  if (argc != 2) {
      cerr << "One argument, the unit to be converted, must be given" << endl;
      abort();
    }
  
  RefUnits unit = new Units(argv[1]);
  cout << indent << "Conversion between " << unit->string_rep() << " and atomic units:" << endl;
  cout << setprecision(10);
  cout << indent << "From atomic units: " << unit->from_atomic_units() << endl;
  cout << indent << "To atomic units: " << unit->to_atomic_units() << endl;

  return 0;

}


// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
