//
// kvopt.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
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
#include <math/optimize/opt.h>
#include <util/keyval/keyval.h>

int
main(int argc, char** argv)
{
  if (argc != 3) {
      cerr << scprintf("Usage: %s inputfile keyword\n", argv[0]);
      exit(1);
    }
  char* inputfile = argv[1];
  char* keyword = argv[2];

  FILE*fp = fopen(inputfile,"r");
  if (!fp) {
      cerr << scprintf("%s: error opening input file \"%s\"\n",
                       argv[0], inputfile);
      perror("fopen");
      exit(1);
    }

  RefParsedKeyVal keyval = new ParsedKeyVal(fp);

  if (!keyval->exists(keyword)) {
      cerr << scprintf("%s: keyword \"%s\" doesn't exist in file \"%s\"\n",
                       argv[0], keyword, inputfile);
      exit(1);
    }

  if (!keyval->classname(keyword)) {
      cerr << scprintf("%s: keyword \"%s\" in file \"%s\" is not an object\n",
                       argv[0], keyword, inputfile);
      exit(1);
    }

  RefDescribedClass dc = keyval->describedclassvalue(keyword);

  if (dc.null()) {
      cerr << scprintf("%s: keyword \"%s\" in file \"%s\" could not be"
                       " converted into a DescribedClass object\n",
                       argv[0], keyword, inputfile);
      exit(1);
    }

  RefOptimize opt = dc;

  if (opt.null()) {
      cerr << scprintf("%s: keyword \"%s\" in file \"%s\" could not be"
                       " converted into an Optimize object\n",
                       argv[0], keyword, inputfile);
      exit(1);
    }

  opt->optimize();

  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
// End:
