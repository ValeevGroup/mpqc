//
// repltest.cc
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
#include <util/group/messshm.h>
#include <math/scmat/repl.h>

ClassDesc* f0 = &ShmMessageGrp::class_desc_;
#ifdef HAVE_MPI
#include <util/group/messmpi.h>
ClassDesc* f1 = &MPIMessageGrp::class_desc_;
#endif

void matrixtest(RefSCMatrixKit kit, RefKeyVal keyval,
                RefSCDimension d1,RefSCDimension d2,RefSCDimension d3);

main(int argc, char** argv)
{
  RefKeyVal keyval = new ParsedKeyVal(SRCDIR "/matrixtest.in");

  RefMessageGrp msg = MessageGrp::initial_messagegrp(argc, argv);

  if (msg.null()) {
      msg = keyval->describedclassvalue("messagegrp");

      if (msg.null()) {
          cerr << indent << "Couldn't initialize MessageGrp\n";
          abort();
        }
    }

  MessageGrp::set_default_messagegrp(msg);

  RefSCMatrixKit kit = new ReplSCMatrixKit;
  RefSCDimension d1(keyval->describedclassvalue("d1"));
  RefSCDimension d2(keyval->describedclassvalue("d2"));
  RefSCDimension d3(keyval->describedclassvalue("d3"));

  matrixtest(kit,keyval,d1,d2,d3);

  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// eval: (c-set-style "CLJ")
// End:
