//
// repltest.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
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
#include <util/group/pregtime.h>
#include <math/scmat/repl.h>

#include <util/group/linkage.h>

using namespace sc;

void matrixtest(Ref<SCMatrixKit> kit, Ref<KeyVal> keyval,
                RefSCDimension d1,RefSCDimension d2,RefSCDimension d3,
                bool have_svd);

int
main(int argc, char** argv)
{
  char *infile = SRCDIR "/matrixtest.in";
  
  if (argc > 1)
      infile = argv[1];

  Ref<KeyVal> keyval = new ParsedKeyVal(infile);

  Ref<MessageGrp> msg = MessageGrp::initial_messagegrp(argc, argv);

  if (msg == 0) {
      msg << keyval->describedclassvalue("messagegrp");

      if (msg == 0) {
          std::cerr << indent << "Couldn't initialize MessageGrp\n";
          abort();
        }
    }

  MessageGrp::set_default_messagegrp(msg);
  Ref<RegionTimer> tim = new ParallelRegionTimer(msg,"matrixtest",1,1);
  RegionTimer::set_default_regiontimer(tim);

  SCFormIO::set_printnode(0);
  if (msg->n() > 1)
    SCFormIO::init_mp(msg->me());

  Ref<SCMatrixKit> kit = new ReplSCMatrixKit;
  RefSCDimension d1; d1 << keyval->describedclassvalue("d1");
  RefSCDimension d2; d2 << keyval->describedclassvalue("d2");
  RefSCDimension d3; d3 << keyval->describedclassvalue("d3");

  matrixtest(kit,keyval,d1,d2,d3,true);

  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
