//
// disttest.cc
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

#include <iostream.h>
#include <math.h>

#include <util/misc/bug.h>
#include <util/misc/formio.h>
#include <util/group/messshm.h>
#include <math/scmat/dist.h>

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

  RefDebugger d = keyval->describedclassvalue("debugger");
  if (d.nonnull()) {
      d->set_prefix(msg->me());
      d->set_exec(argv[0]);
    }

  // test the blocklist send and receive
  if (msg->n() > 1) {
      RefSCMatrixBlockList l = new SCMatrixBlockList;
      l->append(new SCMatrixRectBlock(0,5,0,2));
      l->append(new SCMatrixRectBlock(0,3,0,11));
      l->append(new SCMatrixLTriBlock(7,13));
      if (msg->me() == 0) {
          StateSend out(msg);
          out.target(1);
          out.copy_references();
          l.save_state(out);
          out.flush();
        }
      else if (msg->me() == 1) {
          StateRecv in(msg);
          in.source(0);
          in.copy_references();
          l.restore_state(in);
        }
    }

  RefSCMatrixKit kit = new DistSCMatrixKit;
  RefSCDimension d1(keyval->describedclassvalue("d1"));
  RefSCDimension d2(keyval->describedclassvalue("d2"));
  RefSCDimension d3(keyval->describedclassvalue("d3"));

  int nblocks = (int)sqrt(msg->n());

  // replace dimensions with dimensions that have subblocks
  d1 = new SCDimension(d1.n(), nblocks);
  d2 = new SCDimension(d2.n(), nblocks);
  d3 = new SCDimension(d3.n(), nblocks);

  matrixtest(kit,keyval,d1,d2,d3);

  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
