//
// blockedtest.cc
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

#include <util/keyval/keyval.h>
#include <math/scmat/local.h>
#include <math/scmat/blocked.h>

using namespace sc;

void matrixtest(Ref<SCMatrixKit>, Ref<KeyVal>,
                RefSCDimension d1,RefSCDimension d2,RefSCDimension d3,
                bool have_svd);

int
main(int argc, char **argv)
{
  int i;
  int nblks;

  Ref<KeyVal> keyval = new ParsedKeyVal(SRCDIR "/matrixtest.in");

  Ref<MessageGrp> msg = MessageGrp::initial_messagegrp(argc, argv);

  if (msg.null()) {
      msg << keyval->describedclassvalue("messagegrp");

      if (msg.null()) {
          std::cerr << indent << "Couldn't initialize MessageGrp\n";
          abort();
        }
    }

  MessageGrp::set_default_messagegrp(msg);

  Ref<SCMatrixKit> subkit = new LocalSCMatrixKit;
  Ref<BlockedSCMatrixKit> kit = new BlockedSCMatrixKit(subkit);

  nblks = keyval->intvalue("nblocks");
  if (!nblks)
    nblks=3;

  RefSCDimension sd1; sd1 << keyval->describedclassvalue("d1");
  RefSCDimension sd2; sd2 << keyval->describedclassvalue("d2");
  RefSCDimension sd3; sd3 << keyval->describedclassvalue("d2");
  
  RefSCDimension d1(new SCDimension(sd1.n()*nblks,nblks));
  RefSCDimension d2(new SCDimension(sd2.n()*nblks,nblks));
  RefSCDimension d3(new SCDimension(sd3.n()*nblks,nblks));

  for (i=0; i < nblks; i++) {
    d1->blocks()->set_subdim(i,sd1);
    d2->blocks()->set_subdim(i,sd2);
    d3->blocks()->set_subdim(i,sd3);
  }

  matrixtest(kit,keyval,d1,d2,d3,false);
  
  d1 = new SCDimension(sd1.n(),1);
  d2 = new SCDimension(sd2.n(),1);
  d3 = new SCDimension(sd3.n(),1);

  d1->blocks()->set_subdim(0,sd1);
  d2->blocks()->set_subdim(0,sd2);
  d3->blocks()->set_subdim(0,sd3);

  matrixtest(kit,keyval,d1,d2,d3,false);
  d1=0;

  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
