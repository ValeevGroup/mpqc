//
// rentest.cc
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

#include <iostream>

#include <util/keyval/keyval.h>
#include <util/render/oogl.h>
#include <util/render/object.h>
#include <util/render/sphere.h>

#include <util/render/linkage.h>

using namespace std;
using namespace sc;

int
main(int argc, char* argv[])
{
  Ref<KeyVal> keyval = new ParsedKeyVal(SRCDIR "/rentest.in");
  cout << "getting render" << endl << flush;
  Ref<Render> render; render << keyval->describedclassvalue("render");
  cout << "getting object" << endl << flush;
  Ref<RenderedObject> object; object << keyval->describedclassvalue("object");

  cout << "rendering object" << endl << flush;
  render->render(object);
  cout << "rendered object" << endl << flush;

  cout << "getting rid of keyval" << endl << flush;
  keyval = 0;
  cout << "getting rid of render" << endl << flush;
  render = 0;
  cout << "getting rid of object" << endl << flush;
  object = 0;

  cout << "main done" << endl << flush;

  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
