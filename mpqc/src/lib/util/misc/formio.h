//
// formio.h
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

#ifndef _util_misc_formio_h
#define _util_misc_formio_h

#include <iostream.h>
#include <fstream.h>

class SCFormIO {
  private:
    static char *default_basename_;
    static int  ready_;
    static long nindent_;
    static long indent_size_;
    static long skip_indent_;
    static long verbose_;
    static int node_to_print_;
    static int debug_;
    static int parallel_;
    static int me_;
    static ofstream nullstream_;
    static void init();
  public:
    static ios& indent(ios&o);
    static ios& decindent(ios&o);
    static ios& incindent(ios&o);
    static ios& skipnextindent(ios&o);
    static ostream& node0(ostream&o);

    static void setverbose(ios&o, long v);
    static long getverbose(ios&o);
    static void setindent(ios&o, long column);
    static long getindent(ios&o);
    static void set_printnode(int);
    static void set_debug(int);
    static void init_mp(int me);

    static void set_default_basename(const char *);
    static const char *default_basename();
    static char *fileext_to_filename(const char *extension);

    static ostream& license(ostream&);
    static ostream& warranty(ostream&);
    static ostream& copyright(ostream&);
};

ios& indent(ios&);

ios& decindent(ios&);

ios& incindent(ios&);

ios& skipnextindent(ios&);

ostream& node0(ostream&);

/////////////////////////////////////////////////////////////////////////////

class scprintf {
  private:
    char str[1024];

  public:
    scprintf(const char*,...);
    friend ostream& operator<<(ostream&, const scprintf&);
};

ostream& operator<<(ostream&, const scprintf&);

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
