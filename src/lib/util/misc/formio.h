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

#include <iostream>
#include <fstream>

namespace sc {

/** This utility class is used to print only on node 0 and to
    provide attractive indentation of output. */
class SCFormIO {
  private:
    static char *default_basename_;
    static int  ready_;
    static int  xalloc_inited_;
    static long nindent_;
    static long indent_size_;
    static long skip_indent_;
    static long verbose_;
    static long initialized_;
    static int node_to_print_;
    static int debug_;
    static int parallel_;
    static int me_;
    static void init();
  public:
    static std::ios& indent(std::ios&o);
    static std::ios& decindent(std::ios&o);
    static std::ios& incindent(std::ios&o);
    static std::ios& skipnextindent(std::ios&o);

    static void setverbose(std::ios&o, long v);
    static long getverbose(std::ios&o);
    static void setindent(std::ios&o, long column);
    static long getindent(std::ios&o);
    static int  set_printnode(int);
    static int  get_printnode() { return node_to_print_; }
    static void set_debug(int);
    static int  get_debug() { return debug_; }
    static void init_mp(int me);
    static int  get_node() { return me_; }

    static void set_default_basename(const char *);
    static const char *default_basename();
    static char *fileext_to_filename(const char *extension);
    static std::string fileext_to_filename_string(const char *extension);

    static void init_ostream(std::ostream &);

    static std::ostream& license(std::ostream&);
    static std::ostream& warranty(std::ostream&);
    static std::ostream& copyright(std::ostream&);
};

std::ios& indent(std::ios&);

std::ios& decindent(std::ios&);

std::ios& incindent(std::ios&);

std::ios& skipnextindent(std::ios&);

// ///////////////////////////////////////////////////////////////////////////

class scprintf;
std::ostream& operator<<(std::ostream&, const scprintf&);

/** This class allows <tt>printf</tt>-like output to be sent
    to an <tt>ostream</tt>. */
class scprintf {
  private:
    char str_[1024];

  public:
    scprintf(const char*,...);
    const char* str() const { return str_; }
};

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
