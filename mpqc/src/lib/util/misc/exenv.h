//
// exenv.h
//
// Copyright (C) 1997 Limit Point Systems, Inc.
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

#ifdef __GNUG__
#pragma interface
#endif

#ifndef _util_misc_exenv_h
#define _util_misc_exenv_h

/** The ExEnv class is used to find out about how
    the program is being run. */
class ExEnv {
  protected:
    static int *argc_;
    static char ***argv_;
  public:
    /// Set the argument count and vector.
    static void set_args(int &argcref, char **&argvref);
    /// Return an reference to the argument count.
    static int &argc() { return *argc_; }
    /// Return an reference to the argument vector.
    static char **&argv() { return *argv_; }
    /// Return argv[0] with the path removed.
    static const char *program_name();
};

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
