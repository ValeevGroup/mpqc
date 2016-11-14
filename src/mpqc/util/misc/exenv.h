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

#ifndef _util_misc_exenv_h
#define _util_misc_exenv_h

#include <cstdlib>
#include <iostream>
#include <memory>

#include <madness/world/world.h>

#include "mpqc/mpqc_config.h"
#include "mpqc/util/misc/formio.h"

namespace mpqc {

/// @addtogroup Init
/// @{

/** \brief Describes the execution environment of the program.

    \note This is a singleton.
 */
class ExEnv {
  protected:
    static int initialized_;
    static int *argc_;
    static char ***argv_;
    static char hostname_[256];
    static char username_[32];

    static std::ostream* out_;
    static std::unique_ptr<std::ostream> nullstream_;
  public:
    /// Set the argument count and vector.
    static void init(int &argcref, char **&argvref);
    /// Set the stdout stream
    static void set_out(std::ostream *o) { mpqc::FormIO::init_ostream(*o);out_=o; }
    /// Return nonzero if ExEnv has been initialized.
    static int initialized() { return argc_ != 0; }

    /// Return an reference to the argument count.
    static int &argc() { return *argc_; }
    /// Return an reference to the argument vector.
    static char **&argv() { return *argv_; }
    /// Return argv[0] with the path removed.
    static const char *program_name();
    /// Return the host name.
    static const char *hostname() { return hostname_; }
    /// Return the user name.
    static const char *username() { return username_; }

    /** Return the value of an environment variable. If it does not
        exist, then an empty string is returned. */
    static std::string getenv(const std::string& name);

    /// Return an ostream that writes from all processes.
    static std::ostream& outn();
    /// Return an ostream for error messages that writes from all processes.
    static std::ostream& errn();
    /// Return an ostream that writes from process 0 of the default World.
    /// \note FormIO::set_printnode() can be used to specify which process this
    /// prints from
    static std::ostream& out0();
    /// Return an ostream for error messages that writes from process 0 of the
    /// default World.
    /// \note FormIO::set_printnode() can be used to specify which process this
    /// prints from
    static std::ostream& err0();
    /// Return an ostream that writes from process 0 of the given World.
    static std::ostream& out0(madness::World& world);
    /// Return an ostream for error messages that writes from process 0 of the
    /// given World.
    static std::ostream& err0(madness::World& world);
};

/// @}
// end of addtogroup Init

}

#endif

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
