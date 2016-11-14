//
// exenv.cpp
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

#include <unistd.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>

#include "mpqc/mpqc_config.h"

#ifdef HAVE_PWD_H
#include <pwd.h>
#endif
#ifdef HAVE_SYS_TYPES_H
#include <sys/types.h>
#endif

#include <cstring>

#include "mpqc/util/misc/exenv.h"
#include "mpqc/util/external/c++/memory"

using namespace std;
using namespace mpqc;

int ExEnv::initialized_ = 0;
int *ExEnv::argc_ = 0;
char ***ExEnv::argv_ = 0;
char ExEnv::hostname_[256] = { '\0' };
char ExEnv::username_[32] = { '\0' };
ostream *ExEnv::out_ = nullptr;
std::unique_ptr<ostream> ExEnv::nullstream_;

const char *
ExEnv::program_name()
{
  if (argc_ == 0 || *argc_ == 0) return 0;
  char *start = strrchr((*argv_)[0],'/');
  if (!start) start = (*argv_)[0];
  else start++;
  return start;
}

std::string
ExEnv::getenv(const std::string& name)
{
  const char *env = std::getenv(name.c_str());
  std::string value;
  if (env) {
      value = env;
      std::string::size_type offset = value.find('=');
      if (offset != std::string::npos) {
          value.erase(0,offset+1);
        }
    }
  return value;
}

std::ostream &ExEnv::out0() {
  if (!FormIO::get_debug() && FormIO::get_printnode() != FormIO::get_node()) {
    if (!nullstream_) {
      nullstream_ = std::make_unique<std::ofstream>("/dev/null");
    }
    return *nullstream_;
  }

  return outn();
}

void
ExEnv::init(int &argcref, char **&argvref)
{
  argc_ = &argcref;
  argv_ = &argvref;

#ifdef HAVE_GETHOSTNAME
  gethostname(hostname_, 256);
#else
  strcpy(hostname_, "UNKNOWN");
#endif

  memset(username_, 0, 9);
#if defined(HAVE_GETPWUID) && defined(HAVE_GETEUID)
  struct passwd *pw = getpwuid(geteuid());
  if (pw && pw->pw_name) {
      strncpy(username_, pw->pw_name, 9);
      username_[8] = 0;
    }
  else {
      strcpy(username_,"UNKNOWN");
    }
#else
  strcpy(username_,"UNKNOWN");
#endif

  initialized_ = 1;
}

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
