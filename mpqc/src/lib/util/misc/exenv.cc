//
// exenv.cc
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
#pragma implementation
#endif

#include <iostream.h>
#include <util/misc/exenv.h>
#include <string.h>

#ifdef HAVE_NIAMA
#include <niama.h>
#include <niama_impl.h>
#endif

int ExEnv::initialized_ = 0;
unsigned long ExEnv::mem_ = 0;
int ExEnv::nproc_ = 0;
int *ExEnv::argc_ = 0;
char ***ExEnv::argv_ = 0;

void
ExEnv::err()
{
  cout << "ExEnv: attempted to use before initialized" << endl;
  abort();
}

const char *
ExEnv::program_name()
{
  if (*argc_ == 0) return 0;
  char *start = strrchr((*argv_)[0],'/');
  if (!start) start = (*argv_)[0];
  else start++;
  return start;
}

void
ExEnv::init(int &argcref, char **&argvref)
{
  argc_ = &argcref;
  argv_ = &argvref;

  initialized_ = 1;
#ifdef HAVE_NIAMA
#if 0
  using namespace NIAMA;

  CORBA::ORB_var orb = CORBA::ORB_init(*argc_, *argv_, "mico-local-orb");
  CORBA::BOA_var boa = orb->BOA_init(*argc_, *argv_, "mico-local-boa");
  CORBA::Object_var obj = orb->bind("IDL:NIAMA/Machine:1.0");
  if (CORBA::is_nil (obj)) {
      cout << "could not bind to NIAMA server ... giving up" << endl;
      return;
    }
  Machine_var machine = Machine::_narrow (obj);
  if (CORBA::is_nil(machine)) {
      return;
    }

  nproc_ = machine->n_processor();
  mem_ = machine->memory();

  cout << "ExEnv::init: NIAMA: nproc = " << nproc_ << endl;
  cout << "ExEnv::init: NIAMA: memory = " << mem_ << endl;
#else
  using namespace NIAMA;
  // init ORB
  CORBA::ORB_var orb = CORBA::ORB_init(*argc_, *argv_, "mico-local-orb");

  // server side
  Machine_impl* machine = new Machine_impl;

  nproc_ = machine->n_processor();
  mem_ = machine->memory();

  CORBA::release(machine);
  
#endif
#endif
}

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
