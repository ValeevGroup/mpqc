//
// scls.cc
//
// Copyright (C) 1997 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@limitpt.com>
// Maintainer: LPS
//
// This file is part of the SC Toolkit utilities.
//
// The SC Toolkit utilities are free software; you can redistribute them
// and/or modify them under the terms of the GNU General Public License as
// published by the Free Software Foundation; either version 2, or (at your
// option) any later version.
//
// The SC Toolkit utilities are distributed in the hope that they will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with the SC Toolkit utilities; see the file COPYING.  If not, write to
// the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
//
// The U.S. Government is granted a limited license as per AL 91-7.
//

#ifdef HAVE_CONFIG_H
#include <mpqc_config.h>
#endif

#include <new>
#include <unistd.h>

#include <util/misc/scexception.h>
#include <util/misc/exenv.h>
#include <util/misc/formio.h>
#include <util/group/mstate.h>
#include <util/group/message.h>
#include <util/group/thread.h>
#include <util/group/memory.h>

// Force linkages:
#include <util/state/linkage.h>

#ifdef HAVE_MPI
#define MPICH_SKIP_MPICXX
#include <mpi.h>
#endif

using namespace std;
using namespace sc;

//////////////////////////////////////////////////////////////////////////

// This forces the exception classes to be linked in.  Otherwise they won't
// for single pass linkage of static libraries.  This is because of library
// ordering and interdependency issues.
static SCException ex;

static void
clean_up(void)
{
  MessageGrp::set_default_messagegrp(0);
}

static void
out_of_memory()
{
  cerr << "ERROR: out of memory" << endl;
  abort();
}

int
main(int argc, char *argv[])
{
  atexit(clean_up);
  std::set_new_handler(out_of_memory);

  ExEnv::init(argc, argv);
  ExEnv::set_out(&cout);

#ifdef HAVE_MPI
  // MPI is allowed wait until MPI_Init to fill in argc and argv,
  // so we may have to call MPI_Init before we even know that we
  // want an MPIMessageGrp.  The command name is used to let scls
  // know that an early init is needed.
  if (!strcmp(ExEnv::program_name(), "scls-mpi")) {
    MPI_Init(&argc, &argv);
  }
#endif

  int i;
  int debug = 0;
  int version = 0;
  int warranty = 0;
  int license = 0;
  int help = 0;
  const char *working_dir = 0;
  char **files = 0;
  int nfile = 0;
  for (i=1; i<argc; i++) {
      char *arg = argv[i];
      if (!strcmp(arg,"-messagegrp")) i++;
      else if (!strcmp(arg,"-memorygrp")) i++;
      else if (!strcmp(arg,"-threadgrp")) i++;
      else if (!strcmp(arg,"-W")) working_dir = argv[++i];
      else if (!strcmp(arg,"-d")) debug = 1;
      else if (!strcmp(arg,"-h")) help = 1;
      else if (!strcmp(arg,"-v")) version = 1;
      else if (!strcmp(arg,"-l")) SCFormIO::setverbose(cout,1);
      else if (!strcmp(arg,"-w")) warranty = 1;
      else if (!strcmp(arg,"-L")) license = 1;
      else {
          char **newfiles = new char *[nfile+1];
          memcpy(newfiles, files, sizeof(char*)*nfile);
          delete[] files;
          files = newfiles;
          files[nfile++] = arg;
        }
    }

  if (help) {
      ExEnv::out0()
           << indent << "scls version " << MPQC_VERSION << endl
           << SCFormIO::copyright << endl
           << indent << "-memorygrp <$val> (which memory group to use)" << endl
           << indent << "-threadgrp <$val> (which thread group to use)" << endl
           << indent << "-messagegrp <$val> (which message group to use)"<<endl
           << indent << "-W <$val> (set the working directory)" << endl
           << indent << "-d (turn on debugging)" << endl
           << indent << "-v (print the version)" << endl
           << indent << "-w (print the warranty)" << endl
           << indent << "-l (detailed list of objects)" << endl
           << indent << "-L (print the license)" << endl
           << indent << "-h (print this help)" << endl;
      exit(0);
    }
  
  if (version) {
      ExEnv::out0()
         << indent << "scls version " << MPQC_VERSION << endl
         << SCFormIO::copyright;
    exit(0);
  }
  
  if (warranty) {
      ExEnv::out0()
         << indent << "scls version " << MPQC_VERSION << endl
         << SCFormIO::copyright << endl
         << SCFormIO::warranty;
    exit(0);
  }
  
  if (license) {
      ExEnv::out0()
         << indent << "scls version " << MPQC_VERSION << endl
         << SCFormIO::copyright << endl
         << SCFormIO::license;
    exit(0);
  }

  // set the working dir
  if (working_dir && strcmp(working_dir,"."))
    chdir(working_dir);

  // get the message group.  first try the commandline and environment
  Ref<MessageGrp> grp = MessageGrp::initial_messagegrp(argc, argv);
  if (grp.nonnull())
    MessageGrp::set_default_messagegrp(grp);
  else
    grp = MessageGrp::get_default_messagegrp();

  // get the thread group.  first try the commandline and environment
  Ref<ThreadGrp> thread = ThreadGrp::initial_threadgrp(argc, argv);
  if (thread.nonnull())
    ThreadGrp::set_default_threadgrp(thread);
  else
    thread = ThreadGrp::get_default_threadgrp();

  // set up output classes
  SCFormIO::setindent(ExEnv::outn(), 0);
  SCFormIO::setindent(ExEnv::errn(), 0);
  SCFormIO::setindent(cout, 0);
  SCFormIO::setindent(cerr, 0);

  SCFormIO::set_printnode(0);
  if (grp->n() > 1)
    SCFormIO::init_mp(grp->me());

  if (debug)
    SCFormIO::set_debug(1);


  for (i=0; i<nfile; i++) {
      ExEnv::out0() << indent << files[i] << ":" << endl;
      BcastStateInBin s(grp,files[i]);
      ExEnv::out0() << incindent;
      s.list_objects(ExEnv::out0());
      if (s.has_directory() && !s.seekable()) {
          ExEnv::out0() << "(objects cannot be listed since cannot seek file)" << endl;
        }
      ExEnv::out0() << decindent;
    }

  delete[] files;

  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
