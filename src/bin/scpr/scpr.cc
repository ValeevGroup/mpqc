//
// scpr.cc
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

#include <unistd.h>
#include <new>

#include <util/misc/exenv.h>
#include <util/misc/formio.h>
#include <util/group/mstate.h>
#include <util/group/message.h>
#include <util/group/thread.h>
#include <util/group/memory.h>

// Force linkages:
#include <chemistry/qc/wfn/linkage.h>
#include <chemistry/qc/dft/linkage.h>
#include <chemistry/qc/mbpt/linkage.h>
#ifdef HAVE_LIBINT2
#  include <chemistry/qc/libint2/linkage.h>
#  include <chemistry/qc/mbptr12/linkage.h>
#  include <chemistry/qc/ccr12/linkage.h>
#endif
#ifdef HAVE_PSI3
#  include <chemistry/qc/psi/linkage.h>
#endif
#include <util/state/linkage.h>

using namespace std;
using namespace sc;

#ifdef HAVE_MPI
#define OMPI_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#include <mpi.h>
#endif

//////////////////////////////////////////////////////////////////////////

static void
clean_up(void)
{
  MessageGrp::set_default_messagegrp(0);
}

static void
out_of_memory()
{
  ExEnv::errn() << "ERROR: out of memory" << endl;
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
  // want an MPIMessageGrp.  The command name is used to let scpr
  // know that an early init is needed.
  if (!strcmp(ExEnv::program_name(), "scpr-mpi")) {
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
  char **objects = 0;
  int nfile = 0;
  int nobject = 0;
  for (i=1; i<argc; i++) {
      char *arg = argv[i];
      if (!strcmp(arg,"-messagegrp")) i++;
      else if (!strcmp(arg,"-memorygrp")) i++;
      else if (!strcmp(arg,"-threadgrp")) i++;
      else if (!strcmp(arg,"-W")) working_dir = argv[++i];
      else if (!strcmp(arg,"-d")) debug = 1;
      else if (!strcmp(arg,"-h")) help = 1;
      else if (!strcmp(arg,"-l")) SCFormIO::setverbose(ExEnv::outn(),1);
      else if (!strcmp(arg,"-v")) version = 1;
      else if (!strcmp(arg,"-w")) warranty = 1;
      else if (!strcmp(arg,"-L")) license = 1;
      else if (!strcmp(arg,"-o")) {
          if (argc > i+1) {
              char **newobjects = new char *[nobject+1];
              memcpy(newobjects, objects, sizeof(char*)*nobject);
              delete[] objects;
              objects = newobjects;
              objects[nobject++] = argv[++i];
            }
          else help = 1;
        }
      else {
          char **newfiles = new char *[nfile+1];
          memcpy(newfiles, files, sizeof(char*)*nfile);
          delete[] files;
          files = newfiles;
          files[nfile++] = arg;
        }
    }

  if (help || nobject == 0 || nfile == 0) {
      ExEnv::out0()
           << indent << "scpr version " << MPQC_VERSION << endl
           << SCFormIO::copyright << endl
           << indent << "usage: " << argv[0] << " [options] file ..." << endl
           << indent << "where options are chosen from:" << endl
           << indent << "-o <$val> (print the object with the name $val)"<<endl
           << indent << "-memorygrp <$val> (which memory group to use)" << endl
           << indent << "-threadgrp <$val> (which thread group to use)" << endl
           << indent << "-messagegrp <$val> (which message group to use)"<<endl
           << indent << "-W <$val> (set the working directory)" << endl
           << indent << "-d (turn on debugging)" << endl
           << indent << "-l (verbose printing)" << endl
           << indent << "-v (print the version)" << endl
           << indent << "-w (print the warranty)" << endl
           << indent << "-L (print the license)" << endl
           << indent << "-h (print this help)" << endl;

      ExEnv::out0() << endl
           << indent << "object names take the form classname:ordinal_number"
           << endl
           << indent << "at least one file and object name must be given"
           << endl;
      exit(0);
    }

  if (version) {
      ExEnv::out0()
         << indent << "scpr version " << MPQC_VERSION << endl
         << SCFormIO::copyright;
    exit(0);
  }

  if (warranty) {
      ExEnv::out0()
         << indent << "scpr version " << MPQC_VERSION << endl
         << SCFormIO::copyright << endl
         << SCFormIO::warranty;
    exit(0);
  }

  if (license) {
      ExEnv::out0()
         << indent << "scpr version " << MPQC_VERSION << endl
         << SCFormIO::copyright << endl
         << SCFormIO::license;
    exit(0);
  }

  // set the working dir
  if (working_dir && strcmp(working_dir,"."))
    chdir(working_dir);

  // get the message group.  first try the commandline and environment
  Ref<MessageGrp> grp = MessageGrp::initial_messagegrp(argc, argv);
  if (grp)
    MessageGrp::set_default_messagegrp(grp);
  else
    grp = MessageGrp::get_default_messagegrp();

  // get the thread group.  first try the commandline and environment
  Ref<ThreadGrp> thread = ThreadGrp::initial_threadgrp(argc, argv);
  if (thread)
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
      if (nfile > 1) {
          ExEnv::out0() << indent << files[i] << ":" << endl;
          ExEnv::out0() << incindent;
        }
      BcastStateInBin s(grp,files[i]);
      for (int j=0; j<nobject; j++) {
          if (nobject > 1) {
              ExEnv::out0() << indent << objects[j] << ":" << endl;
              ExEnv::out0() << incindent;
            }
          Ref<SavableState> o;
          o << SavableState::dir_restore_state(s,objects[j]);
          o->print(ExEnv::out0());
          if (nobject > 1) ExEnv::out0() << decindent;
        }
      if (nfile > 1) ExEnv::out0() << decindent;
    }

  delete[] files;

  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
