//
// file.cc
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

#include <sstream>
#include <util/misc/string.h>

#ifdef HAVE_CONFIG_H
#include <mpqc_config.h>
#endif


#include <util/misc/formio.h>
#include <util/misc/scexception.h>
#include <util/group/file.h>

#include <util/group/fileproc.h>
#if defined(HAVE_MPI)
#  include <util/group/messmpi.h>
//#  include <util/group/filemtmpi.h>
#endif

using namespace std;
using namespace sc;

//////////////////////////////////////////////////////////////////////
// FileGrp members

static ClassDesc FileGrp_cd(
  typeid(FileGrp),"FileGrp",1,"public DescribedClass",
  0, 0, 0);

FileGrp::FileGrp()
{
  debug_ = 0;

  datafile_ = 0;
  filename_ = 0;
  offsets_ = 0;

  init_locks();
}

FileGrp::FileGrp(const Ref<KeyVal>& keyval)
{
  debug_ = keyval->intvalue("debug");

  datafile_ = 0;
  filename_ = 0;
  offsets_ = 0;

  init_locks();
}

FileGrp::~FileGrp()
{
  delete[] offsets_;
  delete[] locks_;
  datafile_ = 0;
  filename_ = 0;
}

void
FileGrp::open()
{
  if (filename_) {
    datafile_ = open(filename_, O_RDWR);
  }
  else {
    throw ProgrammingError("open() called but filename has not been set",
                           __FILE__, __LINE__, class_desc());
  }
}

void
FileGrp::close()
{
  close(datafile_);
}

void
FileGrp::set_filename(char *name)
{
  if (filename_) {
    throw ProgrammingError("set_filename() called but filename has been set already",
                           __FILE__, __LINE__, class_desc());
  }
  else
    filename_ = strdup(name);
}

const char*
FileGrp::get_filename() const
{
  return filename_;
}

void
FileGrp::init_locks()
{
  Ref<ThreadGrp> thgrp = ThreadGrp::get_default_threadgrp();
  nlock_ = 2 * thgrp->nthread();
  locks_ = new Ref<ThreadLock>[nlock_];
  for (int i=0; i<nlock_; i++) locks_[i] = thgrp->new_lock();
}

FileGrp *
FileGrp::initial_filegrp()
{
  int argc = 0;
  return initial_filegrp(argc,0);
}

FileGrp *
FileGrp::initial_filegrp(int &argc, char *argv[])
{
  FileGrp *grp = 0;

  char *keyval_string = 0;

  // see if a file group is given on the command line
  if (argc && argv) {
      for (int i=0; i<argc; i++) {
	  if (argv[i] && !strcmp(argv[i], "-filegrp")) {
              char *filegrp_string = argv[i];
              i++;
              if (i >= argc) {
                  ExEnv::errn() << "-filegrp must be following by an argument"
                       << endl;
                  abort();
                }
              keyval_string = argv[i];
              // move the filegrp arguments to the end of argv
              int j;
              for (j=i+1; j<argc; j++) {
                  argv[j-2] = argv[j];
                }
              argv[j++] = filegrp_string;
              argv[j++] = keyval_string;
              // decrement argc to hide the last two arguments
              argc -= 2;
              break;
            }
        }
    }

  if (!keyval_string) {
      // find out if the environment gives the containing file group
      keyval_string = getenv("FILEGRP");
      if (keyval_string) {
          if (!strncmp("FILEGRP=", keyval_string, 11)) {
              keyval_string = strchr(keyval_string, '=');
            }
          if (*keyval_string == '=') keyval_string++;
        }
    }

  // if keyval input for a file group was found, then
  // create it.
  if (keyval_string) {
      //ExEnv::outn() << "Creating FileGrp from \"" << keyval_string << "\"" << endl;
      Ref<ParsedKeyVal> strkv = new ParsedKeyVal();
      strkv->parse_string(keyval_string);
      Ref<DescribedClass> dc = strkv->describedclassvalue();
      grp = dynamic_cast<FileGrp*>(dc.pointer());
      if (dc.null()) {
          ExEnv::errn() << "initial_filegrp: couldn't find a FileGrp in "
               << keyval_string << endl;
          abort();
        }
      else if (!grp) {
          ExEnv::errn() << "initial_filegrp: wanted FileGrp but got "
               << dc->class_name() << endl;
          abort();
        }
      // prevent an accidental delete
      grp->reference();
      strkv = 0;
      dc = 0;
      // accidental delete not a problem anymore since all smart pointers
      // to grp are dead
      grp->dereference();
      return grp;
    }

  return grp;
}

void
FileGrp::activate()
{
}

void
FileGrp::deactivate()
{
}

void
FileGrp::print(ostream&o) const
{
  o << scprintf("FileGrp (node %d):\n", me());
  o << scprintf("%d: n = %d\n", me(), n());
  for (int i=0; i<=n_; i++) {
      o << scprintf("%d: offset[%d] = %5d\n", me(), i, offsets_[i]);
    }
}

void
FileGrp::sum_reduction(double *data, distsize_t doffset, int dlength)
{
  distsize_t offset = doffset * sizeof(double);
  int length = dlength * sizeof(double);

  if (offset + length > totalsize()) {
      ExEnv::errn() << scprintf("FileGrp::sum_reduction: arg out of range\n");
      abort();
    }

  double *source_data = (double*) obtain_readwrite(offset, length);

  for (int i=0; i<dlength; i++) {
      source_data[i] += data[i];
    }

  release_readwrite((void*) source_data, offset, length);
}

void
FileGrp::sum_reduction_on_node(double *data, size_t doffset, int dlength,
                                 int node)
{
  if (node == -1) node = me();

  sum_reduction(data, doffset + offset(node)/sizeof(double),
                dlength);
}

void
FileGrp::catchup()
{
  return;
}

void
FileGrp::obtain_local_lock(size_t start, size_t fence)
{
  distsize_t locked_region_size = 1 + localsize()/nlock_;
  int lstart = start/locked_region_size;
  int llast = fence/locked_region_size;
  for (int i=lstart; i<=llast; i++) {
      locks_[i]->lock();
    }
}

void
FileGrp::release_local_lock(size_t start, size_t fence)
{
  distsize_t locked_region_size = 1 + localsize()/nlock_;
  int lstart = start/locked_region_size;
  int llast = fence/locked_region_size;
  for (int i=lstart; i<=llast; i++) {
      locks_[i]->unlock();
    }
}

static Ref<FileGrp> default_filegrp;

void
FileGrp::set_default_filegrp(const Ref<FileGrp>& grp)
{
  default_filegrp = grp;
}

FileGrp*
FileGrp::get_default_filegrp()
{
  if (default_filegrp) return default_filegrp.pointer();

  Ref<MessageGrp> msg = MessageGrp::get_default_messagegrp();

#if defined(HAVE_MPI) && defined(DEFAULT_MTMPI)
  Ref<ThreadGrp> thr = ThreadGrp::get_default_threadgrp();
//  default_filegrp = new MTMPIFileGrp(msg,thr);
  return default_filegrp.pointer();
#endif


  if (msg.null()) {
      ExEnv::errn() << scprintf("FileGrp::get_default_filegrp: requires default MessageGrp if default behavior not configured\n");
      abort();
    }
#if defined(HAVE_MPI)
  else if (msg->class_desc() == ::class_desc<MPIMessageGrp>()) {
      Ref<ThreadGrp> thr = ThreadGrp::get_default_threadgrp();
//      default_filegrp = new MTMPIFileGrp(msg,thr);
      return default_filegrp.pointer();
    }
#endif
  else if (msg->n() == 1) {
      default_filegrp = new ProcFileGrp();
      return default_filegrp.pointer();
    }
  else {
      ExEnv::errn() << scprintf("FileGrp::get_default_filegrp: cannot create "
              "default for \"%s\"\n.", msg->class_name());
      abort();
    }

  if (default_filegrp.null()) {
      ExEnv::err0() << scprintf("WARNING: FileGrp::get_default_filegrp(): failed\n");
      default_filegrp = new ProcFileGrp;
    }
  return default_filegrp.pointer();
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
