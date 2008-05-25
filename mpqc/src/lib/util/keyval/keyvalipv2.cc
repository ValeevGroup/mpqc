//
// keyvalipv2.cc
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

#include <scconfig.h>
#include <iostream>
#ifdef HAVE_SSTREAM
#  include <sstream>
#else
#  include <strstream.h>
#endif
#include <fstream>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>

#include <util/misc/string.h>
#include <util/misc/formio.h>
#include <util/keyval/ipv2.h>
#include <util/keyval/keyval.h>

using namespace std;
using namespace sc;

ParsedKeyVal::ParsedKeyVal(IPV2*i):
nfile(0),
file(0),
nfp(0)
{
  ipv2 = i;
}

ParsedKeyVal::ParsedKeyVal():
nfile(0),
file(0),
nfp(0)
{
  ipv2 = new IPV2;
}

ParsedKeyVal::ParsedKeyVal(const char* name):
nfile(0),
file(0),
nfp(0)
{
  ipv2 = new IPV2;
  read(name);
}

ParsedKeyVal::ParsedKeyVal(istream& fp):
nfile(0),
file(0),
nfp(0)
{
  ipv2 = new IPV2;
  read(fp);
}

ParsedKeyVal::ParsedKeyVal(const char* keyprefix, const Ref<KeyVal>& keyval):
nfile(0),
file(0),
nfp(0)
{
  ipv2 = new IPV2;

  char* filespec = new char[strlen(keyprefix)+6];
  strcpy(filespec,keyprefix);
  strcat(filespec,"files");

  char* dirspec = new char[strlen(keyprefix)+6];
  strcpy(dirspec,keyprefix);
  strcat(dirspec,"dir");

  std::string directory = keyval->stringvalue(dirspec);
  if (directory.empty()) {
      directory = ExEnv::getenv_string("SCLIBDIR");
      if (directory.empty()) {
          struct stat sb;
          const char *dir = INSTALLED_SCLIBDIR;
#ifdef SRC_SCLIBDIR
          if (stat(dir, &sb) != 0) {
              ExEnv::out0() << indent << "WARNING: could not find "
                   << dir << endl;
              dir = SRC_SCLIBDIR;
            }
#endif
          directory = dir;
        }
    }

  int nfiles = keyval->count(filespec);
  for (int i=0; i<nfiles; i++) {
      std::string filename = keyval->stringvalue(filespec,i);
      std::string fullname;
      if (!directory.empty()) {
          fullname = directory + filename;
        }
      else {
          fullname = filename;
        }
      read(fullname.c_str());
    }

  delete[] dirspec;
  delete[] filespec;

}

void
ParsedKeyVal::cat_files(const char* keyprefix, const Ref<KeyVal>& keyval,
                        ostream &ostr)
{
  std::string filespec = keyprefix;
  filespec += "files";

  std::string dirspec = keyprefix;
  dirspec += "dir";

  std::string directory = keyval->stringvalue(dirspec.c_str());
  if (directory.empty()) {
      directory = ExEnv::getenv_string("SCLIBDIR");
      if (directory.empty()) {
          struct stat sb;
          const char *dir = INSTALLED_SCLIBDIR;
#ifdef SRC_SCLIBDIR
          if (stat(dir, &sb) != 0) {
              ExEnv::out0() << indent << "WARNING: could not find "
                   << dir << endl;
              dir = SRC_SCLIBDIR;
            }
#endif
          directory = dir;
        }
    }

  int nfiles = keyval->count(filespec.c_str());
  for (int i=0; i<nfiles; i++) {
      std::string filename = keyval->stringvalue(filespec.c_str(),i);
      std::string fullname;
      if (!directory.empty()) {
          fullname = directory + filename;
        }
      else {
          fullname = filename;
        }
      ifstream is(fullname.c_str());
      is >> ostr.rdbuf();
    }
}

void
ParsedKeyVal::read(const std::string &name)
{
  std::string tmp(name);
  read(tmp.c_str());
}

void
ParsedKeyVal::read(const char* name)
{
  ifstream infp(name,ios::in);
  if (infp.bad()) {
    ExEnv::errn() << "ParsedKeyVal couldn't open " << name << endl;
    exit(1);
    }

  int i;
  char**newfile = new char*[nfile+1];
  for (i=0; i<nfile; i++) newfile[i] = file[i];
  if (file) delete[] file;
  file = newfile;
  newfile[nfile] = strdup(name);
  nfile++;

  read(infp);
  nfp--; // read(infp) will incr nfp, but this isn't desired so undo
}

void ParsedKeyVal::read(istream&infp)
{
  nfp++;
  ipv2->read(infp,ExEnv::errn(),"<stream>");
}

void
ParsedKeyVal::parse_string(const char *str)
{
#ifdef HAVE_SSTREAM
  istringstream in(str);
#else
  istrstream in(str);
#endif
  ipv2->read(in,ExEnv::errn(),"<string>");
}

ParsedKeyVal::~ParsedKeyVal()
{
  delete ipv2;
  for (int i=0; i<nfile; i++) free(file[i]);
  delete[] file;
}

static KeyVal::KeyValError maperr(IPV2::Status err)
  {
  if (err == IPV2::OK            ) return KeyVal::OK;
  if (err == IPV2::KeyNotFound ) return KeyVal::UnknownKeyword;
  if (err == IPV2::OutOfBounds ) return KeyVal::UnknownKeyword;
  if (err == IPV2::Malloc        ) return KeyVal::OperationFailed;
  if (err == IPV2::NotAnArray  ) return KeyVal::UnknownKeyword;
  if (err == IPV2::NotAScalar  ) return KeyVal::HasNoValue;
  if (err == IPV2::Type          ) return KeyVal::WrongType;
  if (err == IPV2::HasNoValue  ) return KeyVal::HasNoValue;
  if (err == IPV2::ValNotExpd  ) return KeyVal::OperationFailed;
  return KeyVal::OperationFailed;
  }

std::string ParsedKeyVal::stringrep(const char* key)
{
  const char* result;
  seterror(maperr(ipv2->value_v((char *)key,&result,0,0)));
  std::string retval;
  if (error() == OK) {
      retval = result;
    }
  return retval;
}

const char*
ParsedKeyVal::classname(const char* key)
{
  const char* result;
  seterror(maperr(ipv2->classname_v((char *)key,&result,0,0)));
  return result;
}

const char*
ParsedKeyVal::truekeyword(const char*key)
{
  const char* result;
  seterror(maperr(ipv2->truekeyword_v((char *)key,&result,0,0)));
  if (!result && error() == OK) return key;
  else return result;
}

void ParsedKeyVal::errortrace(ostream&fp)
{
  fp << indent << "ParsedKeyVal: error: \"" << errormsg() << "\"" << endl;
  if (nfp) {
      fp << indent
         << "    reading from " << nfp << " files with unknown names" << endl;
    }
  for (int i=0; i<nfile; i++) {
      fp << indent << "    reading from \"" << file[i] << "\"" << endl;
    }
}

void ParsedKeyVal::dump(ostream&fp)
{
  fp << indent << "ParsedKeyVal: error: \"" << errormsg() << "\"" << endl;
  if (nfp) {
      fp << indent
         << "    reading from " << nfp << " files with unknown names" << endl;
    }
  for (int i=0; i<nfile; i++) {
      fp << indent << "    reading from \"" << file[i] << "\"" << endl;
    }
  fp << indent << "The IPV2 tree:" << endl;
  ipv2->print_tree(fp);

}

void ParsedKeyVal::print_unseen(ostream&fp)
{
  ipv2->print_unseen(fp);
}

int ParsedKeyVal::have_unseen()
{
  return ipv2->have_unseen();
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
