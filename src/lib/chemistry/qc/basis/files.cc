//
// files.cc
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

#include <ctype.h>
#include <string.h>
#include <stdio.h>
#include <fstream>

#include <mpqc_config.h>
#include <sstream>

#include <sys/stat.h>

#include <util/misc/formio.h>
#include <util/group/message.h>
#include <chemistry/qc/basis/files.h>

using namespace std;
using namespace sc;

BasisFileSet::BasisFileSet(const Ref<KeyVal>& keyval)
{
  nbasissets_ = 0;
  basissets_ = 0;

  dir_[0] = keyval->stringvalue("basisdir");
  dir_[1] = ExEnv::getenv_string("MPQC_DATA_PATH");
  if (dir_[1].empty()) {
      struct stat sb;
      const char *bdir = MPQCDATAPATH "/basis";
#ifdef SRC_MPQC_DATA_PATH
      if (stat(bdir, &sb) != 0) {
          bdir = SRC_MPQC_DATA_PATH "/basis";
        }
#endif
      dir_[1] = bdir;
    }
  else {
      dir_[1] += "/basis";
    }
}

BasisFileSet::~BasisFileSet()
{
  int i;
  for (i=0; i<nbasissets_; i++) delete[] basissets_[i];
  delete[] basissets_;
}

Ref<KeyVal>
BasisFileSet::keyval(const Ref<KeyVal> &keyval, const char *basisname)
{
  int i;

  // check if the basisname has already been located
  for (i=0; i<nbasissets_; i++) {
      if (!strcmp(basisname, basissets_[i])) return keyval;
    }

  char *filename = new char[strlen(basisname)+1];

  for (i=0; basisname[i] != '\0'; i++) {
      if (basisname[i] >= 'A' && basisname[i] <= 'Z') {
          filename[i] = tolower(basisname[i]);
        }
      else if (basisname[i] == ','
               || basisname[i] == ' ') {
          filename[i] = '_';
        }
      else if (basisname[i] == '+') {
          filename[i] = 'P';
        }
      else if (basisname[i] == '*') {
          filename[i] = 'S';
        }
      else if (basisname[i] == '(') {
          filename[i] = 'L';
        }
      else if (basisname[i] == ')') {
          filename[i] = 'R';
        }
      else if (basisname[i] == '/') {
          filename[i] = 'X';
        }
      else {
          filename[i] = basisname[i];
        }
    }
  filename[i] = '\0';

  Ref<MessageGrp> grp = MessageGrp::get_default_messagegrp();

  // find the basis file
  Ref<KeyVal> newkeyval(keyval);
  for (i=0; i<2; i++) {
      if (dir_[i].empty()) continue;
      if (grp->me() == 0) {
          char *path = new char[dir_[i].size() + strlen(filename) + 5];
          strcpy(path, dir_[i].c_str());
          strcat(path, "/");
          strcat(path, filename);
          strcat(path, ".kv");

          // test to see if the file can be opened read only.
          ifstream is(path);
          if (is.good()) {
              int status = 1;
              ExEnv::out0() << indent << "Reading file " << path << "." << endl;
              grp->bcast(status);
              ostringstream ostrs;
              is >> ostrs.rdbuf();
              int n = 1 + strlen(ostrs.str().c_str());
              char *in_char_array = strcpy(new char[n],ostrs.str().c_str());
              grp->bcast(n);
              grp->bcast(in_char_array, n);
              Ref<ParsedKeyVal> parsedkv = new ParsedKeyVal;
              parsedkv->parse_string(in_char_array);
              delete[] in_char_array;
              Ref<KeyVal> libkeyval = parsedkv.pointer();
              newkeyval = new AggregateKeyVal(keyval,libkeyval);
              delete[] path;
              break;
            }
          else {
              int status = 0;
              grp->bcast(status);
            }
          delete[] path;
        }
      else {
          int status;
          grp->bcast(status);
          if (status) {
              int n;
              grp->bcast(n);
              char *in_char_array = new char[n];
              grp->bcast(in_char_array, n);
              Ref<ParsedKeyVal> parsedkv = new ParsedKeyVal;
              parsedkv->parse_string(in_char_array);
              delete[] in_char_array;
              Ref<KeyVal> libkeyval = parsedkv.pointer();
              newkeyval = new AggregateKeyVal(keyval,libkeyval);
              break;
            }
        }
    }

  // add the current basis set to basissets_
  char **newbasissets = new char*[nbasissets_+1];
  for (i=0; i<nbasissets_; i++) newbasissets[i] = basissets_[i];
  newbasissets[nbasissets_] = strcpy(new char[strlen(basisname)+1],
                                     basisname);
  nbasissets_++;
  delete[] basissets_;
  basissets_ = newbasissets;

  delete[] filename;

  return newkeyval;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
