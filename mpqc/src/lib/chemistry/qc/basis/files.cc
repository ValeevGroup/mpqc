//
// files.cc
//
// Copyright (C) 1996 Limit Point Systems, Inc.
//
// Author: Curtis Janssen <cljanss@ca.sandia.gov>
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
#include <chemistry/qc/basis/files.h>
#include <stdio.h>

BasisFileSet::BasisFileSet(const RefKeyVal& keyval)
{
  nbasissets_ = 0;
  basissets_ = 0;

  dir_[0] = keyval->pcharvalue("basisdir");
  dir_[1] = getenv("SCLIBDIR");
  if (dir_[1]) {
      char *tmp = strchr(dir_[1],'=');
      if (!tmp) tmp = dir_[1];
      else tmp = &tmp[1];
      dir_[1] = strcpy(new char[strlen(tmp)+6+1], tmp);
      strcat(dir_[1], "/basis");
    }
  else {
      dir_[1] = strcpy(new char[strlen(BASISDIR)+1], BASISDIR);
    }
}

BasisFileSet::~BasisFileSet()
{
  int i;
  for (i=0; i<nbasissets_; i++) delete[] basissets_[i];
  delete[] basissets_;
  delete[] dir_[0];
  delete[] dir_[1];
}

RefKeyVal
BasisFileSet::keyval(const RefKeyVal &keyval, const char *basisname)
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
      else {
          filename[i] = basisname[i];
        }
    }
  filename[i] = '\0';

  // find the basis file
  RefKeyVal newkeyval(keyval);
  for (i=0; i<2; i++) {
      if (!dir_[i]) continue;
      char *path = new char[strlen(dir_[i]) + strlen(filename) + 5];
      strcpy(path, dir_[i]);
      strcat(path, "/");
      strcat(path, filename);
      strcat(path, ".kv");

      // test to see if the file can be opened read only.
      FILE *fp = fopen(path, "r");
      if (fp) {
          fclose(fp);
          RefKeyVal libkeyval = new ParsedKeyVal(path);
          newkeyval = new AggregateKeyVal(keyval,libkeyval);
          delete[] path;
          break;
        }
      delete[] path;
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
// eval: (c-set-style "CLJ")
// End:
