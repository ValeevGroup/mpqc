//
// state_file.cc
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

#include <util/class/class.h>
#include <util/state/state.h>

#include <util/state/statenumSet.h>
#include <util/state/classdImplMap.h>

StateOutFile::StateOutFile() :
  opened_(0), fp_(stdout)
{
}

StateOutFile::StateOutFile(FILE* fp) :
  opened_(0), fp_(fp)
{
}

StateOutFile::StateOutFile(const char * path, const char * mode) :
  opened_(1)
{
  fp_ = fopen(path,mode);
}

StateOutFile::~StateOutFile()
{
  if (opened_) close();
}

void StateOutFile::flush() { fflush(fp_); }
void StateOutFile::close()
{
  if(opened_) fclose(fp_);
  opened_=0; fp_=0;

  _classidmap->clear(); _nextclassid=0;
  forget_references();
}

void StateOutFile::rewind() { if(fp_) fseek(fp_,0,0); }

int StateOutFile::open(const char *path, const char * mode)
{
  if (opened_) close();

  if ((fp_ = fopen(path,mode))==0) {
      cerr << "StateOutFile::open(" << path << "," << mode << ") failed"
           << endl;
      return -1;
    }

  opened_ = 1;
  return 0;
}

////////////////////////////////////

StateInFile::StateInFile() :
  opened_(0), fp_(stdin)
{
}

StateInFile::StateInFile(FILE* fp) :
  opened_(0), fp_(fp)
{
}

StateInFile::StateInFile(const char * path, const char * mode) :
  opened_(1)
{
  fp_ = fopen(path,mode);
}

StateInFile::~StateInFile()
{
  if (opened_) close();
}

void StateInFile::flush() { fflush(fp_); }
void StateInFile::close()
{
  if(opened_) fclose(fp_);
  opened_=0; fp_=0;

  _cd.clear();
  forget_references();
}
void StateInFile::rewind() { if(fp_) fseek(fp_,0,0); }

int StateInFile::open(const char *path, const char * mode)
{
  if (opened_) close();

  if ((fp_ = fopen(path,mode))==0) {
      cerr << "StateInFile::open(" << path << "," << mode << ") failed"
           << endl;
      return -1;
    }

  opened_ = 1;
  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
