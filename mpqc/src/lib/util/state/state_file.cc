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

#include <fstream.h>

#include <util/class/class.h>
#include <util/state/state.h>

#include <util/state/statenumSet.h>
#include <util/state/classdImplMap.h>
#include <util/state/stateptrImplSet.h>

StateOutFile::StateOutFile() :
  opened_(0), buf_(cout.rdbuf())
{
}

StateOutFile::StateOutFile(ostream& s) :
  opened_(0), buf_(s.rdbuf())
{
}

StateOutFile::StateOutFile(const char * path)
{
  opened_ = 0;
  open(path);
}

StateOutFile::~StateOutFile()
{
  close();
}

void StateOutFile::flush()
{
  // ostream needed due to out-of-date streambuf implementations
  ostream o(buf_);
  o.flush();
}
void StateOutFile::close()
{
  if(opened_) delete buf_;
  opened_=0; buf_=0;

  _classidmap->clear();
  _nextclassid=0;

  ps_->clear();
  next_pointer_number = 1;
}

void StateOutFile::rewind() { if(buf_) buf_->seekoff(0,ios::beg); }

int StateOutFile::open(const char *path)
{
  if (opened_) close();

  filebuf *fbuf = new filebuf();
  fbuf->open(path, ios::out);
  if (!fbuf->is_open()) {
      cerr << "ERROR: StateOutFile: problems opening " << path << endl;
      abort();
    }
  buf_ = fbuf;

  opened_ = 1;
  return 0;
}

////////////////////////////////////

StateInFile::StateInFile() :
  opened_(0), buf_(cin.rdbuf())
{
}

StateInFile::StateInFile(istream& s) :
  opened_(0), buf_(s.rdbuf())
{
}

StateInFile::StateInFile(const char * path) :
  opened_(1)
{
  opened_ = 0;
  open(path);
}

StateInFile::~StateInFile()
{
  close();
}

void StateInFile::close()
{
  if(opened_) delete buf_;
  opened_=0; buf_=0;

  _cd.clear();
  ps_->clear();
}
void StateInFile::rewind() { if(buf_) buf_->seekoff(0,ios::beg); }

int StateInFile::open(const char *path)
{
  if (opened_) close();

  filebuf *fbuf = new filebuf();
  fbuf->open(path, ios::in);
  if (!fbuf->is_open()) {
      cerr << "ERROR: StateInFile: problems opening " << path << endl;
      abort();
    }
  buf_ = fbuf;

  opened_ = 1;
  return 0;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
