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

StateOutFile::StateOutFile() :
  opened_(0), buf_(cout.rdbuf()), stream_(cout.rdbuf())
{
}

StateOutFile::StateOutFile(ostream& s) :
  opened_(0), buf_(s.rdbuf()), stream_(s.rdbuf())
{
}

static streambuf*
getoutbuf(const char *path)
{
  filebuf *fbuf = new filebuf();
  fbuf->open(path,ios::out);
  if (!fbuf->is_open()) {
      cerr << "ERROR: StateOutFile: problems opening " << path << endl;
      abort();
    }
  return fbuf;
}

StateOutFile::StateOutFile(const char * path) :
  opened_(1), stream_(getoutbuf(path))
{
  buf_ = stream_.rdbuf();
}

StateOutFile::~StateOutFile()
{
  if (opened_) {
      delete buf_;
    }
}

void StateOutFile::flush() { stream_.flush(); }
void StateOutFile::close()
{
  if(opened_) delete buf_;
  opened_=0; buf_=0;

  _classidmap->clear(); _nextclassid=0;
  forget_references();
}

void StateOutFile::rewind() { if(buf_) buf_->seekoff(0,ios::beg); }

int StateOutFile::open(const char *path)
{
//  if (opened_) close();
//
//  filebuf *fbuf = new filebuf();
//  fbuf->open(path, ios::out);
//  if (!fbuf->is_open()) {
//      cerr << "ERROR: StateOutFile: problems opening " << path << endl;
//      abort();
//    }
//  buf_ = fbuf;
//  stream_.rdbuf(buf_);
//
//  opened_ = 1;
//  return 0;
  cerr << "StateOutFile::open is disabled to maintain compatibility" << endl
       << "with older C++ iostream libraries." << endl;
  abort();
  return -1;
}

////////////////////////////////////

StateInFile::StateInFile() :
  opened_(0), buf_(cin.rdbuf()), stream_(cin.rdbuf())
{
}

StateInFile::StateInFile(istream& s) :
  opened_(0), buf_(s.rdbuf()), stream_(s.rdbuf())
{
}

static streambuf*
getinbuf(const char *path)
{
  filebuf *fbuf = new filebuf();
  fbuf->open(path,ios::in);
  if (!fbuf->is_open()) {
      cerr << "ERROR: StateInFile: problems opening " << path << endl;
      abort();
    }
  return fbuf;
}

StateInFile::StateInFile(const char * path) :
  opened_(1), stream_(getinbuf(path))
{
  buf_ = stream_.rdbuf();
}

StateInFile::~StateInFile()
{
  if (opened_) {
      delete buf_;
    }
}

void StateInFile::close()
{
  if(opened_) delete buf_;
  opened_=0; buf_=0;

  _cd.clear();
  forget_references();
}
void StateInFile::rewind() { if(buf_) buf_->seekoff(0,ios::beg); }

int StateInFile::open(const char *path)
{
//  if (opened_) close();
//
//  filebuf *fbuf = new filebuf();
//  fbuf->open(path, ios::in);
//  if (!fbuf->is_open()) {
//      cerr << "ERROR: StateInFile: problems opening " << path << endl;
//      abort();
//    }
//  buf_ = fbuf;
//  stream_.rdbuf(buf_);
//
//  opened_ = 1;
//  return 0;
  cerr << "StateInFile::open is disabled to maintain compatibility" << endl
       << "with older C++ iostream libraries." << endl;
  abort();
  return -1;
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
