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

#ifdef __GNUC__
#pragma implementation
#endif

#include <fstream>

#include <util/state/state_file.h>

using namespace std;
using namespace sc;

static ClassDesc StateOutFile_cd(
    typeid(StateOutFile),"StateOutFile",1,"public StateOut");

StateOutFile::StateOutFile() :
  opened_(0), buf_(ExEnv::outn().rdbuf())
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

  classidmap_.clear();
  nextclassid_=0;

  ps_.clear();
  next_object_number_ = 1;
}

int StateOutFile::open(const char *path)
{
  if (opened_) close();

  filebuf *fbuf = new filebuf();
  fbuf->open(path, ios::out);
  if (!fbuf->is_open()) {
      ExEnv::errn() << "ERROR: StateOutFile: problems opening " << path << endl;
      abort();
    }
  buf_ = fbuf;

  opened_ = 1;
  return 0;
}

////////////////////////////////////

static ClassDesc StateInFile_cd(
    typeid(StateInFile),"StateInFile",1,"public StateIn");

StateInFile::StateInFile() :
  opened_(0), buf_(cin.rdbuf())
{
}

StateInFile::StateInFile(istream& s) :
  opened_(0), buf_(s.rdbuf())
{
}

StateInFile::StateInFile(const char * path)
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

  classidmap_.clear();
  nextclassid_ = 0;
  classdatamap_.clear();
  ps_.clear();
}

int StateInFile::open(const char *path)
{
  if (opened_) close();

  filebuf *fbuf = new filebuf();
  fbuf->open(path, ios::in);
  if (!fbuf->is_open()) {
      ExEnv::errn() << "ERROR: StateInFile: problems opening " << path << endl;
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
