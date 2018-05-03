//
// formio.cpp
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

#include "formio.h"
#include "exenv.h"

#include <cstdio> // for vsprintf
#include <cwchar> // for vswprintf
#include <cstdlib>
#include <cstring>
#include <cstdarg>

using namespace mpqc;

char *FormIO::default_basename_ = nullptr;
char *FormIO::default_work_dir_ = nullptr;
int  FormIO::ready_ = 0;
int  FormIO::xalloc_inited_ = 0;
long FormIO::nindent_ = 0;
long FormIO::indent_size_ = 0;
long FormIO::skip_indent_ = 0;
long FormIO::verbose_ = 0;
long FormIO::initialized_ = 0;
int FormIO::node_to_print_ = 0;
int FormIO::debug_ = 0;
int FormIO::parallel_ = 0;
int FormIO::me_ = 0;

char *
FormIO::fileext_to_filename(const char *ext)
{
  const char *basename;

  if (default_basename_) basename = default_basename_;
  else basename = "mpqc";

  char * res = new char[strlen(basename) + strlen(ext) + 1];
  strcpy(res, basename);
  strcat(res, ext);

  return res;
}

std::string
FormIO::fileext_to_filename_string(const char *ext)
{
  std::string basename;

  if (default_basename_) basename = default_basename_;
  else basename = "mpqc";

  std::string res = basename + ext;

  return res;
}

std::string
FormIO::fileext_to_fullpathname_string(const char *ext){

  std::string fullpathname;

  fullpathname = default_work_dir_;
  fullpathname += '/';
  fullpathname += default_basename_;
  fullpathname += ext;

  return fullpathname;

}

void
FormIO::set_default_basename(const char *basename)
{
  if (default_basename_) delete[] default_basename_;

  if (basename)
      default_basename_ = strcpy(new char[strlen(basename)+1], basename);
  else
      default_basename_ = 0;
}

void
FormIO::set_default_work_dir(const char* work_dir)
{
  if (default_work_dir_) delete[] default_work_dir_;

  if (work_dir)
    default_work_dir_ = strcpy(new char[strlen(work_dir)+1], work_dir);
  else
    default_work_dir_ = 0;
}

const char *
FormIO::default_basename()
{
  return default_basename_;
}

const char *
FormIO::default_work_dir()
{
  return default_work_dir_;
}

int
FormIO::set_printnode(int n)
{
  int r = node_to_print_;
  node_to_print_ = n;
  return r;
}

void
FormIO::set_debug(int n)
{
  debug_ = n;
}

void
FormIO::init_mp(int me)
{
  if (!ready_) init();
  me_ = me;
  parallel_=1;
}

void
FormIO::init_ostream(std::ostream &o)
{
  if (!xalloc_inited_) {
      xalloc_inited_ = 1;
      nindent_ = std::ios::xalloc();
      indent_size_ = std::ios::xalloc();
      skip_indent_ = std::ios::xalloc();
      verbose_ = std::ios::xalloc();
      initialized_ = std::ios::xalloc();
    }

  if (o.iword(initialized_)) return;

  o.iword(skip_indent_) = 0;
  o.iword(indent_size_) = 0;
  o.iword(nindent_) = 2;
  o.iword(verbose_) = 0;
  o.iword(initialized_) = 1;
}

void
FormIO::init()
{
  ready_ = 1;

  init_ostream(std::cout);
  init_ostream(std::cerr);
}

std::ios&
FormIO::indent(std::ios&o)
{
  if (!ready_) init();
  long &skip = o.iword(skip_indent_);
  if (skip) {
      skip--;
      return o;
    }
  if (debug_ && parallel_) {
      char nn[24];
      sprintf(nn,"node %5d:",me_);
      for (size_t i=0; i < strlen(nn); i++) o.rdbuf()->sputc(nn[i]);
    }
  long n = o.iword(nindent_);
  for (int i=0; i<n; i++) o.rdbuf()->sputc(' ');
  return o;
}

std::ios&
FormIO::incindent(std::ios&o)
{
  if (!ready_) init();
  long &n = o.iword(nindent_);
  long size = o.iword(indent_size_);
  if (size == 0) size = 2;
  else if (size < 0) size = 0;
  n += size;
  return o;
}

std::ios&
FormIO::decindent(std::ios&o)
{
  if (!ready_) init();
  long &n = o.iword(nindent_);
  long size = o.iword(indent_size_);
  if (size == 0) size = 2;
  else if (size < 0) size = 0;
  n -= size;
  if (n<0) n=0;
  return o;
}

long
FormIO::getindent(std::ios&o)
{
  if (!ready_) init();
  return o.iword(nindent_);
}

void
FormIO::setindent(std::ios&o, long n)
{
  if (!ready_) init();
  o.iword(nindent_) = n;
}

long
FormIO::getverbose(std::ios&o)
{
  if (!ready_) init();
  return o.iword(verbose_);
}

void
FormIO::setverbose(std::ios&o, long n)
{
  if (!ready_) init();
  o.iword(verbose_) = n;
}

std::ios&
FormIO::skipnextindent(std::ios&o)
{
  if (!ready_) init();
  o.iword(skip_indent_)++;
  return o;
}

std::ostream&
FormIO::copyright(std::ostream& o)
{
  o << indent
    << "Copyright (C) 2014-2017 Virginia Tech."
    << std::endl;

  return o;
}

std::ostream&
FormIO::license(std::ostream& o)
{
  o << indent
    << "This program is open-source software; you can redistribute it and/or modify"
    << std::endl << indent
    << "it under the terms of the GNU General Public License as published by"
    << std::endl << indent
    << "the Free Software Foundation; either version 3 of the License, or"
    << std::endl << indent
    << "(at your option) any later version."
    << std::endl;

  return o;
}

std::ostream&
FormIO::warranty(std::ostream& o)
{
  o << indent
    << "This program is distributed in the hope that it will be useful,"
    << std::endl << indent
    << "but WITHOUT ANY WARRANTY; without even the implied warranty of"
    << std::endl << indent
    << "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the"
    << std::endl << indent
    << "GNU General Public License for more details."
    << std::endl;

  return o;
}

std::ios&
mpqc::indent(std::ios& o)
{
  return FormIO::indent(o);
}

std::ios&
mpqc::decindent(std::ios& o)
{
  return FormIO::decindent(o);
}

std::ios&
mpqc::incindent(std::ios& o)
{
  return FormIO::incindent(o);
}

std::ios&
mpqc::skipnextindent(std::ios& o)
{
  return FormIO::skipnextindent(o);
}

/////////////////////////////////////////////////////////////////////////////

namespace mpqc {

template <>
mpqcprintf<char>::mpqcprintf(const char* fmt, ...) {
  va_list args;

  va_start(args, fmt);

  str_[0] = '\0';

  // hopefully this won't overflow
  if (fmt && fmt[0] != '\0') {
    if (std::vsprintf(str_, fmt, args) > 1023) {
      ExEnv::errn() << indent << "mpqcprintf overflow\n";
      abort();
    }
  }

  va_end(args);
}

template <>
mpqcprintf<wchar_t>::mpqcprintf(const wchar_t* fmt, ...) {
  va_list args;

  va_start(args, fmt);

  str_[0] = L'\0';

  // hopefully this won't overflow
  if (fmt && fmt[0] != '\0') {
    if (std::vswprintf(str_, sizeof(str_)/sizeof(wchar_t), fmt, args) > 1023) {
      ExEnv::errn() << indent << "mpqcprintf overflow\n";
      abort();
    }
  }

  va_end(args);
}
}  // namespace mpqc

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
