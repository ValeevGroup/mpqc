//
// ipv2_error.cc
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

#include <stdlib.h>
#include <stdarg.h>
#include <string.h>

#include <util/misc/formio.h>
#include <util/keyval/ipv2.h>

using namespace std;
using namespace sc;

/* Returns some text for an errcod. */
const char*
IPV2::error_message(IPV2::Status errcod)
{
  const char *ipe_ok = "No problem has been detected.";
  const char *ipe_key_not_found = "No match was found for the given keyword.";
  const char *ipe_out_of_bounds = "An array index is out of bounds.";
  const char *ipe_malloc = "Memory allocation failed.";
  const char *ipe_not_an_array = "An index was given for a scalar quantity.";
  const char *ipe_not_a_scalar = "Expected a scalar, but found an array.";
  const char *ipe_type = "The datum is not of the appropiate type.";
  const char *ipe_has_no_value = "The keyword has no value.";
  const char *ipe_val_not_expd = "A value was not expected for the keyword.";
  const char *huh = "The nature of the problem is unknown.";
  
  if (errcod == OK) return ipe_ok;
  if (errcod == KeyNotFound) return ipe_key_not_found;
  if (errcod == OutOfBounds) return ipe_out_of_bounds;
  if (errcod == Malloc) return ipe_malloc;
  if (errcod == NotAnArray) return ipe_not_an_array;
  if (errcod == NotAScalar) return ipe_not_a_scalar;
  if (errcod == Type) return ipe_type;
  if (errcod == HasNoValue) return ipe_has_no_value;
  if (errcod == ValNotExpd) return ipe_val_not_expd;
  return huh;
}

void
IPV2::error(const char *msg)
{
  ExEnv::errn() << "IPV2::error: ";
  ExEnv::errn() << msg;
  ExEnv::errn() << endl;
  showpos();
  exit(1);
}

void
IPV2::warn(const char *msg)
{
  char *newmsg;
  const char *poskey;
  
  /* If msg has a %k in it, then substitute in the last keyword. */
#if defined(NCUBE)||defined(DEC)||defined(I860)
  /* The NCUBE and DEC are missing the strstr function. */
  for (poskey=msg; *poskey!='\0'; poskey++) {
      if (poskey[0] == '%' && poskey[1] == 'k') break;
    }
  if (poskey[0] != '%') poskey = NULL;
#else
  poskey = ::strstr(msg,"%k");
#endif
  if (poskey) {
      newmsg = (char *) malloc(strlen(msg)-1 + strlen(lastkeyword));
      strcpy(newmsg,msg);
      newmsg[poskey-msg] = '\0';
      strcat(newmsg,lastkeyword);
      strcat(newmsg,&poskey[2]);
      ExEnv::errn() << "IPV2::warn: ";
      ExEnv::errn() << newmsg;
      ExEnv::errn() << endl;
      if (poskey) free(newmsg);
    }
  else {
      ExEnv::errn() << "IPV2::warn: ";
      ExEnv::errn() << msg;
      ExEnv::errn() << endl;
    }
}

void
IPV2::ip_lastkeyword(const char*keyword)
{
  strncpy(lastkeyword,keyword, KEYWORD_LENGTH-1);
}

void
IPV2::ip_lastkeywordtree(ip_keyword_tree_t*kt)
{
  lastkeyword[0] = '\0';
  ip_lastkeyword_(kt);
}

void
IPV2::ip_lastkeyword_(ip_keyword_tree_t*kt)
{
  if (kt->up) ip_lastkeyword_(kt->up);
  if (strlen(lastkeyword) + strlen(kt->keyword) + 2 > KEYWORD_LENGTH) {
      ExEnv::errn() << "IPV2: keyword too big" << endl;
      abort();
    }
  strcat(lastkeyword,":");
  strcat(lastkeyword,kt->keyword);
}

/////////////////////////////////////////////////////////////////////////////

// Local Variables:
// mode: c++
// c-file-style: "CLJ"
// End:
