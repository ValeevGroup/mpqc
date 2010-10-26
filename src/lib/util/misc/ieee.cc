/*
 * ieee.c
 *
 * Copyright (C) 1996 Limit Point Systems, Inc.
 *
 * Author: Curtis Janssen <cljanss@ca.sandia.gov>
 * Maintainer: LPS
 *
 * This file is part of the SC Toolkit.
 *
 * The SC Toolkit is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Library General Public License as published by
 * the Free Software Foundation; either version 2, or (at your option)
 * any later version.
 *
 * The SC Toolkit is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public License
 * along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
 * the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 * The U.S. Government is granted a limited license as per AL 91-7.
 */

#include <util/misc/ieee.h>

#ifdef SGI
#include <sys/fpu.h>
namespace sc {
void
ieee_trap_errors()
{
  union fpc_csr fc;

  fc.fc_word = get_fpc_csr();
  fc.fc_struct.en_divide0 = 1;
  fc.fc_struct.en_invalid = 1;
  fc.fc_struct.en_overflow = 1;
  set_fpc_csr(fc.fc_word);
}
}

#else
namespace sc {
void
ieee_trap_errors()
{
}
} 
#endif
