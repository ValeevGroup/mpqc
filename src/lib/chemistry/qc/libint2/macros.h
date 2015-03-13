//
// macros.h
//
// Copyright (C) 2007 Edward Valeev
//
// Author: Edward Valeev <evaleev@vt.edu>
// Maintainer: EV
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

#include <libint2.h>

// see libint2's configure for the hardwired values
# if LIBINT2_CGSHELL_ORDERING == LIBINT2_CGSHELL_ORDERING_STANDARD
#  include <chemistry/qc/basis/macros_cca.h>
# elif LIBINT2_CGSHELL_ORDERING == LIBINT2_CGSHELL_ORDERING_INTV3
#  include <chemistry/qc/intv3/macros.h>
# elif LIBINT2_CGSHELL_ORDERING == LIBINT2_CGSHELL_ORDERING_GAMESS
#  include <chemistry/qc/libint2/macros_gamess.h>
# else
#  error "This version of Libint2 uses unsupported ordering of functions in shells"
# endif
