#
# Makefile
#
# Copyright (C) 1996 Limit Point Systems, Inc.
#
# Author: Curtis Janssen <cljanss@ca.sandia.gov>
# Maintainer: LPS
#
# This file is part of the SC Toolkit.
#
# The SC Toolkit is free software; you can redistribute it and/or modify
# it under the terms of the GNU Library General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# The SC Toolkit is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General Public License
# along with the SC Toolkit; see the file COPYING.LIB.  If not, write to
# the Free Software Foundation, 675 Mass Ave, Cambridge, MA 02139, USA.
#
# The U.S. Government is granted a limited license as per AL 91-7.
#

TOPDIR=../../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = libarray
BIN_OR_LIB = LIB

SGENFLAGS = -ifndef _math_array_libmath_h

CSRC = diag.c flin.c matrix_utl.c mxmm.c mxmv.c print.c

SGENFILE = matrix.sgen

SGENSRC = matrixbc0.c matrixrbc0.c matrixsbc0.c \
          matrixfree.c matrixallc.c matrixinit.c \
          matrixprnt.c matrixiseq.c matrixasgn.c matrixsnd0.c \
          matrixrcv0.c matrixzero.c

INC = math_lib.h
SGENINC = $(SGENSRC:%.c=%.h)
GENINC = $(CSRC:%.c=%.gbl)

LIBSRC = $(CSRC) $(SGENSRC)

LIBOBJ = $(LIBSRC:%.c=%.$(OBJSUF))

DISTFILES = $(CSRC) $(SGENFILE) Makefile $(INC) LIBS.h
DEPENDINCLUDE = $(GENINC) $(SGENINC) $(INC)

############################################################3

default:: $(DEPENDINCLUDE)
interface:: $(DEPENDINCLUDE)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d)
endif

