TOPDIR = ../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = libarray
BIN_OR_LIB = LIB

ifeq ($(ARCH),SUN4)
override COPTIONS = -O2
endif

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
default:: libarray.a

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d)
endif

