TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=.
endif
include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

BIN_OR_LIB = LIB
TARGET_TO_MAKE = newmat7

#CXXFLAGS = $(CXXFLAGS_DEBUG)

#
OBJ = example.o

LIBSRC  = inlines.cc \
        cholesky.cc mcholesky.cc evalue.cc fft.cc hholder.cc jacobi.cc \
        newmat1.cc newmat2.cc newmat3.cc newmat4.cc newmat5.cc  \
        newmat6.cc newmat7.cc newmat8.cc newmatrm.cc           \
        sort.cc submat.cc svd.cc bandmat.cc except.cc newmatex.cc \
        print.cc cblas.cc linalg.c

LIBOBJ = $(LIBSRC:%.cc=%.o)
LIBOBJ := $(LIBOBJ:%.c=%.o)

default:: newmat7.a

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(OBJ:.o=.d) $(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(OBJ:.o=.d) $(LIBOBJ:.o=.d)
endif

example:      $(OBJ) newmat7.a
	      $(CXX) $(CXXFLAGS) -o $@ $(OBJ) newmat7.a -lm

clean::
	/bin/rm -f $(LIBOBJ) $(OBJ) example newmat7.a core

