TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

LD=$(CXX)

CXXSRC = volume.cc tess.cc dtess.cc shape.cc surf.cc isosurf.cc
CSRC = ge_error.c dirichlet.c gauss.c

LIBOBJ = $(CXXSRC:%.cc=%.o) $(CSRC:%.c=%.o)

INC = volume.h shape.h tess.h ge_error.h dirichlet.h pgen_objects.h

DEPENDINCLUDE = $(INC)

BIN_OR_LIB = LIB
TARGET_TO_MAKE = isosurf

DISTFILES = isotest.cc $(CXXSRC) $(CSRC) $(INC) Makefile LIBS.h

default:: $(DEPENDINCLUDE)
default:: isosurf.a

isotest: isotest.o isosurf.a util/options/liboptions.a \
	$(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)
	$(LD) $(LDFLAGS) -o isotest $^ $(SYSLIBS)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d) isotest.d
endif

