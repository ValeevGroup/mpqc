TOPDIR=../../..
ifndef SRCDIR
  SRCDIR=.
endif
include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

REQUIREDCXXSRC = mstate.cc \
         message.cc messagei.cc messproc.cc messint.cc \
         reduce.cc topology.cc hcube.cc picl.cc \
         rnglock.cc memory.cc memmsg.cc pool.cc memoryi.cc \
         globcnt.cc memamsg.cc memproc.cc memmid.cc

OPTIONALCXXSRC = messshm.cc messpgon.cc messpvm.cc messmpi.cc memshm.cc \
                 mempgon.cc memipgon.cc memmpl.cc

CXXSRC := $(REQUIREDCXXSRC)
ifneq ($(ARCH),SUNMOS)
ifneq ($(ARCH),PARAGON)
  CXXSRC += messshm.cc memshm.cc
endif
endif
ifeq ($(ARCH),PARAGON)
  CXXSRC += messpgon.cc mempgon.cc memipgon.cc
endif
ifeq ($(ARCH),SUNMOS)
  CXXSRC += messpgon.cc memipgon.cc
endif
ifeq ($(HAVE_PVM),yes)
  CXXSRC += messpvm.cc
endif
ifeq ($(HAVE_MPI),yes)
  CXXSRC += messmpi.cc
endif
ifeq ($(HAVE_MPL),yes)
  CXXSRC += memmpl.cc
endif

INC = pool.h message.h mstate.h topology.h picl.h \
      messshm.h messpvm.h messmpi.h \
      rnglock.h memory.h memmsg.h memshm.h pool.h memamsg.h memproc.h \
      memmpl.h

LIBOBJ = $(CXXSRC:%.cc=%.o)

BIN_OR_LIB = LIB
TARGET_TO_MAKE = libgroup

TESTCXXSRC = pooltest.cc messtest.cc rnglocktest.cc
DISTFILES = $(REQUIREDCXXSRC) $(OPTIONALCXXSRC) $(INC) \
            Makefile LIBS.h $(TESTCXXSRC)

DEPENDINCLUDE =

default:: $(DEPENDINCLUDE)
default:: $(TOPDIR)/lib/libgroup.a
ifeq ($(ELF),yes)
default:: $(TOPDIR)/lib/libgroup.so
endif

memtest.o: memtest.cc
	$(CXX) $(CXXFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<

memtest: memtest.o \
	libgroup.a \
	libkeyval.a \
	libstate.a \
	libclass.a \
	libcontainer.a \
	libref.a
	$(LD) $(LDFLAGS) -o memtest $^ $(SYSLIBS)

pooltest: pooltest.o libgroup.a
	$(CXX) $(CXXFLAGS) -o pooltest $^

rnglocktest: rnglocktest.o libgroup.a
	$(CXX) $(CXXFLAGS) -o rnglocktest $^

messtest.o: messtest.cc
	$(CXX) $(CXXFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<

messtest: messtest.o \
	libgroup.a \
	libkeyval.a \
	libmisc.a \
	libstate.a \
	libclass.a \
	libcontainer.a \
	libref.a
	$(LD) $(LDFLAGS) -o messtest $^ $(SYSLIBS)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(TESTCXXSRC:%.cc=%.d) $(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(TESTCXXSRC:%.cc=%.d) $(LIBOBJ:.o=.d)
endif

