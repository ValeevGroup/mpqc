TOPDIR=../../..
ifndef SRCDIR
  BUILDING_IN_SRCDIR=yes
  SRCDIR=$(shell pwd)
endif

LOCALMAKEFILE_OPTIONAL = yes
include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

ifeq ($(LOCALMAKEFILE_FOUND),yes)

include $(TOPDIR)/lib/Makedirlist
INCLUDE += -I. -I$(SRCDIR) -I$(SRCTOPDIR)/src/bin/mpqc
CXXINCLUDE += -I. -I$(SRCDIR)

CXXSRC = pt2r12.cc moinfo.cc extern_pt2r12.cc

TARGET_TO_MAKE = pt2r12
default:: $(TARGET_TO_MAKE)

.PHONY: gamess
gamess : $(TARGET_TO_MAKE) 
gamess: CPPFLAGS += -DPT2R12GAMESS

.PHONY: gamesslib
gamesslib : libpt2r12gamess.$(LIBSUF)
gamesslib: CPPFLAGS += -DPT2R12GAMESS

INC =
DEPENDINCLUDE = $(INC)

BINOBJ = $(CXXSRC:%.cc=%.$(OBJSUF))

DISTFILES = $(CXXSRC) $(CSRC) $(INC) Makefile LIBS.h

LIBS :=  $(TOPDIR)/lib/libmpqcinit.a $(shell $(LISTLIBS) $(INCLUDE) $(DEFINES) $(SRCDIR)/LIBS.h)

#################################################################

default:: $(DEPENDINCLUDE)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(TARGET_TO_MAKE):$(BINOBJ) $(LIBS)
	echo "$(LISTLIBS) $(INCLUDE) $(DEFINES) $(XDEF) $(SRCDIR)/LIBS.h"
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(XLIBS) $(SYSLIBS) $(LTLINKBINOPTS)

BINOBJ_LIBGAMESS = pt2r12_not_main.o moinfo.o extern_pt2r12.o
libpt2r12gamess.$(LIBSUF): $(BINOBJ_LIBGAMESS)
	$(AR) $(ARFLAGS) $@ $^
	$(RANLIB) $@

pt2r12_not_main.o: pt2r12.cc
	$(CXX) $(CPPFLAGS) -DSKIP_MAIN $(CXXFLAGS) -c pt2r12.cc -o $@

ifeq ($(ENABLE_DEPLIBS),yes)
MPQCLTDEPLIBS=$(LIBS)
MPQCLTSYSLIBS=$(SYSLIBS)
else
MPQCLTDEPLIBS=
MPQCLTSYSLIBS=
endif

install:: $(TARGET_TO_MAKE)
	test -d $(installroot)$(bindir) || $(MKINSTALLDIRS) $(installroot)$(bindir)
	$(LTINST) $(INSTALL) $(INSTALLBINOPT) $< $(installroot)$(bindir)

clean::
	/bin/rm -f $(TARGET_TO_MAKE)

distclean:: clean

#################################################################

$(CXXSRC:.cc=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(CXXSRC:.cc=.d)
endif

endif
