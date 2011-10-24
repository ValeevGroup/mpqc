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

TARGET_TO_MAKE = pt2r12
CXXSRC = pt2r12.cc moinfo.cc \
pt2r12_molcas.cc \
#pt2r12_gamess.cc \
 

INC =
DEPENDINCLUDE = $(INC)

BINOBJ = $(CXXSRC:%.cc=%.$(OBJSUF))

DISTFILES = $(CXXSRC) $(CSRC) $(INC) Makefile LIBS.h

default:: $(TARGET_TO_MAKE)

LIBS :=  $(TOPDIR)/lib/libmpqcinit.a $(shell $(LISTLIBS) $(INCLUDE) $(DEFINES) $(SRCDIR)/LIBS.h)

#################################################################

default:: $(DEPENDINCLUDE)

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(TARGET_TO_MAKE): $(BINOBJ) $(LIBS)
	echo "$(LISTLIBS) $(INCLUDE) $(DEFINES) $(XDEF) $(SRCDIR)/LIBS.h"
	$(LTLINK) $(LD) $(LDFLAGS) -o $@ $^ $(XLIBS) $(SYSLIBS) $(LTLINKBINOPTS)

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

$(BINOBJ:.$(OBJSUF)=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(BINOBJ:.$(OBJSUF)=.d)
endif

endif
