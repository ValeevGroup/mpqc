TOPDIR=../../../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

TARGET_TO_MAKE = libmbpt
BIN_OR_LIB = LIB

CXXSRCS = mp2grad.cc cphf.cc ffo.cc opt2.cc opt2_v1.cc opt2_v2.cc opt2v2lb.cc \
          gmat.cc \
          mbpt.cc csgrad.cc csgmat.cc cscphf.cc \
          hsosv1.cc
          # hsosv1.cc hsosv2.cc hsosv2lb.cc

CSRCS = opt2_fock.c bzerofast.c gcollect.c

LIBOBJ= $(CXXSRCS:%.cc=%.$(OBJSUF)) $(CSRCS:%.c=%.$(OBJSUF))

default:: $(DEPENDINCLUDE)
default:: $(TOPDIR)/lib/libmbpt.a
ifeq ($(ELF),yes)
default:: $(TOPDIR)/lib/libmbpt.so
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d)
endif

LIBS = $(shell $(LISTLIBS) $(INCLUDE) $(SRCDIR)/LIBS.h)

mbpttest:: mbpttest.o $(LIBS)
	$(CXX) $(LDFLAGS) -o mbpttest $^ $(SYSLIBS)

mbpttest.o: mbpttest.cc
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -DSRCDIR=\"$(SRCDIR)\" -c $<
