TOPDIR=../../../../..
ifndef SRCDIR
  SRCDIR=.
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile

#CXXFLAGS = $(CXXFLAGS_DEBUG)

TARGET_TO_MAKE = libintegral
BIN_OR_LIB = LIB

SRCS = integral.cc integralv2.cc ortho.cc intiter.cc

LIBOBJ= $(SRCS:%.cc=%.$(OBJSUF))

default:: $(DEPENDINCLUDE)
default:: $(TOPDIR)/lib/libintegral.a
ifeq ($(ELF),yes)
default:: $(TOPDIR)/lib/libintegral.so
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalRules

$(LIBOBJ:.o=.d): $(DEPENDINCLUDE)
ifneq ($(DODEPEND),no)
include $(LIBOBJ:.o=.d)
endif

