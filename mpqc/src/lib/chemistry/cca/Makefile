TOPDIR=../../../..
ifndef SRCDIR
  SRCDIR=$(shell pwd)
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalMakefile
include $(TOPDIR)/lib/Makedirlist

SUBDIRS = client int server 
ifeq ($(ENABLESHARED),yes)
  SUBDIRS += module
endif

include $(SRCDIR)/$(TOPDIR)/lib/GlobalSubDirs

