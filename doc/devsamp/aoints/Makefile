# Change this to the path to your installed sc-config script.
SCCONFIG = /usr/local/mpqc/current/bin/sc-config
CXX := $(shell $(SCCONFIG) --cxx)
CXXFLAGS := $(shell $(SCCONFIG) --cxxflags)
CC := $(shell $(SCCONFIG) --cc)
CCFLAGS := $(shell $(SCCONFIG) --cflags)
F77 := $(shell $(SCCONFIG) --f77)
F90 := $(F77)
FFLAGS := $(shell $(SCCONFIG) --f77flags)
F90FLAGS := $(FFLAGS)
CPPFLAGS := $(shell $(SCCONFIG) --cppflags)
LIBS := $(shell $(SCCONFIG) --libs)
LIBDIR  := $(shell $(SCCONFIG) --libdir)
LTCOMP := $(shell $(SCCONFIG) --ltcomp)
LTLINK := $(shell $(SCCONFIG) --ltlink)
LTLINKBINOPTS := $(shell $(SCCONFIG) --ltlinkbinopts)

aoints: aoints.o
	$(LTLINK) $(CXX) $(CXXFLAGS) -o $@ $^ -L$(LIBDIR) $(LIBS) $(LTLINKBINOPTS)

c_aoints: c_aoints.o
	$(LTLINK) $(CXX) $(CXXFLAGS) -o $@ $^ -L$(LIBDIR) -lSCexternqc $(LIBS) $(LTLINKBINOPTS)

f90_aoints: f90_aoints.o
	$(LTLINK) $(CXX) $(CXXFLAGS) -o $@ $^ -L$(LIBDIR) -lSCexternqc $(LIBS) $(LTLINKBINOPTS)

f90_aoints.o: f90_aoints.f90
	$(LTCOMP) $(F90) $(F90FLAGS) $(CPPFLAGS) -c $<

clean:
	-rm -f aoints.o aoints c_aoints.o c_aoints f90_aoints.o f90_aoints
