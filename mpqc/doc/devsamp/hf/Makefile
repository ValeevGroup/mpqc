# Change this to the path to your installed sc-config script.
SCCONFIG = /usr/local/mpqc/current/bin/sc-config
CXX := $(shell $(SCCONFIG) --cxx)
CXXFLAGS := $(shell $(SCCONFIG) --cxxflags)
CPPFLAGS := $(shell $(SCCONFIG) --cppflags)
LIBS := $(shell $(SCCONFIG) --libs)
LIBDIR  := $(shell $(SCCONFIG) --libdir)
LTLINK := $(shell $(SCCONFIG) --ltlink)
LTLINKBINOPTS := $(shell $(SCCONFIG) --ltlinkbinopts)

hf: hf.o
	$(LTLINK) $(CXX) $(CXXFLAGS) -o $@ $^ -L$(LIBDIR) -lmpqcinit $(LIBS) $(LTLINKBINOPTS)

clean:
	/bin/rm -f hf hf.o hf.wfn hf.wfn.tmp
