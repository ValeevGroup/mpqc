#
CPLUSPLUS = g++
CPLUSPLUSFLAGS = 
OBJ =   tmt.o tmt1.o tmt2.o tmt3.o tmt4.o tmt5.o tmt6.o tmt7.o tmt8.o          \
  tmt9.o tmta.o tmtb.o tmtc.o tmtd.o tmtf.o tmtg.o tmth.o tmti.o

tmt:      $(OBJ) ../newmat7.a
	      $(CPLUSPLUS) $(CPLUSPLUSFLAGS) -o $@ $(OBJ) ../newmat7.a -lm

%.o:          %.cxx
	      $(CPLUSPLUS) $(CPLUSPLUSFLAGS) -c $*.cxx

newmatxx:     include.h ../include.h newmat.h ../newmat.h boolean.h \
		../boolean.h except.h ../except.h
	      rm -f newmatxx
	      echo "main .h files uptodate?" > newmatxx

tmt.o:        newmatxx newmatap.h ../newmatap.h tmt.cxx 

tmt1.o:       newmatxx newmatap.h ../newmatap.h tmt1.cxx 

tmt2.o:       newmatxx newmatap.h ../newmatap.h tmt2.cxx 

tmt3.o:       newmatxx newmatap.h ../newmatap.h tmt3.cxx 

tmt4.o:       newmatxx newmatap.h ../newmatap.h tmt4.cxx 

tmt5.o:       newmatxx newmatap.h ../newmatap.h tmt5.cxx 

tmt6.o:       newmatxx newmatap.h ../newmatap.h tmt6.cxx 

tmt7.o:       newmatxx newmatap.h ../newmatap.h tmt7.cxx 

tmt8.o:       newmatxx newmatap.h ../newmatap.h tmt8.cxx 

tmt9.o:       newmatxx newmatap.h ../newmatap.h tmt9.cxx 

tmta.o:       newmatxx newmatap.h ../newmatap.h tmta.cxx 

tmtb.o:       newmatxx newmatap.h ../newmatap.h tmtb.cxx 

tmtc.o:       newmatxx newmatap.h ../newmatap.h tmtc.cxx 

tmtd.o:       newmatxx newmatap.h ../newmatap.h tmtd.cxx 

tmte.o:       newmatxx newmatap.h ../newmatap.h tmte.cxx 

tmtf.o:       newmatxx newmatap.h ../newmatap.h tmtf.cxx 

tmtg.o:       newmatxx newmatap.h ../newmatap.h tmtg.cxx 

tmth.o:       newmatxx newmatap.h ../newmatap.h tmth.cxx

tmti.o:       newmatxx newmatap.h ../newmatap.h tmti.cxx

.cxx.o:
	$(CPLUSPLUS) $(CPLUSPLUSFLAGS) -g -c $<

clean:
	/bin/rm -f $(OBJ) tmt core

