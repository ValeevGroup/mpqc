OBJ = fft.o evalue.o submat.o cholesky.o hholder.o sort.o newmatrm.o     \
  jacobi.o tmtf.o svd.o tmte.o tmtd.o newmat8.o tmtc.o tmtb.o            \
  newmat7.o newmat6.o newmat5.o newmat3.o newmat4.o newmat2.o newmat1.o  \
  tmt.o tmt1.o tmt2.o tmt3.o tmt4.o tmt5.o tmt6.o tmt7.o tmt8.o          \
  tmt9.o tmta.o tmtg.o tmth.o tmti.o bandmat.o except.o newmatex.o

tmt:          $(OBJ)
	      g++ -o $@ $(OBJ) -lm

%.o:          %.cxx
	      g++ -c $*.cxx

newmatxx:     include.h newmat.h boolean.h except.h
	      rm -f newmatxx
	      echo "main .h files uptodate?" > newmatxx

except.o:     except.h except.cxx

newmatex.o:   newmatxx newmatex.cxx

example.o:    newmatxx newmatap.h example.cxx

cholesky.o:   newmatxx cholesky.cxx

evalue.o:     newmatxx newmatrm.h precisio.h evalue.cxx

fft.o:        newmatxx newmatap.h fft.cxx

hholder.o:    newmatxx newmatap.h hholder.cxx

jacobi.o:     newmatxx precisio.h newmatrm.h jacobi.cxx

bandmat.o:    newmatxx newmatrc.h controlw.h bandmat.cxx

newmat1.o:    newmatxx newmat1.cxx

newmat2.o:    newmatxx newmatrc.h controlw.h newmat2.cxx

newmat3.o:    newmatxx newmatrc.h controlw.h newmat3.cxx

newmat4.o:    newmatxx newmatrc.h controlw.h newmat4.cxx

newmat5.o:    newmatxx newmatrc.h controlw.h newmat5.cxx

newmat6.o:    newmatxx newmatrc.h controlw.h newmat6.cxx

newmat7.o:    newmatxx newmatrc.h controlw.h newmat7.cxx

newmat8.o:    newmatxx newmatap.h newmat8.cxx

newmat9.o:    newmatxx newmatrc.h controlw.h newmatio.h newmat9.cxx

newmatrm.o:   newmatxx newmatrm.h newmatrm.cxx

sort.o:       newmatxx newmatap.h sort.cxx

submat.o:     newmatxx newmatrc.h controlw.h submat.cxx

svd.o:        newmatxx newmatrm.h precisio.h svd.cxx

tmt.o:        newmatxx newmatap.h tmt.cxx 

tmt1.o:       newmatxx newmatap.h tmt1.cxx 

tmt2.o:       newmatxx newmatap.h tmt2.cxx 

tmt3.o:       newmatxx newmatap.h tmt3.cxx 

tmt4.o:       newmatxx newmatap.h tmt4.cxx 

tmt5.o:       newmatxx newmatap.h tmt5.cxx 

tmt6.o:       newmatxx newmatap.h tmt6.cxx 

tmt7.o:       newmatxx newmatap.h tmt7.cxx 

tmt8.o:       newmatxx newmatap.h tmt8.cxx 

tmt9.o:       newmatxx newmatap.h tmt9.cxx 

tmta.o:       newmatxx newmatap.h tmta.cxx 

tmtb.o:       newmatxx newmatap.h tmtb.cxx 

tmtc.o:       newmatxx newmatap.h tmtc.cxx 

tmtd.o:       newmatxx newmatap.h tmtd.cxx 

tmte.o:       newmatxx newmatap.h tmte.cxx 

tmtf.o:       newmatxx newmatap.h tmtf.cxx 

tmtg.o:       newmatxx newmatap.h tmtg.cxx 

tmth.o:       newmatxx newmatap.h tmth.cxx

tmti.o:       newmatxx newmatap.h tmti.cxx


