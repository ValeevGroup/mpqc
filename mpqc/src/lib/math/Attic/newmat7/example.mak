
OBJ  =  example.o                                          \
        cholesky.o evalue.o fft.o hholder.o jacobi.o       \
        newmat1.o newmat2.o newmat3.o newmat4.o newmat5.o  \
        newmat6.o newmat7.o newmat8.o newmatrm.o           \
        sort.o submat.o svd.o bandmat.o except.o newmatex.o

example:      $(OBJ)
	      g++ -o $@ $(OBJ) -lm

%.o:          %.cxx
	      g++ -D__GNUG__ -c $*.cxx

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












