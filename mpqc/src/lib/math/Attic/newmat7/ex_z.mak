C = ztc -c -ml $*.cxx

OBJ = fft.obj evalue.obj submat.obj cholesky.obj hholder.obj          \
  sort.obj newmatrm.obj jacobi.obj svd.obj example.obj                \
  newmat8.obj newmat7.obj newmat6.obj                                 \
  newmat5.obj newmat3.obj newmat4.obj newmat2.obj newmat1.obj         \
  bandmat.obj except.obj newmatex.obj

ex_z.exe:       $(OBJ) ex_z.lnk
                blink @ex_z.lnk

ex_z.lnk:       ex_z.mak
	        echo newmat1.obj+newmat2.obj+newmat3.obj+    > ex_z.lnk
	        echo newmat4.obj+svd.obj+newmat5.obj+       >> ex_z.lnk
	        echo newmat6.obj+newmat7.obj+newmat8.obj+   >> ex_z.lnk
	        echo cholesky.obj+hholder.obj+sort.obj+     >> ex_z.lnk
	        echo submat.obj+jacobi.obj+newmatrm.obj+    >> ex_z.lnk
	        echo fft.obj+evalue.obj+bandmat.obj+        >> ex_z.lnk
	        echo newmatex.obj+except.obj+example.obj    >> ex_z.lnk
	        echo ex_z.exe                               >> ex_z.lnk

newmatxx:       include.h newmat.h boolean.h except.h
	        echo "main .h files uptodate?" > newmatxx

except.obj:     except.h except.cxx
	        $C

newmatex.obj:   newmatxx newmatex.cxx
	        $C

example.obj:    newmatxx newmatap.h example.cxx
	        $C

cholesky.obj:   newmatxx cholesky.cxx
	        $C

evalue.obj:     newmatxx newmatrm.h precisio.h evalue.cxx
	        $C

fft.obj:        newmatxx newmatap.h fft.cxx
	        $C

hholder.obj:    newmatxx newmatap.h hholder.cxx
	        $C

jacobi.obj:     newmatxx precisio.h newmatrm.h jacobi.cxx
	        $C

bandmat.obj:    newmatxx newmatrc.h controlw.h bandmat.cxx
	        $C

newmat1.obj:    newmatxx newmat1.cxx
	        $C

newmat2.obj:    newmatxx newmatrc.h controlw.h newmat2.cxx
	        $C

newmat3.obj:    newmatxx newmatrc.h controlw.h newmat3.cxx
	        $C

newmat4.obj:    newmatxx newmatrc.h controlw.h newmat4.cxx
	        $C

newmat5.obj:    newmatxx newmatrc.h controlw.h newmat5.cxx
	        $C

newmat6.obj:    newmatxx newmatrc.h controlw.h newmat6.cxx
	        $C

newmat7.obj:    newmatxx newmatrc.h controlw.h newmat7.cxx
	        $C

newmat8.obj:    newmatxx newmatap.h newmat8.cxx
	        $C

newmat9.obj:    newmatxx newmatrc.h controlw.h newmatio.h newmat9.cxx
	        $C

newmatrm.obj:   newmatxx newmatrm.h newmatrm.cxx
	        $C

sort.obj:       newmatxx newmatap.h sort.cxx
	        $C

submat.obj:     newmatxx newmatrc.h controlw.h submat.cxx
	        $C

svd.obj:        newmatxx newmatrm.h precisio.h svd.cxx
	        $C

example.obj:    newmatxx newmatap.h example.cxx 
	        $C

