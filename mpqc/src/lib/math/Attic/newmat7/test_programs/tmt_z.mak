C = ztc -c -ml $*.cxx -o

OBJ = fft.obj evalue.obj submat.obj cholesky.obj hholder.obj          \
  sort.obj newmatrm.obj jacobi.obj tmtf.obj svd.obj tmte.obj          \
  tmtd.obj newmat8.obj tmtc.obj tmtb.obj newmat7.obj newmat6.obj      \
  newmat5.obj newmat3.obj newmat4.obj newmat2.obj newmat1.obj         \
  tmt.obj tmt1.obj tmt2.obj tmt3.obj tmt4.obj tmt5.obj tmt6.obj       \
  tmt7.obj tmt8.obj tmt9.obj tmta.obj tmtg.obj tmth.obj tmti.obj      \
  bandmat.obj except.obj newmatex.obj

tmt_z.exe:      $(OBJ) tmt_z.lnk
                blink @tmt_z.lnk

tmt_z.lnk:      tmt_z.mak
	        echo newmat1.obj+newmat2.obj+newmat3.obj+    > tmt_z.lnk
	        echo newmat4.obj+svd.obj+newmat5.obj+       >> tmt_z.lnk
	        echo newmat6.obj+newmat7.obj+newmat8.obj+   >> tmt_z.lnk
	        echo tmt.obj+tmt1.obj+tmt2.obj+tmt3.obj+    >> tmt_z.lnk
	        echo tmt4.obj+tmt5.obj+tmt6.obj+tmt7.obj+   >> tmt_z.lnk
	        echo tmt8.obj+tmt9.obj+tmta.obj+tmti.obj+   >> tmt_z.lnk
	        echo tmtb.obj+tmtc.obj+tmtd.obj+tmte.obj+   >> tmt_z.lnk
	        echo tmtf.obj+tmtg.obj+tmth.obj+            >> tmt_z.lnk
	        echo cholesky.obj+hholder.obj+sort.obj+     >> tmt_z.lnk
	        echo submat.obj+jacobi.obj+newmatrm.obj+    >> tmt_z.lnk
	        echo fft.obj+evalue.obj+bandmat.obj+        >> tmt_z.lnk
	        echo newmatex.obj+except.obj                >> tmt_z.lnk
	        echo tmt_z.exe                              >> tmt_z.lnk

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

tmt.obj:        newmatxx newmatap.h tmt.cxx 
	        $C

tmt1.obj:       newmatxx newmatap.h tmt1.cxx 
	        $C

tmt2.obj:       newmatxx newmatap.h tmt2.cxx 
	        $C

tmt3.obj:       newmatxx newmatap.h tmt3.cxx 
	        $C

tmt4.obj:       newmatxx newmatap.h tmt4.cxx 
	        $C

tmt5.obj:       newmatxx newmatap.h tmt5.cxx 
	        $C

tmt6.obj:       newmatxx newmatap.h tmt6.cxx 
	        $C

tmt7.obj:       newmatxx newmatap.h tmt7.cxx 
	        $C

tmt8.obj:       newmatxx newmatap.h tmt8.cxx 
	        $C

tmt9.obj:       newmatxx newmatap.h tmt9.cxx 
	        $C

tmta.obj:       newmatxx newmatap.h tmta.cxx 
	        $C

tmtb.obj:       newmatxx newmatap.h tmtb.cxx 
	        $C

tmtc.obj:       newmatxx newmatap.h tmtc.cxx 
	        $C

tmtd.obj:       newmatxx newmatap.h tmtd.cxx 
	        $C

tmte.obj:       newmatxx newmatap.h tmte.cxx 
	        $C

tmtf.obj:       newmatxx newmatap.h tmtf.cxx 
	        $C

tmtg.obj:       newmatxx newmatap.h tmtg.cxx 
	        $C

tmth.obj:       newmatxx newmatap.h tmth.cxx
	        $C

tmti.obj:       newmatxx newmatap.h tmti.cxx
	        $C


