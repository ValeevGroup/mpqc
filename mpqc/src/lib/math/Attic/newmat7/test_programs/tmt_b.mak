.AUTODEPEND

#		*Translator Definitions*
CC = bcc +TMT_B.CFG
TASM = TASM
TLIB = tlib
TLINK = tlink
LIBPATH = C:\BORLANDC\LIB
INCLUDEPATH = C:\BORLANDC\INCLUDE


#		*Implicit Rules*
.c.obj:
  $(CC) -c {$< }

.cpp.obj:
  $(CC) -c {$< }

#		*List Macros*


EXE_dependencies =  \
 tmti.obj \
 except.obj \
 newmatex.obj \
 tmth.obj \
 bandmat.obj \
 fft.obj \
 newmat9.obj \
 evalue.obj \
 submat.obj \
 cholesky.obj \
 hholder.obj \
 sort.obj \
 newmatrm.obj \
 jacobi.obj \
 tmtf.obj \
 svd.obj \
 tmte.obj \
 tmtd.obj \
 newmat8.obj \
 tmtc.obj \
 tmtb.obj \
 newmat7.obj \
 newmat6.obj \
 newmat5.obj \
 newmat3.obj \
 newmat4.obj \
 newmat2.obj \
 newmat1.obj \
 tmt.obj \
 tmt1.obj \
 tmt2.obj \
 tmt3.obj \
 tmt4.obj \
 tmt5.obj \
 tmt6.obj \
 tmt7.obj \
 tmt8.obj \
 tmt9.obj \
 tmta.obj \
 tmtg.obj

#		*Explicit Rules*
tmt_b.exe: tmt_b.cfg $(EXE_dependencies)
  $(TLINK) /v/x/c/d/P-/L$(LIBPATH) @&&|
c0l.obj+
tmti.obj+
except.obj+
newmatex.obj+
tmth.obj+
bandmat.obj+
fft.obj+
newmat9.obj+
evalue.obj+
submat.obj+
cholesky.obj+
hholder.obj+
sort.obj+
newmatrm.obj+
jacobi.obj+
tmtf.obj+
svd.obj+
tmte.obj+
tmtd.obj+
newmat8.obj+
tmtc.obj+
tmtb.obj+
newmat7.obj+
newmat6.obj+
newmat5.obj+
newmat3.obj+
newmat4.obj+
newmat2.obj+
newmat1.obj+
tmt.obj+
tmt1.obj+
tmt2.obj+
tmt3.obj+
tmt4.obj+
tmt5.obj+
tmt6.obj+
tmt7.obj+
tmt8.obj+
tmt9.obj+
tmta.obj+
tmtg.obj
tmt_b
		# no map file
fp87.lib+
mathl.lib+
cl.lib
|


#		*Individual File Dependencies*
tmti.obj: tmt_b.cfg tmti.cxx 
	$(CC) -c tmti.cxx

except.obj: tmt_b.cfg except.cxx 
	$(CC) -c except.cxx

newmatex.obj: tmt_b.cfg newmatex.cxx 
	$(CC) -c newmatex.cxx

tmth.obj: tmt_b.cfg tmth.cxx 
	$(CC) -c tmth.cxx

bandmat.obj: tmt_b.cfg bandmat.cxx 
	$(CC) -c bandmat.cxx

fft.obj: tmt_b.cfg fft.cxx 
	$(CC) -c fft.cxx

newmat9.obj: tmt_b.cfg newmat9.cxx 
	$(CC) -c newmat9.cxx

evalue.obj: tmt_b.cfg evalue.cxx 
	$(CC) -c evalue.cxx

submat.obj: tmt_b.cfg submat.cxx 
	$(CC) -c submat.cxx

cholesky.obj: tmt_b.cfg cholesky.cxx 
	$(CC) -c cholesky.cxx

hholder.obj: tmt_b.cfg hholder.cxx 
	$(CC) -c hholder.cxx

sort.obj: tmt_b.cfg sort.cxx 
	$(CC) -c sort.cxx

newmatrm.obj: tmt_b.cfg newmatrm.cxx 
	$(CC) -c newmatrm.cxx

jacobi.obj: tmt_b.cfg jacobi.cxx 
	$(CC) -c jacobi.cxx

tmtf.obj: tmt_b.cfg tmtf.cxx 
	$(CC) -c tmtf.cxx

svd.obj: tmt_b.cfg svd.cxx 
	$(CC) -c svd.cxx

tmte.obj: tmt_b.cfg tmte.cxx 
	$(CC) -c tmte.cxx

tmtd.obj: tmt_b.cfg tmtd.cxx 
	$(CC) -c tmtd.cxx

newmat8.obj: tmt_b.cfg newmat8.cxx 
	$(CC) -c newmat8.cxx

tmtc.obj: tmt_b.cfg tmtc.cxx 
	$(CC) -c tmtc.cxx

tmtb.obj: tmt_b.cfg tmtb.cxx 
	$(CC) -c tmtb.cxx

newmat7.obj: tmt_b.cfg newmat7.cxx 
	$(CC) -c newmat7.cxx

newmat6.obj: tmt_b.cfg newmat6.cxx 
	$(CC) -c newmat6.cxx

newmat5.obj: tmt_b.cfg newmat5.cxx 
	$(CC) -c newmat5.cxx

newmat3.obj: tmt_b.cfg newmat3.cxx 
	$(CC) -c newmat3.cxx

newmat4.obj: tmt_b.cfg newmat4.cxx 
	$(CC) -c newmat4.cxx

newmat2.obj: tmt_b.cfg newmat2.cxx 
	$(CC) -c newmat2.cxx

newmat1.obj: tmt_b.cfg newmat1.cxx 
	$(CC) -c newmat1.cxx

tmt.obj: tmt_b.cfg tmt.cxx 
	$(CC) -c tmt.cxx

tmt1.obj: tmt_b.cfg tmt1.cxx 
	$(CC) -c tmt1.cxx

tmt2.obj: tmt_b.cfg tmt2.cxx 
	$(CC) -c tmt2.cxx

tmt3.obj: tmt_b.cfg tmt3.cxx 
	$(CC) -c tmt3.cxx

tmt4.obj: tmt_b.cfg tmt4.cxx 
	$(CC) -c tmt4.cxx

tmt5.obj: tmt_b.cfg tmt5.cxx 
	$(CC) -c tmt5.cxx

tmt6.obj: tmt_b.cfg tmt6.cxx 
	$(CC) -c tmt6.cxx

tmt7.obj: tmt_b.cfg tmt7.cxx 
	$(CC) -c tmt7.cxx

tmt8.obj: tmt_b.cfg tmt8.cxx 
	$(CC) -c tmt8.cxx

tmt9.obj: tmt_b.cfg tmt9.cxx 
	$(CC) -c tmt9.cxx

tmta.obj: tmt_b.cfg tmta.cxx 
	$(CC) -c tmta.cxx

tmtg.obj: tmt_b.cfg tmtg.cxx 
	$(CC) -c tmtg.cxx

#		*Compiler Configuration File*
tmt_b.cfg: tmt_b.mak
  copy &&|
-ml
-3
-f287
-N
-H=TMT_B.SYM
-weas
-wpre
-I$(INCLUDEPATH)
-L$(LIBPATH)
-P.C
| tmt_b.cfg


