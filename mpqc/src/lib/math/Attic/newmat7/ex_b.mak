.AUTODEPEND

#		*Translator Definitions*
CC = bcc +EX_B.CFG
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
 bandmat.obj \
 cholesky.obj \
 evalue.obj \
 example.obj \
 except.obj \
 fft.obj \
 hholder.obj \
 jacobi.obj \
 newmat1.obj \
 newmat2.obj \
 newmat3.obj \
 newmat4.obj \
 newmat5.obj \
 newmat6.obj \
 newmat7.obj \
 newmat8.obj \
 newmat9.obj \
 newmatex.obj \
 newmatrm.obj \
 sort.obj \
 submat.obj \
 svd.obj

#		*Explicit Rules*
ex_b.exe: ex_b.cfg $(EXE_dependencies)
  $(TLINK) /v/x/c/P-/L$(LIBPATH) @&&|
c0l.obj+
bandmat.obj+
cholesky.obj+
evalue.obj+
example.obj+
except.obj+
fft.obj+
hholder.obj+
jacobi.obj+
newmat1.obj+
newmat2.obj+
newmat3.obj+
newmat4.obj+
newmat5.obj+
newmat6.obj+
newmat7.obj+
newmat8.obj+
newmat9.obj+
newmatex.obj+
newmatrm.obj+
sort.obj+
submat.obj+
svd.obj
ex_b
		# no map file
emu.lib+
mathl.lib+
cl.lib
|


#		*Individual File Dependencies*
bandmat.obj: ex_b.cfg bandmat.cxx 
	$(CC) -c bandmat.cxx

cholesky.obj: ex_b.cfg cholesky.cxx 
	$(CC) -c cholesky.cxx

evalue.obj: ex_b.cfg evalue.cxx 
	$(CC) -c evalue.cxx

example.obj: ex_b.cfg example.cxx 
	$(CC) -c example.cxx

except.obj: ex_b.cfg except.cxx 
	$(CC) -c except.cxx

fft.obj: ex_b.cfg fft.cxx 
	$(CC) -c fft.cxx

hholder.obj: ex_b.cfg hholder.cxx 
	$(CC) -c hholder.cxx

jacobi.obj: ex_b.cfg jacobi.cxx 
	$(CC) -c jacobi.cxx

newmat1.obj: ex_b.cfg newmat1.cxx 
	$(CC) -c newmat1.cxx

newmat2.obj: ex_b.cfg newmat2.cxx 
	$(CC) -c newmat2.cxx

newmat3.obj: ex_b.cfg newmat3.cxx 
	$(CC) -c newmat3.cxx

newmat4.obj: ex_b.cfg newmat4.cxx 
	$(CC) -c newmat4.cxx

newmat5.obj: ex_b.cfg newmat5.cxx 
	$(CC) -c newmat5.cxx

newmat6.obj: ex_b.cfg newmat6.cxx 
	$(CC) -c newmat6.cxx

newmat7.obj: ex_b.cfg newmat7.cxx 
	$(CC) -c newmat7.cxx

newmat8.obj: ex_b.cfg newmat8.cxx 
	$(CC) -c newmat8.cxx

newmat9.obj: ex_b.cfg newmat9.cxx 
	$(CC) -c newmat9.cxx

newmatex.obj: ex_b.cfg newmatex.cxx 
	$(CC) -c newmatex.cxx

newmatrm.obj: ex_b.cfg newmatrm.cxx 
	$(CC) -c newmatrm.cxx

sort.obj: ex_b.cfg sort.cxx 
	$(CC) -c sort.cxx

submat.obj: ex_b.cfg submat.cxx 
	$(CC) -c submat.cxx

svd.obj: ex_b.cfg svd.cxx 
	$(CC) -c svd.cxx

#		*Compiler Configuration File*
ex_b.cfg: ex_b.mak
  copy &&|
-ml
-wpro
-weas
-wpre
-I$(INCLUDEPATH)
-L$(LIBPATH)
-P
| ex_b.cfg


