ORIGIN = PWB
ORIGIN_VER = 2.0
PROJ = EX_MS
PROJFILE = EX_MS.MAK
DEBUG = 0

CC  = cl
CFLAGS_G  = /W2 /BATCH
CFLAGS_D  = /f /Zi /Od
CFLAGS_R  = /f- /Ot /Oi /Ol /Oe /Og /Gs
CXX  = cl
CXXFLAGS_G  = /AL /W2 /BATCH
CXXFLAGS_D  = /f /Od /Zi
CXXFLAGS_R  = /f- /Ot /Gs
MAPFILE_D  = NUL
MAPFILE_R  = NUL
LFLAGS_G  = /NOI /BATCH /ONERROR:NOEXE
LFLAGS_D  = /CO /FAR /PACKC
LFLAGS_R  = /EXE /FAR /PACKC
LINKER	= link
ILINK  = ilink
LRF  = echo > NUL
ILFLAGS  = /a /e

FILES  = BANDMAT.CXX CHOLESKY.CXX EXAMPLE.CXX EXCEPT.CXX HHOLDER.CXX\
	NEWMAT1.CXX NEWMAT2.CXX NEWMAT3.CXX NEWMAT4.CXX NEWMAT5.CXX\
	NEWMAT6.CXX NEWMAT7.CXX NEWMAT8.CXX NEWMATEX.CXX SUBMAT.CXX SVD.CXX\
	NEWMATRM.CXX
OBJS  = BANDMAT.obj CHOLESKY.obj EXAMPLE.obj EXCEPT.obj HHOLDER.obj\
	NEWMAT1.obj NEWMAT2.obj NEWMAT3.obj NEWMAT4.obj NEWMAT5.obj\
	NEWMAT6.obj NEWMAT7.obj NEWMAT8.obj NEWMATEX.obj SUBMAT.obj SVD.obj\
	NEWMATRM.obj

all: $(PROJ).exe

.SUFFIXES:
.SUFFIXES:
.SUFFIXES: .obj .cxx

BANDMAT.obj : BANDMAT.CXX include.h newmat.h newmatrc.h boolean.h except.h\
	controlw.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /FoBANDMAT.obj BANDMAT.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /FoBANDMAT.obj BANDMAT.CXX
<<
!ENDIF

CHOLESKY.obj : CHOLESKY.CXX include.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /FoCHOLESKY.obj CHOLESKY.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /FoCHOLESKY.obj CHOLESKY.CXX
<<
!ENDIF

EXAMPLE.obj : EXAMPLE.CXX include.h newmatap.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /FoEXAMPLE.obj EXAMPLE.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /FoEXAMPLE.obj EXAMPLE.CXX
<<
!ENDIF

EXCEPT.obj : EXCEPT.CXX include.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /FoEXCEPT.obj EXCEPT.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /FoEXCEPT.obj EXCEPT.CXX
<<
!ENDIF

HHOLDER.obj : HHOLDER.CXX include.h newmatap.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /FoHHOLDER.obj HHOLDER.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /FoHHOLDER.obj HHOLDER.CXX
<<
!ENDIF

NEWMAT1.obj : NEWMAT1.CXX include.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /FoNEWMAT1.obj NEWMAT1.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /FoNEWMAT1.obj NEWMAT1.CXX
<<
!ENDIF

NEWMAT2.obj : NEWMAT2.CXX include.h newmat.h newmatrc.h boolean.h except.h\
	controlw.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /FoNEWMAT2.obj NEWMAT2.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /FoNEWMAT2.obj NEWMAT2.CXX
<<
!ENDIF

NEWMAT3.obj : NEWMAT3.CXX include.h newmat.h newmatrc.h boolean.h except.h\
	controlw.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /FoNEWMAT3.obj NEWMAT3.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /FoNEWMAT3.obj NEWMAT3.CXX
<<
!ENDIF

NEWMAT4.obj : NEWMAT4.CXX include.h newmat.h newmatrc.h boolean.h except.h\
	controlw.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /FoNEWMAT4.obj NEWMAT4.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /FoNEWMAT4.obj NEWMAT4.CXX
<<
!ENDIF

NEWMAT5.obj : NEWMAT5.CXX include.h newmat.h newmatrc.h boolean.h except.h\
	controlw.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /FoNEWMAT5.obj NEWMAT5.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /FoNEWMAT5.obj NEWMAT5.CXX
<<
!ENDIF

NEWMAT6.obj : NEWMAT6.CXX include.h newmat.h newmatrc.h boolean.h except.h\
	controlw.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /FoNEWMAT6.obj NEWMAT6.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /FoNEWMAT6.obj NEWMAT6.CXX
<<
!ENDIF

NEWMAT7.obj : NEWMAT7.CXX include.h newmat.h newmatrc.h boolean.h except.h\
	controlw.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /FoNEWMAT7.obj NEWMAT7.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /FoNEWMAT7.obj NEWMAT7.CXX
<<
!ENDIF

NEWMAT8.obj : NEWMAT8.CXX include.h newmatap.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /FoNEWMAT8.obj NEWMAT8.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /FoNEWMAT8.obj NEWMAT8.CXX
<<
!ENDIF

NEWMATEX.obj : NEWMATEX.CXX include.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /FoNEWMATEX.obj NEWMATEX.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /FoNEWMATEX.obj NEWMATEX.CXX
<<
!ENDIF

SUBMAT.obj : SUBMAT.CXX include.h newmat.h newmatrc.h boolean.h except.h\
	controlw.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /FoSUBMAT.obj SUBMAT.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /FoSUBMAT.obj SUBMAT.CXX
<<
!ENDIF

SVD.obj : SVD.CXX include.h newmat.h newmatrm.h precisio.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /FoSVD.obj SVD.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /FoSVD.obj SVD.CXX
<<
!ENDIF

NEWMATRM.obj : NEWMATRM.CXX include.h newmat.h newmatrm.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /FoNEWMATRM.obj NEWMATRM.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /FoNEWMATRM.obj NEWMATRM.CXX
<<
!ENDIF


$(PROJ).exe : $(OBJS)
!IF $(DEBUG)
	$(LRF) @<<$(PROJ).lrf
$(RT_OBJS: = +^
) $(OBJS: = +^
)
$@
$(MAPFILE_D)
$(LIBS: = +^
) +
$(LLIBS_G: = +^
) +
$(LLIBS_D: = +^
)
$(DEF_FILE) $(LFLAGS_G) $(LFLAGS_D);
<<
!ELSE
	$(LRF) @<<$(PROJ).lrf
$(RT_OBJS: = +^
) $(OBJS: = +^
)
$@
$(MAPFILE_R)
$(LIBS: = +^
) +
$(LLIBS_G: = +^
) +
$(LLIBS_R: = +^
)
$(DEF_FILE) $(LFLAGS_G) $(LFLAGS_R);
<<
!ENDIF
	$(LINKER) @$(PROJ).lrf


.cxx.obj :
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /Fo$@ $<
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /Fo$@ $<
<<
!ENDIF


run: $(PROJ).exe
	$(PROJ).exe $(RUNFLAGS)

debug: $(PROJ).exe
	CV $(CVFLAGS) $(PROJ).exe $(RUNFLAGS)
