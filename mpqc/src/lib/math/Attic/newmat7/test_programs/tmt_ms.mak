ORIGIN = PWB
ORIGIN_VER = 2.0
PROJ = TMT_MS
PROJFILE = TMT_MS.MAK
DEBUG = 0

CC  = cl
CFLAGS_G  = /AL /G2 /BATCH
CFLAGS_D  = /f /Od /Zi
CFLAGS_R  = /f- /Ot
CXX  = cl
CXXFLAGS_G  = /AL /BATCH
CXXFLAGS_D  = /f /Od /Zi
CXXFLAGS_R  = /f- /Ot
MAPFILE_D  = NUL
MAPFILE_R  = NUL
LFLAGS_G  = /NOI /STACK:4096 /BATCH /ONERROR:NOEXE
LFLAGS_D  = /CO /FAR /PACKC
LFLAGS_R  = /EXE /FAR /PACKC
LINKER	= link
ILINK  = ilink
LRF  = echo > NUL
ILFLAGS  = /a /e
BSCMAKE  = bscmake
SBRPACK  = sbrpack
NMAKEBSC1  = set
NMAKEBSC2  = nmake

FILES  = EXCEPT.CXX BANDMAT.CXX CHOLESKY.CXX EVALUE.CXX FFT.CXX HHOLDER.CXX\
	JACOBI.CXX NEWMAT1.CXX NEWMAT2.CXX NEWMAT3.CXX NEWMAT4.CXX NEWMAT5.CXX\
	NEWMAT6.CXX NEWMAT7.CXX NEWMAT8.CXX NEWMAT9.CXX NEWMATEX.CXX\
	NEWMATRM.CXX SORT.CXX SUBMAT.CXX SVD.CXX TMT.CXX TMT1.CXX TMT2.CXX\
	TMT3.CXX TMT4.CXX TMT5.CXX TMT6.CXX TMT7.CXX TMT8.CXX TMT9.CXX\
	TMTA.CXX TMTB.CXX TMTC.CXX TMTD.CXX TMTE.CXX TMTF.CXX TMTG.CXX\
	TMTH.CXX TMTI.CXX
OBJS  = EXCEPT.obj BANDMAT.obj CHOLESKY.obj EVALUE.obj FFT.obj HHOLDER.obj\
	JACOBI.obj NEWMAT1.obj NEWMAT2.obj NEWMAT3.obj NEWMAT4.obj NEWMAT5.obj\
	NEWMAT6.obj NEWMAT7.obj NEWMAT8.obj NEWMAT9.obj NEWMATEX.obj\
	NEWMATRM.obj SORT.obj SUBMAT.obj SVD.obj TMT.obj TMT1.obj TMT2.obj\
	TMT3.obj TMT4.obj TMT5.obj TMT6.obj TMT7.obj TMT8.obj TMT9.obj\
	TMTA.obj TMTB.obj TMTC.obj TMTD.obj TMTE.obj TMTF.obj TMTG.obj\
	TMTH.obj TMTI.obj
SBRS  = EXCEPT.sbr BANDMAT.sbr CHOLESKY.sbr EVALUE.sbr FFT.sbr HHOLDER.sbr\
	JACOBI.sbr NEWMAT1.sbr NEWMAT2.sbr NEWMAT3.sbr NEWMAT4.sbr NEWMAT5.sbr\
	NEWMAT6.sbr NEWMAT7.sbr NEWMAT8.sbr NEWMAT9.sbr NEWMATEX.sbr\
	NEWMATRM.sbr SORT.sbr SUBMAT.sbr SVD.sbr TMT.sbr TMT1.sbr TMT2.sbr\
	TMT3.sbr TMT4.sbr TMT5.sbr TMT6.sbr TMT7.sbr TMT8.sbr TMT9.sbr\
	TMTA.sbr TMTB.sbr TMTC.sbr TMTD.sbr TMTE.sbr TMTF.sbr TMTG.sbr\
	TMTH.sbr TMTI.sbr

all: $(PROJ).exe

.SUFFIXES:
.SUFFIXES:
.SUFFIXES: .obj .sbr .cxx

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

EXCEPT.sbr : EXCEPT.CXX include.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_D) /FREXCEPT.sbr EXCEPT.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_R) /FREXCEPT.sbr EXCEPT.CXX
<<
!ENDIF

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

BANDMAT.sbr : BANDMAT.CXX include.h newmat.h newmatrc.h boolean.h except.h\
	controlw.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_D) /FRBANDMAT.sbr BANDMAT.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_R) /FRBANDMAT.sbr BANDMAT.CXX
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

CHOLESKY.sbr : CHOLESKY.CXX include.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_D) /FRCHOLESKY.sbr CHOLESKY.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_R) /FRCHOLESKY.sbr CHOLESKY.CXX
<<
!ENDIF

EVALUE.obj : EVALUE.CXX include.h newmat.h newmatrm.h precisio.h boolean.h\
	except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /FoEVALUE.obj EVALUE.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /FoEVALUE.obj EVALUE.CXX
<<
!ENDIF

EVALUE.sbr : EVALUE.CXX include.h newmat.h newmatrm.h precisio.h boolean.h\
	except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_D) /FREVALUE.sbr EVALUE.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_R) /FREVALUE.sbr EVALUE.CXX
<<
!ENDIF

FFT.obj : FFT.CXX include.h newmatap.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /FoFFT.obj FFT.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /FoFFT.obj FFT.CXX
<<
!ENDIF

FFT.sbr : FFT.CXX include.h newmatap.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_D) /FRFFT.sbr FFT.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_R) /FRFFT.sbr FFT.CXX
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

HHOLDER.sbr : HHOLDER.CXX include.h newmatap.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_D) /FRHHOLDER.sbr HHOLDER.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_R) /FRHHOLDER.sbr HHOLDER.CXX
<<
!ENDIF

JACOBI.obj : JACOBI.CXX include.h newmat.h precisio.h newmatrm.h boolean.h\
	except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /FoJACOBI.obj JACOBI.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /FoJACOBI.obj JACOBI.CXX
<<
!ENDIF

JACOBI.sbr : JACOBI.CXX include.h newmat.h precisio.h newmatrm.h boolean.h\
	except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_D) /FRJACOBI.sbr JACOBI.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_R) /FRJACOBI.sbr JACOBI.CXX
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

NEWMAT1.sbr : NEWMAT1.CXX include.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_D) /FRNEWMAT1.sbr NEWMAT1.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_R) /FRNEWMAT1.sbr NEWMAT1.CXX
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

NEWMAT2.sbr : NEWMAT2.CXX include.h newmat.h newmatrc.h boolean.h except.h\
	controlw.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_D) /FRNEWMAT2.sbr NEWMAT2.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_R) /FRNEWMAT2.sbr NEWMAT2.CXX
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

NEWMAT3.sbr : NEWMAT3.CXX include.h newmat.h newmatrc.h boolean.h except.h\
	controlw.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_D) /FRNEWMAT3.sbr NEWMAT3.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_R) /FRNEWMAT3.sbr NEWMAT3.CXX
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

NEWMAT4.sbr : NEWMAT4.CXX include.h newmat.h newmatrc.h boolean.h except.h\
	controlw.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_D) /FRNEWMAT4.sbr NEWMAT4.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_R) /FRNEWMAT4.sbr NEWMAT4.CXX
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

NEWMAT5.sbr : NEWMAT5.CXX include.h newmat.h newmatrc.h boolean.h except.h\
	controlw.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_D) /FRNEWMAT5.sbr NEWMAT5.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_R) /FRNEWMAT5.sbr NEWMAT5.CXX
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

NEWMAT6.sbr : NEWMAT6.CXX include.h newmat.h newmatrc.h boolean.h except.h\
	controlw.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_D) /FRNEWMAT6.sbr NEWMAT6.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_R) /FRNEWMAT6.sbr NEWMAT6.CXX
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

NEWMAT7.sbr : NEWMAT7.CXX include.h newmat.h newmatrc.h boolean.h except.h\
	controlw.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_D) /FRNEWMAT7.sbr NEWMAT7.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_R) /FRNEWMAT7.sbr NEWMAT7.CXX
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

NEWMAT8.sbr : NEWMAT8.CXX include.h newmatap.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_D) /FRNEWMAT8.sbr NEWMAT8.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_R) /FRNEWMAT8.sbr NEWMAT8.CXX
<<
!ENDIF

NEWMAT9.obj : NEWMAT9.CXX include.h newmat.h newmatrc.h newmatio.h boolean.h\
	except.h controlw.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /FoNEWMAT9.obj NEWMAT9.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /FoNEWMAT9.obj NEWMAT9.CXX
<<
!ENDIF

NEWMAT9.sbr : NEWMAT9.CXX include.h newmat.h newmatrc.h newmatio.h boolean.h\
	except.h controlw.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_D) /FRNEWMAT9.sbr NEWMAT9.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_R) /FRNEWMAT9.sbr NEWMAT9.CXX
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

NEWMATEX.sbr : NEWMATEX.CXX include.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_D) /FRNEWMATEX.sbr NEWMATEX.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_R) /FRNEWMATEX.sbr NEWMATEX.CXX
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

NEWMATRM.sbr : NEWMATRM.CXX include.h newmat.h newmatrm.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_D) /FRNEWMATRM.sbr NEWMATRM.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_R) /FRNEWMATRM.sbr NEWMATRM.CXX
<<
!ENDIF

SORT.obj : SORT.CXX include.h newmatap.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /FoSORT.obj SORT.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /FoSORT.obj SORT.CXX
<<
!ENDIF

SORT.sbr : SORT.CXX include.h newmatap.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_D) /FRSORT.sbr SORT.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_R) /FRSORT.sbr SORT.CXX
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

SUBMAT.sbr : SUBMAT.CXX include.h newmat.h newmatrc.h boolean.h except.h\
	controlw.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_D) /FRSUBMAT.sbr SUBMAT.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_R) /FRSUBMAT.sbr SUBMAT.CXX
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

SVD.sbr : SVD.CXX include.h newmat.h newmatrm.h precisio.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_D) /FRSVD.sbr SVD.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_R) /FRSVD.sbr SVD.CXX
<<
!ENDIF

TMT.obj : TMT.CXX include.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /FoTMT.obj TMT.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /FoTMT.obj TMT.CXX
<<
!ENDIF

TMT.sbr : TMT.CXX include.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_D) /FRTMT.sbr TMT.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_R) /FRTMT.sbr TMT.CXX
<<
!ENDIF

TMT1.obj : TMT1.CXX include.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /FoTMT1.obj TMT1.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /FoTMT1.obj TMT1.CXX
<<
!ENDIF

TMT1.sbr : TMT1.CXX include.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_D) /FRTMT1.sbr TMT1.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_R) /FRTMT1.sbr TMT1.CXX
<<
!ENDIF

TMT2.obj : TMT2.CXX include.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /FoTMT2.obj TMT2.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /FoTMT2.obj TMT2.CXX
<<
!ENDIF

TMT2.sbr : TMT2.CXX include.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_D) /FRTMT2.sbr TMT2.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_R) /FRTMT2.sbr TMT2.CXX
<<
!ENDIF

TMT3.obj : TMT3.CXX include.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /FoTMT3.obj TMT3.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /FoTMT3.obj TMT3.CXX
<<
!ENDIF

TMT3.sbr : TMT3.CXX include.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_D) /FRTMT3.sbr TMT3.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_R) /FRTMT3.sbr TMT3.CXX
<<
!ENDIF

TMT4.obj : TMT4.CXX include.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /FoTMT4.obj TMT4.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /FoTMT4.obj TMT4.CXX
<<
!ENDIF

TMT4.sbr : TMT4.CXX include.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_D) /FRTMT4.sbr TMT4.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_R) /FRTMT4.sbr TMT4.CXX
<<
!ENDIF

TMT5.obj : TMT5.CXX include.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /FoTMT5.obj TMT5.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /FoTMT5.obj TMT5.CXX
<<
!ENDIF

TMT5.sbr : TMT5.CXX include.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_D) /FRTMT5.sbr TMT5.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_R) /FRTMT5.sbr TMT5.CXX
<<
!ENDIF

TMT6.obj : TMT6.CXX include.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /FoTMT6.obj TMT6.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /FoTMT6.obj TMT6.CXX
<<
!ENDIF

TMT6.sbr : TMT6.CXX include.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_D) /FRTMT6.sbr TMT6.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_R) /FRTMT6.sbr TMT6.CXX
<<
!ENDIF

TMT7.obj : TMT7.CXX include.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /FoTMT7.obj TMT7.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /FoTMT7.obj TMT7.CXX
<<
!ENDIF

TMT7.sbr : TMT7.CXX include.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_D) /FRTMT7.sbr TMT7.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_R) /FRTMT7.sbr TMT7.CXX
<<
!ENDIF

TMT8.obj : TMT8.CXX include.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /FoTMT8.obj TMT8.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /FoTMT8.obj TMT8.CXX
<<
!ENDIF

TMT8.sbr : TMT8.CXX include.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_D) /FRTMT8.sbr TMT8.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_R) /FRTMT8.sbr TMT8.CXX
<<
!ENDIF

TMT9.obj : TMT9.CXX include.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /FoTMT9.obj TMT9.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /FoTMT9.obj TMT9.CXX
<<
!ENDIF

TMT9.sbr : TMT9.CXX include.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_D) /FRTMT9.sbr TMT9.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_R) /FRTMT9.sbr TMT9.CXX
<<
!ENDIF

TMTA.obj : TMTA.CXX include.h newmatap.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /FoTMTA.obj TMTA.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /FoTMTA.obj TMTA.CXX
<<
!ENDIF

TMTA.sbr : TMTA.CXX include.h newmatap.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_D) /FRTMTA.sbr TMTA.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_R) /FRTMTA.sbr TMTA.CXX
<<
!ENDIF

TMTB.obj : TMTB.CXX include.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /FoTMTB.obj TMTB.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /FoTMTB.obj TMTB.CXX
<<
!ENDIF

TMTB.sbr : TMTB.CXX include.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_D) /FRTMTB.sbr TMTB.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_R) /FRTMTB.sbr TMTB.CXX
<<
!ENDIF

TMTC.obj : TMTC.CXX include.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /FoTMTC.obj TMTC.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /FoTMTC.obj TMTC.CXX
<<
!ENDIF

TMTC.sbr : TMTC.CXX include.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_D) /FRTMTC.sbr TMTC.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_R) /FRTMTC.sbr TMTC.CXX
<<
!ENDIF

TMTD.obj : TMTD.CXX include.h newmatap.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /FoTMTD.obj TMTD.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /FoTMTD.obj TMTD.CXX
<<
!ENDIF

TMTD.sbr : TMTD.CXX include.h newmatap.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_D) /FRTMTD.sbr TMTD.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_R) /FRTMTD.sbr TMTD.CXX
<<
!ENDIF

TMTE.obj : TMTE.CXX include.h newmatap.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /FoTMTE.obj TMTE.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /FoTMTE.obj TMTE.CXX
<<
!ENDIF

TMTE.sbr : TMTE.CXX include.h newmatap.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_D) /FRTMTE.sbr TMTE.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_R) /FRTMTE.sbr TMTE.CXX
<<
!ENDIF

TMTF.obj : TMTF.CXX include.h newmatap.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /FoTMTF.obj TMTF.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /FoTMTF.obj TMTF.CXX
<<
!ENDIF

TMTF.sbr : TMTF.CXX include.h newmatap.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_D) /FRTMTF.sbr TMTF.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_R) /FRTMTF.sbr TMTF.CXX
<<
!ENDIF

TMTG.obj : TMTG.CXX include.h newmatap.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /FoTMTG.obj TMTG.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /FoTMTG.obj TMTG.CXX
<<
!ENDIF

TMTG.sbr : TMTG.CXX include.h newmatap.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_D) /FRTMTG.sbr TMTG.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_R) /FRTMTG.sbr TMTG.CXX
<<
!ENDIF

TMTH.obj : TMTH.CXX include.h newmatap.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /FoTMTH.obj TMTH.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /FoTMTH.obj TMTH.CXX
<<
!ENDIF

TMTH.sbr : TMTH.CXX include.h newmatap.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_D) /FRTMTH.sbr TMTH.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_R) /FRTMTH.sbr TMTH.CXX
<<
!ENDIF

TMTI.obj : TMTI.CXX include.h newmatap.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_D) /FoTMTI.obj TMTI.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/c $(CXXFLAGS_G)
$(CXXFLAGS_R) /FoTMTI.obj TMTI.CXX
<<
!ENDIF

TMTI.sbr : TMTI.CXX include.h newmatap.h newmat.h boolean.h except.h
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_D) /FRTMTI.sbr TMTI.CXX
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_R) /FRTMTI.sbr TMTI.CXX
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

$(PROJ).bsc : $(SBRS)
	$(BSCMAKE) @<<
$(BRFLAGS) $(SBRS)
<<


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

.cxx.sbr :
!IF $(DEBUG)
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_D) /FR$@ $<
<<
!ELSE
	@$(CXX) @<<$(PROJ).rsp
/Zs $(CXXFLAGS_G)
$(CXXFLAGS_R) /FR$@ $<
<<
!ENDIF


run: $(PROJ).exe
	$(PROJ).exe $(RUNFLAGS)

debug: $(PROJ).exe
	CV $(CVFLAGS) $(PROJ).exe $(RUNFLAGS)
