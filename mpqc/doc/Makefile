
PERL=perl
DOCXX=doc++

.PHONY: html
html: sc.dxx include.dxx library.dxx binary.dxx
	$(DOCXX) -d html $<

.PHONY: include.dxx
include.dxx:
	$(PERL) ../bin/mkdoc.pl ../src/lib $@

.PHONY: library.dxx
library.dxx:
	$(PERL) ../bin/mkextradoc.pl ../src/lib $@

.PHONY: binary.dxx
binary.dxx:
	$(PERL) ../bin/mkextradoc.pl ../src/bin $@

.PHONY: veryclean
veryclean:
	/bin/rm -rf html include.dxx library.dxx binary.dxx
	/bin/rm -rf *~

.PHONY: clean
clean:
	/bin/rm -rf include.dxx library.dxx binary.dxx
	/bin/rm -rf *~
