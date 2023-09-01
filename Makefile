# External Makefile

DIRS := Src
DOCSDIRS := Doc

.phony: all docs clean distclean docsclean

all:
	-for dir in $(DIRS); do (cd $$dir; $(MAKE) all); done
	
docs:	
	#cd $(DOCSDIRS); doxygen ProjectPacs.conf; open doxygen/html/index.html
	-for dir in $(DOCSDIRS); do (cd $$dir; doxygen ProjectPacs.conf; open doxygen/html/index.html); done
	
clean: 
	-for dir in $(DIRS); do (cd $$dir; $(MAKE) clean); done

distclean: 
	-for dir in $(DIRS); do (cd $$dir; $(MAKE) distclean); done
	
docsclean:
	-for dir in $(DOCSDIRS); do (cd $$dir; rm -r doxygen); done


