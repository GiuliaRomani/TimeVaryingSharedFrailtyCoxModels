
#DIRS := ${shell find * -maxdepth 0 -type d -print}
DIRS := Src

.phony: all clean

all:
	-for dir in $(DIRS); do (cd $$dir; $(MAKE) all); done

clean: 
	-for dir in $(DIRS); do (cd $$dir; $(MAKE) clean); done

distclean: 
	-for dir in $(DIRS); do (cd $$dir; $(MAKE) distclean); done

