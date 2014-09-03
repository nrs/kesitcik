# Sample makefile
# The shared library and its headers which this example is
# for should be able to be found by the linker and the compiler

CC = gcc
CXX = g++
CXXFLAGS = -g -DSELF_CHECK -O3 -Wno-write-strings -fPIC
INCLUDES = -I/usr/X11R6/include/
LIBS = -L/usr/X11R6/lib -lX11 -lXi -lXmu -lglut -lGL -lGLU

#library as a shared object
obj-suffix := o

# source files
srcfiles := $(wildcard *.C)
includefiles := $(wildcard *.h)

objects	:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles))

.PHONY: clean clobber distclean

# Target:
target := $(shell pwd | xargs basename)
# target := `basename "$(CURID)"`

all: $(target)

# Production rules:  how to make the target - depends on library configuration
$(target): $(objects) $(includefiles)
#$(CXX)  $(INCLUDES) $(CXXFLAGS) $(objects) -c $<
	@echo "Linking "$@"..."
	$(CXX) $(CXXFLAGS) $(LIBS) -lkesitcik $(objects) -o $@

%.o: %.C
	$(CXX) $(INCLUDES) $(CXXFLAGS) -c $< -o $@

%.o: %.c
	$(CC) $(INCLUDES) $(CXXFLAGS) -c $< -o $@

# Useful rules.
clean:
	@rm -f $(objects) ./src/*~ *~ .depend $(target)

clobber:
	@$(MAKE) clean
	@rm -f $(target) out.gmv
	@rm *.vtk dck_* -f

distclean:
	@$(MAKE) clobber
	@rm -f *.o *.g.o *.pg.o .depend

run: $(target)
	@echo "***************************************************************"
	@echo "* Running Example " $(target)
	@echo "***************************************************************"
	@echo " "
	@ ./$(target)
	@echo " "
	@echo "***************************************************************"
	@echo "* Done Running Example " $(target)
	@echo "***************************************************************"

edit:
	emacs $(srcfiles) $(includefiles) &
