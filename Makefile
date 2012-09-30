# The location of the mesh library




###############################################################################
# File management.  This is where the source, header, and object files are
# defined

#
# source files
srcfiles 	:= $(wildcard *.C)  
srcfiles 	+= $(wildcard ./src/*.C)
srcfiles 	+= $(wildcard ./src/*.c)

includefiles 	:= $(wildcard ./include/*.h)

#SOURCE := ./src
obj-suffix	:= o
#
# object files
objects		:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles))
objects		:= $(patsubst %.c, %.$(obj-suffix), $(objects))


###############################################################################
CC		= gcc
CXX		= g++
CXXFLAGS	=  -g -DSELF_CHECK -O3 -Wno-write-strings
INCLUDES	= -Iinclude/ -I/usr/X11R6/include/ 
LIBS 		= -L/usr/X11R6/lib -lX11 -lXi -lXmu -lglut -lGL -lGLU 
#OBJS = $(SRCS:.cxx=.o)

# PREDCXXFLAGS is for compiling J. Shewchuk's predicates. It should
# always be equal to -O0 (no optimization). Otherwise, TetGen may not
# work properly.



.PHONY: clean clobber distclean

###############################################################################
# Target:
#
target 	   := ./pr


all:: $(target)

# Production rules:  how to make the target - depends on library configuration
$(target): $(objects) $(includefiles)
#$(CXX)  $(INCLUDES) $(CXXFLAGS) $(objects) -c $<
	@echo "Linking "$@"..."
	$(CXX) $(CXXFLAGS)  $(LIBS) $(objects) -o $@ 
# @$(libmesh_CXX) $(libmesh_CPPFLAGS) $(libmesh_CXXFLAGS) $(objects) -o $@ $(libmesh_LIBS) $(libmesh_LDFLAGS)

%.o: %.C
	$(CXX)  $(INCLUDES) $(CXXFLAGS) -c $< -o $@

%.o: %.c
	$(CC)  $(INCLUDES) $(CXXFLAGS) -c $< -o $@
# Useful rules.
clean:
	@rm -f $(objects) *~ .depend
	rm *.vtk

distclean:
	@$(MAKE) clobber
	@rm -f *.o *.g.o *.pg.o .depend

run: $(target)
	@echo "***************************************************************"
	@echo "* Running Example " $(LIBMESH_RUN) $(target) $(LIBMESH_OPTIONS)
	@echo "***************************************************************"
	@echo " "
#@$(LIBMESH_RUN) $(target) $(LIBMESH_OPTIONS)
	@$(target) mouse.vtk
	@echo " "
	@echo "***************************************************************"
	@echo "* Done Running Example " $(LIBMESH_RUN) $(target) $(LIBMESH_OPTIONS)
	@echo "***************************************************************"

edit:
	emacs $(srcfiles)  $(includefiles)&
trolo:
	@echo -e "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n"

gnuplot:
	gnuplot -e "plot 'lolfile.txt' using 1:2 w l" -p
