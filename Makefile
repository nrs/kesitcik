include Make.common

default 	= all

# source files
#srcfiles 	:= $(wildcard *.C)  

KESIT_DIR = ./


#SOURCE := ./src
#
# object files
objects		:= $(patsubst %.C, %.$(obj-suffix), $(kesit_srcfiles))
objects		:= $(patsubst %.c, %.$(obj-suffix), $(objects))


###############################################################################

#OBJS = $(SRCS:.cxx=.o)

# PREDCXXFLAGS is for compiling J. Shewchuk's predicates. It should
# always be equal to -O0 (no optimization). Otherwise, TetGen may not
# work properly.



.PHONY: clean clobber distclean

###############################################################################
# Target:
#
target 	   := $(kesit_library)


all:: $(target)

# Production rules:  how to make the target - depends on library configuration
$(target): $(objects) $(kesit_includefiles)
#$(KESIT_CXX)  $(KESIT_INCLUDES) $(KESIT_CXXFLAGS) $(objects) -c $<
	@echo "Linking "$@"..."
	$(KESIT_CXX) $(LIBS) $(objects) -shared -o $@ 
# @$(libmesh_CXX) $(libmesh_CPPFLAGS) $(libmesh_CXXFLAGS) $(objects) -o $@ $(libmesh_LIBS) $(libmesh_LDFLAGS)

%.o: %.C
	$(KESIT_CXX)  $(KESIT_INCLUDES) $(KESIT_CXXFLAGS) -c $< -o $@

%.o: %.c
	$(KESIT_CC)  $(KESIT_INCLUDES) $(KESIT_CXXFLAGS) -c $< -o $@

# Useful rules.
clean:
	@rm -f $(objects) ./src/*~ *~ .depend $(kesit_library)


clobber:
	@$(MAKE) clean
	@rm -f $(target) out.gmv
	rm *.vtk

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
	emacs $(kesit_srcfiles) $(kesit_includefiles)&
trolo:
	@echo -e "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n"

gnuplot:
	gnuplot -e "plot 'interaction.txt' using 1:2 w l" -p
