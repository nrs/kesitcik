# The location of the mesh library
KESIT_DIR = ../../

include $(KESIT_DIR)/Make.common


# source files
srcfiles 	:= $(wildcard *.C)  

objects		:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles))


###############################################################################

.PHONY: clean clobber distclean

###############################################################################
# Target:
#
target 	   := ex6


all:: $(target)

# Production rules:  how to make the target - depends on library configuration
$(target): $(objects) $(includefiles)
#$(CXX)  $(INCLUDES) $(CXXFLAGS) $(objects) -c $<
	@echo "Linking "$@"..."
	$(KESIT_CXX) $(KESIT_CXXFLAGS)  $(KESIT_LIBS) $(kesit_library) $(objects) -o $@ 


%.o: %.C
	$(KESIT_CXX)  $(KESIT_INCLUDES) $(KESIT_CXXFLAGS) -c $< -o $@

%.o: %.c
	$(KESIT_CC)  $(KESIT_INCLUDES) $(KESIT_CXXFLAGS) -c $< -o $@
# Useful rules.
clean:
	@rm -f $(objects) ./src/*~ *~ .depend $(target)


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
	emacs $(local_srcfiles)  $(local_includefiles)&

trolo:
	@echo -e "\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n\n"

plot:
	gnuplot plot.gp -p
