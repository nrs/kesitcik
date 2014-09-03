# AUTOMAKE_OPTIONS = subdir-objects
# SUBDIRS = include
# ACLOCAL_AMFLAGS = -I m4

# Tools in the auxillary directory
AUX_DIST = build-aux/install-sh \
	   build-aux/missing

INCLUDE_DIRS = -I$(abs_srcdir)/include/
AM_CFLAGS = -g -Wall -O2 $(INCLUDE_DIRS)
AM_CPPFLAGS = -g -Wall -O2 $(INCLUDE_DIRS)
AM_LDFLAGS = -L/usr/X11R6/lib -lX11 -lXi -lXmu -lglut -lGL -lGLU
# -L/usr/X11R6/lib -lX11 -lm

# this lists the binaries to produce, the (non-PHONY, binary) targets in
# the previous manual Makefile

lib_LTLIBRARIES  = libkesitcik.la

libkesitcik_la_SOURCES = src/constitutive.C  src/interval.C  src/predicates.C  src/triangle.C \
src/glutMaster.C    src/point.C     src/prgl.C        src/type_vector.C \
src/glutWindow.C    src/pr1.C       src/section.C

libkesitcik_la_HEADERS = include/kesitcik/compare_types.h    include/kesitcik/libmesh_exceptions.h include/kesitcik/constitutive.h     include/kesitcik/ostream_proxy.h include/kesitcik/glutMaster.h       include/kesitcik/point.h include/kesitcik/glutWindow.h       include/kesitcik/pr1.h include/kesitcik/interval.h         include/kesitcik/predicates.h include/kesitcik/libmesh_base.h     include/kesitcik/prgl.h include/kesitcik/libmesh_C_isnan.h  include/kesitcik/section.h include/kesitcik/libmesh_common.h   include/kesitcik/triangle.h include/kesitcik/libmesh_config.h   include/kesitcik/type_vector.h
# include_HEADERS = libkesitcik_la_HEADERS
libkesitcik_ladir = $(includedir)/kesitcik
