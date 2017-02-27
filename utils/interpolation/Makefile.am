noinst_LIBRARIES = libutils.a
bin_PROGRAMS = interpolate

libutils_a_SOURCES = input.f90 grid.f connectivity.f interpolators.f90 solution.f90

interpolate_SOURCES = interpolate.f90 

libutils_a_FFLAGS = -I../../ihc/src/
libutils_a_FCFLAGS = -I../../ihc/src/
interpolate_FFLAGS = -I../../ihc/src/
interpolate_FCFLAGS = -I../../ihc/src/

libutils_a_LIBADD = ../../ihc/src/libihc.a 
interpolate_LDADD = libutils.a ../../ihc/src/libihc.a 

include $(top_srcdir)/config/rules.mk

EXTRA_DIST =
