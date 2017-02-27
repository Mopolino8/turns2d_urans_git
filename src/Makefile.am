noinst_LIBRARIES = libmodules.a liboverturns2d.a
bin_PROGRAMS = overturns2d

libmodules_a_SOURCES = modules.f
liboverturns2d_a_SOURCES = setup.f eulerrhs.f visrhs.f \
                           sources.f bc.f sa_gammatheta.f sst_gammatheta.f samodel2d.f \
                           komegasst.f motion.f tef.f tef_ts.f gettimemetrics_tef_ts.f forces.f objectives.f\
                           objectivef_ts.f

overturns2d_SOURCES = io.f init.f lhs.f time.f step.f connectivity.f Turns_overset.f

overturns2d_FFLAGS = -I../ihc/src/
overturns2d_FCFLAGS = -I../ihc/src/

liboverturns2d_a_LIBADD = libmodules.a
overturns2d_LDADD = libmodules.a liboverturns2d.a ../ihc/src/libihc.a 

include $(top_srcdir)/config/rules.mk

EXTRA_DIST =  sadata.h komega.h gammatheta.h
