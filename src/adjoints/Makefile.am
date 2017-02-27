noinst_LIBRARIES = libadoverturns2d.a
bin_PROGRAMS = adoverturns2d adsensitivity errorestimate

libadoverturns2d_a_SOURCES = adBuffer.f adStack.c \
                        modules_ad.f bc_bq.f bcextp_bq.f bcout_inf_bq.f bcout_trkl_bq.f \
                        bcperiodic_bq.f bcsym_bq.f bctany_bq.f bcwake_bq.f bcwall_bq.f \
                        force2d_bq.f iflux_bq.f lim_minmod_bq.f lim_vanleer_bq.f \
                        muscld_new_bq.f quad_bq.f rhsup_bq.f \
                        roeflx_bq.f roetrklflx_bq.f tridag_bq.f \
                        uv_bq.f visrhs_bq.f weno5_bq.f weno_bq.f \
                        io_ad.f sources_ad.f time_ad.f lhs_ad.f bc_ad.f \
                        connectivity_ad.f \
			c_turm_bq.f convec_bq.f diffus_bq.f lamvis_bq.f rhslhs_bq.f \
			samodel2d_ad.f source_bq.f \
			turbc_bq.f turbc_extpt_bq.f turbc_inf_bq.f turbc_out_bq.f \
			turbc_periodic_bq.f turbc_sym_bq.f turbc_wake_bq.f turbc_wall_bq.f \
			vortic_bq.f\
			do_interpolations_ad.f do_interpolations_tan.f \
			calc_cf_bo.f calc_cp_bo.f calc_norm_bo.f objectivef_bo.f\
                        objectivef_ts_ad.f

adoverturns2d_SOURCES = init_ad.f  step_ad_part1.f step_ad_part2.f forces_ad.f Turns_adjoints.f
adsensitivity_SOURCES = Sensitivity.f
errorestimate_SOURCES = ../step.f step_ad.f forces_ad.f error_estimate.f

libadoverturns2d_a_FFLAGS = -I../
libadoverturns2d_a_FCFLAGS = -I../
adoverturns2d_FFLAGS = -I../ -I../../ihc/src/
adoverturns2d_FCFLAGS = -I../ -I../../ihc/src/
adsensitivity_FFLAGS = -I../ -I../../ihc/src/
adsensitivity_FCFLAGS = -I../ -I../../ihc/src/
errorestimate_FFLAGS = -I../ -I../../ihc/src/
errorestimate_FCFLAGS = -I../ -I../../ihc/src/

libadoverturns2d_a_LIBADD = ../libmodules.a
adoverturns2d_LDADD = ../libmodules.a libadoverturns2d.a ../liboverturns2d.a ../../ihc/src/libihc.a
adsensitivity_LDADD = ../libmodules.a libadoverturns2d.a ../liboverturns2d.a ../../ihc/src/libihc.a
errorestimate_LDADD = ../libmodules.a libadoverturns2d.a ../liboverturns2d.a ../../ihc/src/libihc.a

include $(top_srcdir)/config/rules.mk

EXTRA_DIST = DIFFSIZES.inc
