AUTOMAKE_OPTIONS = foreign


fits_templatesdir = $(prefix)/fits_templates
initdir = $(prefix)

if USE_HEADAS
HEASOFT_CFLAGS= -I@HEADAS@/include
HEASOFT_LFLAGS= -L@HEADAS@/lib -lcfitsio
else
HEASOFT_CFLAGS= @LHEAS_CFLAGS@
HEASOFT_LFLAGS= @LHEAS_LDFLAGS@ -lcfitsio
endif

srcbase = src/src2

dist_fits_templates_DATA = fits_template/untrigdat.txt

mysrccodebase = \
$(srcbase)/DBStringUtilities.cpp \
$(srcbase)/DBStringUtilities.hh \
$(srcbase)/DataRecorder.hh \
$(srcbase)/DynMatrix.h \
$(srcbase)/EdgeSet.hh \
$(srcbase)/PHA_IO.hh \
$(srcbase)/PHA_IO.cpp \
$(srcbase)/PHA2_IO.hh \
$(srcbase)/PHA2_IO.cpp \
$(srcbase)/PHAElements.hh \
$(srcbase)/PHAStructures.cpp \
$(srcbase)/PHAStructures.hh \
$(srcbase)/TTEventTable.hh \
$(srcbase)/VectorSpan.h \
$(srcbase)/gbm_geometry.cpp \
$(srcbase)/gbm_geometry.hh \
$(srcbase)/spoccExeUtilities.cpp \
$(srcbase)/spoccExeUtilities.h \
$(srcbase)/GeoTransform.h \
$(srcbase)/GeoTransform.cpp \
$(srcbase)/Gbm_PositionIO.h \
$(srcbase)/Gbm_PositionIO.cpp \
$(srcbase)/algorithms_t.h \
$(srcbase)/algorithms_t.cpp \
$(srcbase)/TTE_IO.hh \
$(srcbase)/TTE_IO.cpp

ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

mylinkerflags = $(HEASOFT_LFLAGS)

mymacros = \
	 \
	-DPKG_PREFIX="\"$(prefix)\"" \
	-DPKG_LDFLAGS="\"$(mylinkerflags)\"" \
	-DPKG_CFLAGS="\"-Wall $(ACLOCAL_AMFLAGS) $(HEASOFT_CFLAGS)\"" \
	-DUTC_FMT_DEFAULT="\"%Y-%m-%d %H:%M:%S.nnnnnn\"" \
	-DPHA1TEMPLATE=\"$(fits_templatesdir)/pha1template.pha\" \
	-DPHA2TEMPLATE=\"$(fits_templatesdir)/pha2template.pha\" \
	-DTTETEMPLATE_FIT=\"$(fits_templatesdir)/tte.fit\" \
	-DWHICH_RM="\"@WHICH_RM@\"" \
	-DWHICH_ENV="\"@WHICH_ENV@\"" \
    -DTRIGDATATEMPLATE="\"$(fits_templatesdir)/untrigdat.txt\""
	-DTTE_DEFAULT_PAD=10.0 \
    -DUNTRIG_EXE="\"gbmCont2UnTrig\""

AM_CPPFLAGS= -Wall \
	$(mymacros) \
	$(HEASOFT_CFLAGS) \
	-I$(top_srcdir)/$(srcbase)

EXTRA_DIST = README

	
bin_PROGRAMS = gbmCont2UnTrig
gbmCont2UnTrig_LDADD = $(mylinkerflags)

gbmCont2UnTrig_SOURCES = $(mysrccodebase) src/gbmCont2UnTrig.cpp
