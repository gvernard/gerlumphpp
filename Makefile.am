ACLOCAL_AMFLAGS=-I m4

lib_LTLIBRARIES=libgerlumph.la

libgerlumph_la_SOURCES= \
	src/normal_source/fixed_locs.cpp \
	src/normal_source/image.cpp \
	src/normal_source/light_curve.cpp \
	src/normal_source/mpd.cpp \
	src/normal_source/profile.cpp \
	src/normal_source/velocities.cpp



libgerlumph_la_CPPFLAGS = \
	-std=c++11 \
	-fPIC \
	-Iinclude

#libvkl_la_LDFLAGS=-shared -fPIC 
#libvkl_la_LIBADD = -lgfortran -lCCfits -lcfitsio -lgmp -lCGAL
libgerlumph_la_LIBADD = libmagmap.la libcpu.la



noinst_LTLIBRARIES = libcpu.la libmagmap.la

libcpu_la_SOURCES = src/cpu_gpu/cpu_functions.cpp
libcpu_la_CPPFLAGS = \
	-std=c++11 \
	-fPIC \
	-Iinclude
libcpu_la_LIBADD = -lfftw3

libmagmap_la_SOURCES = src/magmap/magnification_map.cpp
libmagmap_la_CPPFLAGS = \
	-std=c++11 \
	-fPIC \
	-Iinclude \
	-DMAP_PATH='"@MAP_PATH@"'

