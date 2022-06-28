ACLOCAL_AMFLAGS=-I m4

lib_LTLIBRARIES=libgerlumph.la

libgerlumph_la_SOURCES= \
	src/normal_source/fixed_locs.cpp \
	src/normal_source/image.cpp \
	src/normal_source/light_curve.cpp \
	src/normal_source/mpd.cpp \
	src/normal_source/profile.cpp \
	src/normal_source/fitsInterface.cpp \
	src/normal_source/rectGrid.cpp \
	src/normal_source/velocities.cpp

libgerlumph_la_CPPFLAGS = -std=c++11 -fPIC -g -frounding-math

include_HEADERS = \
	include/fixed_locs.hpp \
	include/gerlumph.hpp \
	include/image.hpp \
	include/light_curve.hpp \
	include/magnification_map.hpp \
	include/mpd.hpp \
	include/profile.hpp \
	include/fitsInterface.hpp \
	include/rectGrid.hpp \
	include/velocities.hpp



if WITH_GPU

libgerlumph_la_LIBADD = libmagmap.la libgpu.la
noinst_LTLIBRARIES = libmagmap.la libgpu.la

# GPU functions
.cu.o:
	$(NVCC) $(CUDA_FLAGS) -o $@ -c $<
.cu.lo:
	python $(top_srcdir)/cudalt.py $@ $(NVCC) \
	--compiler-options=\"$(CFLAGS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) \" -c $<

LINK = $(LIBTOOL) --mode=link $(CC) -o $@ $(CUDA_LIBS)
libgpu_la_SOURCES = src/cpu_gpu/gpu_functions.cu
libgpu_la_LIBADD  = $(CUDA_LIBS)

else

libgerlumph_la_LIBADD = libmagmap.la libcpu.la
noinst_LTLIBRARIES = libmagmap.la libcpu.la

# CPU functions
libcpu_la_SOURCES = src/cpu_gpu/cpu_functions.cpp
libcpu_la_CXXFLAGS = -std=c++11	-fPIC -g -frounding-math
libcpu_la_LIBADD = -lfftw3

endif


# MAP PATH
libmagmap_la_SOURCES = src/magmap/magnification_map.cpp
libmagmap_la_CPPFLAGS = -std=c++11 -fPIC -g -frounding-math -DMAP_PATH='"@MAP_PATH@"'
