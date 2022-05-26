AC_INIT([libgerlumph], [0.1], [gvernard@astro.rug.nl])
AC_CONFIG_AUX_DIR(aux-dist)
AM_INIT_AUTOMAKE([foreign subdir-objects])
AC_CONFIG_MACRO_DIRS([m4])

AC_PROG_CXX


LT_INIT([disable-static])
AC_FUNC_MALLOC
AC_CHECK_FUNCS([floor pow sqrt])




##########################################################################
LIBS="-lCCfits -lcfitsio -lpng $LIBS"

# introduce the optional configure parameter for a non-standard install prefix of XXX
AC_ARG_WITH([external],
    [AS_HELP_STRING([--with-external=prefix],[try this for a non-standard install prefix of all the external libraries])],
    [EXT_PATHSET=1],
    [EXT_PATHSET=0])
	
# if optional parameter used, extend path flags for compliler and linker
if test $EXT_PATHSET = 1 ; then
    # extend the compiler and linker flags according to the path set
    AM_CXXFLAGS="-I${with_external}cfitsio/include -I${with_external}CCfits/include -I${with_external}libpng/include -I${with_external}fftw/include $AM_CXXFLAGS"
    AM_LDFLAGS="$AM_LDFLAGS -L${with_external}cfitsio/lib -L${with_external}CCfits/lib -L${with_external}libpng/lib -L${with_external}fftw/lib"
fi
##########################################################################


##########################################################################
AC_ARG_WITH([map-path],
    [AS_HELP_STRING([--with-map-path=prefix],[give the path to the gerlumph maps])],
    [MAP_PATHSET=1],
    [MAP_PATHSET=0])
	
# if optional parameter used, extend path flags for compliler and linker
if test $MAP_PATHSET = 1 ; then
    MAP_PATH="$with_map_path"
else
    AC_MSG_ERROR([You need to provide the MAP_PATH to comile the gerlumph++ library!])
fi
##########################################################################


##########################################################################
AC_ARG_ENABLE([gpu],
	[AS_HELP_STRING([--enable-gpu],[replace convolutions on the CPU (fftw) by the GPU (cufft)])],
	[GPU_SET=1],
	[GPU_SET=0])
AM_CONDITIONAL([WITH_GPU], [test $GPU_SET = 1])
    
CUDA_CFLAGS="-std=c++11 --compiler-options '-fPIC' -Wno-deprecated-gpu-targets"
CUDA_LIBS="-lcudart -lcufft"
NVCC="nvcc"

AC_SUBST(CUDA_CFLAGS)
AC_SUBST(CUDA_LIBS)
AC_SUBST(NVCC)
##########################################################################



# distribute the changed variables among the Makefiles
AC_SUBST([MAP_PATH])
AC_SUBST([LIBS])
AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_LDFLAGS])



AC_CONFIG_FILES([Makefile])

AC_OUTPUT
