AC_INIT([vesicle_sim], [1.0], [kayla.sapp@nih.gov])
AC_CONFIG_SRCDIR([src])
AM_INIT_AUTOMAKE

c_compilers="clang gcc"
cxx_compilers="clang g++"

AC_PROG_CC([$c_compilers])
AC_PROG_CXX([$cxx_compilers])
AC_PROG_CPP
AC_CONFIG_FILES([
Makefile
src/Makefile
])

AC_CHECK_LIB(gsl,main,
        [], 
        [AC_MSG_ERROR("Couldn't find or use the GSL libraries.")])

AC_OUTPUT
