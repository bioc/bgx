#  This file is part of BGX, the Bayesian Gene eXpression program.
#  Copyright 2007 Ernest Turro <ernest.turro@ic.ac.uk>
#
#  BGX is free software; you can redistribute it and/or modify it
#  under the terms of the GNU General Public License, version 2, as
#  published by the Free Software Foundation.
#
#  BGX is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


# Initialise.
AC_INIT([bgx], [1.0.4], [ernest.turro@ic.ac.uk])
AC_CONFIG_SRCDIR([src/bgx.cc])
#AC_CONFIG_HEADERS([src/config.h])

# Checks for programs.
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CXX=`"${R_HOME}/bin/R" CMD config CXX`

# Set default flags
#CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`
#CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
#CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`

# Add optimisation flags.
#AX_CC_MAXOPT
# Catering for gnu only
AX_CHECK_COMPILER_FLAGS(-ffast-math,CFLAGS="$CFLAGS -ffast-math")
AX_CHECK_COMPILER_FLAGS(-fomit-frame-pointer,CFLAGS="$CFLAGS -fomit-frame-pointer")
AX_CHECK_COMPILER_FLAGS(-O3,CFLAGS="$CFLAGS -O3")
AX_GCC_ARCHFLAG([no])
BGX_CFLAGS=$CFLAGS
BGX_CXXFLAGS=$CFLAGS

# We are compiling for R.
BGX_CXXFLAGS="$BGX_CXXFLAGS -DUSING_R"

# Get C99 flag(s) (if any), and append to BGX_CFLAGS.
PRIORCC=$CC
AC_PROG_CC_C99 # May change CC
for i in $CC; do 
  if test `expr "$PRIORCC" : ".*$i"` -eq 0;then 
    BGX_CFLAGS="$BGX_CFLAGS $i"
  fi;
done

#Link against math library (sets LIBS)
AC_CHECK_LIB(m, main)

# Checks for Boost libraries.
AX_BOOST([1.31.0], BGX_CPPFLAGS="$BOOST_CPPFLAGS")

# Checks for typedefs, structures, and compiler characteristics.
#AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE

# Finish.
AC_SUBST(BGX_CPPFLAGS)
AC_SUBST(BGX_CXXFLAGS)
AC_SUBST(BGX_CFLAGS)
AC_SUBST(LIBS)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
