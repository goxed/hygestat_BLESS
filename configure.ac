AC_PREREQ([2.69])
AC_INIT([hygestat_bless], [1.3.3], [befongan@utmb.edu], [hygestat_bless],[utmb.edu])
AC_CONFIG_SRCDIR([src/hygestat_windows.cpp])
AC_PREFIX_DEFAULT(${HOME})

AC_PROG_CXX
AC_PROG_AWK
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_MKDIR_P
AC_CHECK_LIB([boost_filesystem], [main])
AC_CHECK_LIB([boost_system], [main])
AC_CHECK_LIB([gsl], [main])
AC_CHECK_LIB([gslcblas], [main])
AC_CHECK_LIB([pthread], [main])
AM_INIT_AUTOMAKE
AC_LANG([C++])
AC_GNU_SOURCE
LIBS=" -lm -lgsl -lgslcblas -lboost_system -lboost_filesystem -lpthread"
CXXFLAGS='-Wno-comment -Wno-sign-compare -std=c++11 -O3 -fopenmp -flto -march=native'
AC_LANG_PUSH([C++])
AC_CHECK_HEADERS([boost/foreach.hpp])
AC_LANG_POP([C++])
AC_TYPE_SIZE_T

# who is installing
USER=`whoami`

# Check for platform-specific settings
case $host in
   *cygwin*)
      AC_DEFINE(Cygwin, [], [Using cygwin])
      MACHID="cygwin"
      CFLAGS="$CFLAGS -Wl,--stack,9000000 "
      ;;
   *darwin*)
      AC_DEFINE(MACOSX, [], [Using MacOS X])
      MACHID="darwin"
      ;;
   *linux*)
      AC_DEFINE(Linux, [], [Using Linux])
      MACHID="linux"
      ;;
   *aix*)
      ac_mpicc_flags="-O3 -qstrict -qarch=pwr3 -qtune=pwr3"
      AC_DEFINE(AIX, [], [Using AIX])
      MACHID="aix"
      if /usr/sbin/lsattr -EHl proc0 | grep POWER >/dev/null 2>&1; then
          MACHID=rs6000
      fi
      ;;
   *cray*)
      ac_mpicc_flags="-O3 -h msglevel_3"
      AC_DEFINE(cray3e, [], [Using cray3e])
      MACHID="cray"
      ;;
   *irix6*)
      ac_mpicc_flags="-O3"
      AC_DEFINE(IRIX64, [], [Using IRIX64])
      MACHID="irix"
      ;;
   *irix*)
      ac_mpicc_flags="-O3"
      MACHID="irix"
      ;;
   *solaris*|*sunos*)
      if test "${CC}" != "gcc"; then
          ac_cc_opt="-O -Xa"
      fi
      AC_DEFINE(SunOS, [], [Using SunOS])
      AC_CHECK_LIB([nsl], [gethostbyname], [], [AC_MSG_ERROR([Nsl library not found])])
      AC_CHECK_LIB([socket], [socket], [], [AC_MSG_ERROR([Socket library not found])])
      LIBS="-lsocket -lnsl -lm"
      MACHID="solaris"
      ;;
esac


####################################
# Check for configuration switches
####################################

# debugging
AC_ARG_ENABLE(
    debug,
    [  --enable-debug           Produce an executable with debugging symbols. Default is disabled.],
    [if test "$enableval" = "yes"; then
        DEBUG_FLAGS="$ac_cc_debug";
        ac_cc_opt="-O0";
     fi],
    [DEBUG_FLAGS=""])

# optimization
AC_ARG_ENABLE(
    opt,
    [  --enable-opt             Produce an executable with optimization. Default compiler dependent.],
    [if test "$enableval" = "yes"; then
       OPT_FLAGS="$ac_cc_opt";
     fi],
    [OPT_FLAGS="$ac_cc_opt"])

# get mpi info, and check if parallel version can be compiled
AC_ARG_WITH(
    mpidir,
    [  --with-mpidir=MPIDIR    Set MPI installation path. Default is path of found mpicc.],
    [mpidir="$withval";
     MPICC="";
     MPIRUN="";
     AC_PATH_PROGS(MPICC, mpicc openmpicc mpcc, [], [$mpidir/bin])
     AC_PATH_PROGS(MPIRUN, mpirun openmpirun, [], [$mpidir/bin])
     MPILIBDIR=$mpidir/lib;
     MPIINC="-DMPI -I$mpidir/include";
     AC_CHECK_LIB([mpi], [MPI_Comm_rank])
     parallel="yes";],
    [
     AC_PATH_PROGS(MPICC, mpicc openmpicc mpcc)
     AC_PATH_PROGS(MPIRUN, mpirun openmpirun)
     if test "$MPICC" = ""; then
        MPILIBDIR="";
        MPIINC="";
        parallel="no";
     else
        mpidir=`echo $MPICC | sed -e 's/\/bin\/.*//'`
        MPILIBDIR="$mpidir/lib"
        MPIINC="-DMPI -I$mpidir/include"
        AC_CHECK_LIB([mpi], [MPI_Comm_rank])
        parallel="yes";
     fi])
AC_ARG_WITH(bowtie,
          [  --with-bowtie=PATH_TO_BOWTIE  or standart path,  configure bowtie.],
          [
                 AC_CHECK_PROG(BOWTIE,bowtie,yes,no,"$with_bowtie")
            if test "$BOWTIE" == "yes" ; then
                 echo "Creating a symlink for bowtie ...."
                 cd tools/
                 rm -rf bowtie
                 ln -s "$with_bowtie/bowtie" bowtie
                 cd ..

            fi
          ])
            if test "$BOWTIE" != "yes" ; then
            	 AC_CHECK_PROG(BOWTIE_,bowtie,yes)
                 if test "$BOWTIE_" != "yes" ; then
                     AC_MSG_ERROR([Please install bowtie before installing hygestat.])
                 fi
            fi
                 AM_CONDITIONAL([BOWTIE_HOME_SET], [test ! -z "$BOWTIE"])

AC_ARG_WITH(samstat,
          [  --with-samstat=PATH_TO_SAMSTAT   you should have samtools,    use samstat for post-mapping statistics.],
          [
               AC_CHECK_PROG(SAMSTAT,samstat,yes,no,"$with_samstat")
                  if test "$SAMSTAT" != "yes" ; then
               AC_MSG_ERROR([Please install samstat before installing hygestat.])
                  fi
	       AC_CHECK_PROG(SAMTOOLS,samtools,yes,no,thirdPartyCodes/samtools-1.2/)
                  if test "$SAMTOOLS" != "yes" ; then
               AC_MSG_ERROR([Please install samtools before installing hygestat.])
                  fi
                 echo "Creating a symlink for samstat ...."
                 cd tools/
                 rm -rf samstat
                 ln -s "$with_samstat/samstat" samstat
                 cd ..
          ])
               AM_CONDITIONAL([SAMSTAT_HOME_SET], [test ! -z "$SAMSTAT"])
AC_ARG_WITH(fastqc,
          [  --with-fastqc=PATH_TO_FASTQC        use fastqc for sequencing quality analysis.],
          [
              AC_CHECK_PROG(FASTQC,fastqc,yes,no,"$with_fastqc")
              if test "$FASTQC" != "yes" ; then
              AC_MSG_ERROR([Please install fastqc before installing hygestat.])
              fi
              echo "Creating a symlink for fastqc ...."
              cd tools/
              rm -rf fastqc
              ln -s "$with_fastqc/fastqc" fastqc
              cd ..
          ])
              AM_CONDITIONAL([FASTQC_HOME_SET], [test ! -z "$FASTQC"])

AC_ARG_ENABLE
AC_ARG_WITH
AC_PROG_INSTALL
AC_COPYRIGHT("This program belongs to the MAGA's LAB @ UTMB Galveston")
AC_CONFIG_FILES([Makefile
                 include/Makefile
                 src/Makefile
                 thirdPartyCodes/Makefile
                 tools/Makefile])
AM_SILENT_RULES([yes])
AC_OUTPUT
