# AC_CHECK_IO_LIB([DEFAULT-ACTION], [MINIMUM-VERSION])
# Autoconf macro to find io_lib.
# If found it defines HAVE_IO_LIB and sets IO_LIB_CPPFLAGS and IO_LIB_LDFLAGS.
#
# DEFAULT-ACTION is the string "yes" or "no", defaulting to "yes" when not
# specified.
# Requires io_lib 1.10.0 or above (for the io_lib-config script).

AC_DEFUN([AC_CHECK_IO_LIB],
[
  AC_ARG_WITH(io_lib, AC_HELP_STRING([--with-io_lib=DIR],
	                             [Look for io_lib root in DIR]),
	      [_io_lib_with=$withval],
	      [_io_lib_with=ifelse([$1],,[yes],[$1])])

  # Defaults to enabled
  if test "$_io_lib_with" != "no"
  then
    # Identify the location of io_lib-config
    if test -d "$_io_lib_with"
    then
      if test -x "$_io_lib_with/bin/io_lib-config"
      then
        _io_lib_config="$_io_lib_with/bin/io_lib-config"
      else
        _io_lib_config=
      fi
    else
      AC_PATH_PROG([_io_lib_config], [io_lib-config])
    fi

    # Check version is sufficient; sneakily entirely in sh syntax
    if test x$_io_lib_config != "x"
    then
      _io_lib_version=`$_io_lib_config --version`
      SAVE_IFS=$IFS; IFS=.
      _val=0
      for v in $_io_lib_version; do _val=`expr $_val '*' 100 + $v`; done
      _io_lib_version=$_val
      IFS=$SAVE_IFS

      _io_lib_wanted=`echo ifelse([$2],,[0],[$2])`
      SAVE_IFS=$IFS; IFS=.
      _val=0
      for v in $_io_lib_wanted; do _val=`expr $_val '*' 100 + $v`; done
      _io_lib_wanted=$_val
      IFS=$SAVE_IFS

      if test $_io_lib_version -ge $_io_lib_wanted
      then
        io_lib_version_ok=yes
      else
        io_lib_version_ok=no
      fi
    else
      io_lib_version_ok=yes; # Have to just guess and hope
    fi

    if test $io_lib_version_ok = "yes" 
    then    
      # Configure IO_LIB_CPPFLAGS and IO_LIB_LDFLAGS
      if test x$_io_lib_config != "x"
      then
	test x"$IO_LIB_CPPFLAGS" = "x" && IO_LIB_CPPFLAGS=`$_io_lib_config --cflags`
        test x"$IO_LIB_LDFLAGS"  = "x" && IO_LIB_LDFLAGS=`$_io_lib_config --libs`
      else
        # defaults when io_lib-config isn't found
        test x"$IO_LIB_CPPFLAGS" = "x" && IO_LIB_CPPFLAGS="-I$withval/include"
        test x"$IO_LIB_LDFLAGS"  = "x" && IO_LIB_LDFLAGS="-L$withval/lib -lread"
      fi
      AC_DEFINE(HAVE_IO_LIB, 1, [Define to 1 if you have a working io_lib])
      AC_SUBST(IO_LIB_CPPFLAGS)
      AC_SUBST(IO_LIB_LDFLAGS)
    fi
  fi
])dnl
