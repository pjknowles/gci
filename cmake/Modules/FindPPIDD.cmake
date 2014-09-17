# - Try to find ppidd
# Once done this will define
#  PPIDD_FOUND - System has ppidd
#  PPIDD_INCLUDE_DIRS - The ppidd include directories
#  PPIDD_LIBRARIES - The libraries needed to use ppidd
#  PPIDD_DEFINITIONS - Compiler switches required for using ppidd

find_package(PkgConfig)
pkg_check_modules(PC_PPIDD QUIET ppidd)
set(PPIDD_DEFINITIONS ${PC_PPIDD_CFLAGS_OTHER})

find_path(PPIDD_INCLUDE_DIR ppidd_c.h
          HINTS ${PC_PPIDD_INCLUDEDIR} ${PC_PPIDD_INCLUDE_DIRS}
          PATH_SUFFIXES ppidd )

find_library(PPIDD_LIBRARY NAMES ppidd
             HINTS ${PC_PPIDD_LIBDIR} ${PC_PPIDD_LIBRARY_DIRS} )

set(PPIDD_LIBRARIES ${PPIDD_LIBRARY} )
set(PPIDD_INCLUDE_DIRS ${PPIDD_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set PPIDD_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(ppidd  DEFAULT_MSG
                                  PPIDD_LIBRARY PPIDD_INCLUDE_DIR)

mark_as_advanced(PPIDD_INCLUDE_DIR PPIDD_LIBRARY )

find_package(GA)
