# - Try to find ga
# Once done this will define
#  GA_FOUND - System has ga
#  GA_INCLUDE_DIRS - The ga include directories
#  GA_LIBRARIES - The libraries needed to use ga
#  GA_DEFINITIONS - Compiler switches required for using ga

find_package(PkgConfig)
pkg_check_modules(PC_GA QUIET ga)
set(GA_DEFINITIONS ${PC_GA_CFLAGS_OTHER})

find_path(GA_INCLUDE_DIR ga.h
          HINTS ${PC_GA_INCLUDEDIR} ${PC_GA_INCLUDE_DIRS}
          PATH_SUFFIXES ga )

find_library(GA_LIBRARY NAMES ga
             HINTS ${PC_GA_LIBDIR} ${PC_GA_LIBRARY_DIRS} )

set(GA_LIBRARIES ${GA_LIBRARY} -larmci )
set(GA_INCLUDE_DIRS ${GA_INCLUDE_DIR} )

include(FindPackageHandleStandardArgs)
# handle the QUIETLY and REQUIRED arguments and set GA_FOUND to TRUE
# if all listed variables are TRUE
find_package_handle_standard_args(ga  DEFAULT_MSG
                                  GA_LIBRARY GA_INCLUDE_DIR)

mark_as_advanced(GA_INCLUDE_DIR GA_LIBRARY )
