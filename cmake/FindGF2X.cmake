# Locate GF2X
# This module defines
# GF2X_LIBRARY
# GF2X_FOUND, if false, do not try to link to OpenAL 
# GF2X_INCLUDE_DIR, where to find the headers
#

FIND_PATH(GF2X_INCLUDE_DIR gf2x.h
  HINTS
  $ENV{GF2XDIR}
  PATH_SUFFIXES GF2X include/GF2X include
  PATHS
  /usr/local
  /usr
)

FIND_LIBRARY(GF2X_LIBRARY
  NAMES gf2x
  HINTS
  $ENV{GF2XDIR}
  PATH_SUFFIXES lib64 lib libs64 libs libs/Win32 libs/Win64
  PATHS
  /usr/local
  /usr
)

MESSAGE(STATUS "GF2X libs: " ${GF2X_LIBRARY})
# handle the QUIETLY and REQUIRED arguments and set NTL_FOUND to TRUE if
# all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(GF2X DEFAULT_MSG  GF2X_LIBRARY GF2X_INCLUDE_DIR)

MARK_AS_ADVANCED(GF2X_LIBRARY GF2X_INCLUDE_DIR)