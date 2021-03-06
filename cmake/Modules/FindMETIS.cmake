# -*-cmake-*-
# - Try to find the libMETIS graph partioning library
# Once done this will define:
#  METIS_FOUND        - system has the libMETIS graph partioning library
#  METIS_INCLUDE_DIR  - incude paths to use libMETIS
#  METIS_LIBRARIES    - Link these to use libMETIS
Include(EwomsMacros)

EwomsSetup("METIS" "METIS" "METIS")

set(MyIncludeSuffixes "METISLib")
#set(MyLibSuffixes "build/Linux-x86_64/" "build/Linux-i686/" "GKlib/builds/Linux-x86_64/" "GKlib/builds/Linux-i686/")

EwomsAddPathSuffixes("${MyIncludeSuffixes}" "${MyLibSuffixes}")
#EwomsAddPathSuffixes("${MyIncludeSuffixes}" "")

EwomsFindIncludeDir("metis.h")
EwomsFindLibrary("metis")
#EwomsFindLibrary("GKlib")

#EwomsRequiredLibsFound("metis" "GKlib")
EwomsRequiredLibsFound("metis")
EwomsIncludeDirsFound()
EwomsCheckFound()
