find_package(BLAS QUIET REQUIRED)
if(NOT BLAS_FOUND AND REQUIRED)
  message("BLAS not found but required for SuperLU")
  return()
endif(NOT BLAS_FOUND AND REQUIRED)

# look for header files
find_path(SUPERLU_INCLUDE_DIR
  NAMES supermatrix.h
  HINTS ${SUPERLU_DIR}
  PATH_SUFFIXES "superlu" "include/superlu" "include" "SRC"
)

# look for library
find_library(SUPERLU_LIBRARY
  NAMES "superlu_4.3" "superlu_4.2" "superlu_4.1" "superlu_4.0" "superlu_3.1" "superlu_3.0" "superlu"
  HINTS ${SUPERLU_DIR}
  PATH_SUFFIXES "lib" "lib64"
)

# check if version is 4.3
include(CheckCSourceCompiles)
include(CheckIncludeFiles)
set(CMAKE_REQUIRED_INCLUDES ${SUPERLU_INCLUDE_DIR})
set(CMAKE_REQUIRED_LIBRARIES ${SUPERLU_LIBRARY})
set(CMAKE_REQUIRED_LIBRARIES ${BLAS_LIBRARY})

# check if slu_ddefs.h is available
CHECK_INCLUDE_FILES("slu_ddefs.h" HAVE_SLU_DDEFS)
if(HAVE_SLU_DDEFS)
  set(my_slu_hdr "slu_ddefs.h")
else(HAVE_SLU_DDEFS)
  set(my_slu_hdr "dsp_defs.h")
endif(HAVE_SLU_DDEFS)

CHECK_C_SOURCE_COMPILES("
#include <${my_slu_hdr}>
int main(void)
{
  return SLU_DOUBLE;
}"
SUPERLU_MIN_VERSION_4_3)

if(SUPERLU_MIN_VERSION_4_3)
  set(SUPERLU_WITH_VERSION "SuperLU >= 4.3")
else()
  set(SUPERLU_WITH_VERSION "SuperLU <= 4.2, post 2005")
endif(SUPERLU_MIN_VERSION_4_3)

# check if the expansions member is available
CHECK_C_SOURCE_COMPILES("
#include <${my_slu_hdr}>
int main()
{
   static mem_usage_t tmp;
   if (sizeof(tmp.expansions))
      return 0;
  return 0;
}
"
HAVE_MEM_USAGE_T_EXPANSIONS)

# reset include paths
set(CMAKE_REQUIRED_INCLUDES "")
set(CMAKE_REQUIRED_LIBRARIES "")

# behave like a CMake module is supposed to behave
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(
  "SuperLU"
  DEFAULT_MSG
  SUPERLU_INCLUDE_DIR
  SUPERLU_LIBRARY)

mark_as_advanced(SUPERLU_INCLUDE_DIR SUPERLU_LIBRARY)

# if both headers and library are found, store results
if(SUPERLU_FOUND)
  set(SUPERLU_INCLUDE_DIRS ${SUPERLU_INCLUDE_DIR})
  set(SUPERLU_LIBRARIES    ${SUPERLU_LIBRARY})
  # log result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeOutput.log
    "Determing location of ${SUPERLU_WITH_VERSION} succeded:\n"
    "Include directory: ${SUPERLU_INCLUDE_DIR}\n"
    "Library directory: ${SUPERLU_LIBRARY}\n\n")
  set(SUPERLU_CPPFLAGS "-I${SUPERLU_INCLUDE_DIRS} -DENABLE_SUPERLU=1")
  set(SUPERLU_LIBS "-L. ${SUPERLU_LIBRARIES} ${BLAS_LIBRARIES}")
else()
  # log errornous result
  file(APPEND ${CMAKE_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/CMakeError.log
    "Determing location of SuperLU failed:\n"
    "Include directory: ${SUPERLU_INCLUDE_DIR}\n"
    "Library directory: ${SUPERLU_LIBRARY}\n\n")
endif()
