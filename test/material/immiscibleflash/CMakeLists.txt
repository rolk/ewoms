# build target for the simple twophase lens problem
add_executable("test_immiscibleflash" test_immiscibleflash.cc)
target_link_libraries("test_immiscibleflash" ${EwomsLinkLibraries})

# add required libraries and includes to the build flags 
link_directories(${EwomsLinkDirectories})
include_directories(${EwomsIncludeDirectories})

