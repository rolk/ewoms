# add build targets
ADD_EXECUTABLE("test_impes" test_impes.cc)
target_link_libraries("test_impes" ${EwomsLinkLibraries})

add_executable("test_impesadaptive" test_impesadaptive.cc)
target_link_libraries("test_impesadaptive" ${EwomsLinkLibraries})

add_executable("test_transport" test_transport.cc)
target_link_libraries("test_transport" ${EwomsLinkLibraries})

add_executable("test_mpfa2p" test_mpfa2p.cc)
target_link_libraries("test_mpfa2p" ${EwomsLinkLibraries})

# add required libraries and includes to the build flags 
link_directories(${EwomsLinkDirectories})
include_directories(${EwomsIncludeDirectories})

# make sure the grids are present in the build directory
add_custom_command(TARGET "test_transport"
                   POST_BUILD
                   COMMAND ${CMAKE_COMMAND} -E
                        copy_directory 
                           "${CMAKE_CURRENT_SOURCE_DIR}/grids"
                           "${CMAKE_CURRENT_BINARY_DIR}/grids")

