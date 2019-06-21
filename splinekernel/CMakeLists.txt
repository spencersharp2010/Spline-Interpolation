cmake_minimum_required( VERSION 3.4.3 )

project( splinekernel CXX )

# Directories to include header files from
include_directories( inc )

# Gather source files into CMake variables
file( GLOB SOURCE_FILES src/*.cpp)
file( GLOB TEST_FILES test/*.cpp)
file( GLOB HEADER_FILES inc/*.h* )

# This enables exporting all symbols to the dll on windows
set( CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON )

# ----------------- Set up spline library -------------------------

# Add library from the collected source files. The headers are given so visual studio displays them
add_library( splinekernel SHARED ${SOURCE_FILES} ${HEADER_FILES} ) 

# specify the relative path the shared library object shall be installed to
if( WIN32 )
  install( TARGETS splinekernel RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX} )
else( )
  install( TARGETS splinekernel LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX} )
endif( )

target_link_libraries( splinekernel linalg )

# ------------------- Set up unit tests ---------------------------

# Add an executable and link to the library
add_executable( splinekernel_testrunner ${TEST_FILES} )

target_link_libraries( splinekernel_testrunner splinekernel linalg )

# specify the relative path the testrunner shall be installed to
install( TARGETS splinekernel_testrunner RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX} )

# To be able to debug in Windows we have to copy the linalg library to the splinekernel debug folder.
# This adds a post build command that after building splinekernel_testrunner copies linalg.dll over.
if( MSVC )
    add_custom_command(TARGET splinekernel_testrunner POST_BUILD        	
        COMMAND ${CMAKE_COMMAND} -E copy 					
            $<TARGET_FILE_DIR:linalg>/linalg.dll  
            $<TARGET_FILE_DIR:splinekernel_testrunner>)
endif( )


