#
# This is the CMake file for the src subfolder.
# It defines the library to build.
#

# Create a library called "t8code_AMD_experiments" which includes the source files.
# The extension is already found. Any number of sources could be listed here.
add_library (t8code_AMD_experiments t8code_AMD_experiments_helloworld.cxx t8code_AMD_experiments_data.cxx t8code_AMD_experiments_data_main.cxx)

# Link library against t8code, p4est, sc
target_link_libraries (t8code_AMD_experiments PRIVATE T8CODE::T8 )

# Make sure the compiler can find include files for our t8code_AMD_experiments library
# when other libraries or executables link to t8code_AMD_experiments
target_include_directories (t8code_AMD_experiments PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

# Add headers files to install in /include folder
install( FILES
    t8code_AMD_experiments_helloworld.hxx
    TYPE INCLUDE
)

# VTK dependent settings.
if( T8CODE_AMD_EXPERIMENTS_ENABLE_VTK )
    set ( T8CODE_AMD_EXPERIMENTS_VTK_VERSION_USED "${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}" )
    if (NOT ${T8CODE_AMD_EXPERIMENTS_VTK_VERSION_USED} STREQUAL ${T8CODE_VTK_VERSION_USED})
        message (FATAL_ERROR "VTK version mismatch. t8code uses ${T8CODE_VTK_VERSION_USED} and this project uses ${T8CODE_AMD_EXPERIMENTS_VTK_VERSION_USED}.")
    endif ()
    # Set CMake variables for lib target

    # Define T8CODE_AMD_EXPERIMENTS_VTK_VERSION_USED for internal and external use in source files.
    # Note, external libraries linking against t8code (i.e. using find_package ( T8CODE REQUIRED )) can
    # use T8CODE_AMD_EXPERIMENTS_VTK_VERSION_USED in source code and in CMake scripts (triggered by config.cmake.in).
    target_compile_definitions( t8code_AMD_experiments PUBLIC T8CODE_AMD_EXPERIMENTS_VTK_VERSION_USED="${T8CODE_AMD_EXPERIMENTS_VTK_VERSION_USED}" )
    target_compile_definitions( t8code_AMD_experiments PUBLIC T8CODE_AMD_EXPERIMENTS_ENABLE_VTK=1 )
    target_include_directories( t8code_AMD_experiments PUBLIC ${VTK_INCLUDE_DIR} )
    target_link_libraries( t8code_AMD_experiments PUBLIC ${VTK_LIBRARIES} )
    # Internal source files that need VTK linkage to compile:
#    target_sources(T8CODE_AMD_EXPERIMENTS PRIVATE
#     put/source/file/here
#    )
#    
    # Installed header files that need VTK linkage:
#    install( FILES
#     put/header/file/here
#    DESTINATION include)
endif()

# Make this package available as a cmake install package, so 
# that other projects can find it with find_package.
include( CMakePackageConfigHelpers )

configure_package_config_file( ${CMAKE_CURRENT_LIST_DIR}/config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/T8CODE_AMD_EXPERIMENTSConfig.cmake
INSTALL_DESTINATION cmake
)

write_basic_package_version_file(
${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/T8CODE_AMD_EXPERIMENTSConfigVersion.cmake
COMPATIBILITY SameMajorVersion
VERSION "1.0.0"
)

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/T8CODE_AMD_EXPERIMENTSConfig.cmake
              ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles/T8CODE_AMD_EXPERIMENTSConfigVersion.cmake
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake )

install( TARGETS t8code_AMD_experiments LIBRARY )

#
