#
# This file was generated from CMakeLists.txt.in
#
# This file was copied and modified from https://cmake.org/examples/
#

cmake_minimum_required( VERSION 3.16 )


project(
  T8CODE_AMD_EXPERIMENTS
  DESCRIPTION "A project that links against the t8code library."
  LANGUAGES C CXX
  VERSION "1.0.0"
)

# rudimentary version handling by setting it manually in this file
set ( T8CODE_AMD_EXPERIMENTS_VERSION_MAJOR 1 )
set ( T8CODE_AMD_EXPERIMENTS_VERSION_MINOR 0 )
set ( T8CODE_AMD_EXPERIMENTS_VERSION_PATCH 0 )

# Find MPI and set compilers
#find_package(MPI 3.0 COMPONENTS C REQUIRED )
#set(CMAKE_C_COMPILER ${MPI_C_COMPILER})
#set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER})

#
# Enable "make test" to run all the tests
#
include( CTest )

# Populate the thirdparty libraries
include (cmake/thirdparty.cmake)

# OpenMP and GPU offloading (provides t8_amd_enable_openmp)
include (cmake/OpenMPOffload.cmake)

# if t8code is linked against VTK, this project will
# automatically detect it and link against VTK as well.
# We set T8CODE_AMD_EXPERIMENTS_ENABLE_VTK to 1 in order to define our own VTK dependent macros in src/CMakeLists.txt
if (T8CODE_ENABLE_VTK)
    message ("t8code links against VTK with version ${T8CODE_VTK_VERSION_USED}")
    set (T8CODE_AMD_EXPERIMENTS_ENABLE_VTK "1")
endif ()

# Recurse into the "src" subdirectory. This does not actually
# cause another cmake executable to run. The same process will walk through
# the project's entire directory structure.
add_subdirectory (src)

# Add executable called "t8code_AMD_experiments_helloworld" that is built from the source file
# "t8code_AMD_experiments_main.cxx". The extensions are automatically found.
add_executable (t8code_AMD_experiments_helloworld t8code_AMD_experiments_main.cxx)

# Link the executable to the t8code_AMD_experiments library. Since the library has
# public include directories we will use those link directories when building
# t8code_AMD_experiments_demo.
target_link_libraries (t8code_AMD_experiments_helloworld LINK_PUBLIC t8code_AMD_experiments T8CODE::T8)
t8_amd_enable_openmp (t8code_AMD_experiments_helloworld)


# Add executable called "t8code_AMD_experiments_helloworld" that is built from the source file
# "t8code_AMD_experiments_main.cxx". The extensions are automatically found.
add_executable (t8code_AMD_experiments_data src/t8code_AMD_experiments_data_main.cxx)

# Link the executable to the t8code_AMD_experiments library. Since the library has
# public include directories we will use those link directories when building
# t8code_AMD_experiments_demo.
target_link_libraries (t8code_AMD_experiments_data LINK_PUBLIC t8code_AMD_experiments T8CODE::T8)
t8_amd_enable_openmp (t8code_AMD_experiments_data)

# Recurse into the "test" subdirectory. This does not actually
# cause another cmake executable to run. The same process will walk through
# the project's entire directory structure.
add_subdirectory (test)
