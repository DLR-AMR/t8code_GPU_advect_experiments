# OpenMP and GPU offloading support.
#
# Cache option:
#   OFFLOAD_ARCH  GPU architecture(s) for --offload-arch (e.g. gfx90a, gfx942, sm_80).
#                 Empty = autodetect, NONE = disable offloading.
#
# Provides:
#   t8_amd_enable_openmp(<target>)  Links <target> against OpenMP and adds
#                                   --offload-arch if a GPU is available.

include_guard(GLOBAL)

find_package(OpenMP REQUIRED)

set(OFFLOAD_ARCH "" CACHE STRING
    "GPU architecture(s) for --offload-arch (e.g. gfx90a, gfx942, sm_80). Empty = autodetect, NONE = disable offloading.")

# Detect the GPUs in this machine with amdgpu-arch / nvptx-arch (shipped with LLVM/ROCm clang).
function(_t8_amd_detect_offload_arch out_var)
  set(_archs "")
  find_program(AMDGPU_ARCH_EXECUTABLE NAMES amdgpu-arch
               HINTS ENV ROCM_PATH /opt/rocm PATH_SUFFIXES llvm/bin bin)
  find_program(NVPTX_ARCH_EXECUTABLE NAMES nvptx-arch
               HINTS ENV ROCM_PATH /opt/rocm PATH_SUFFIXES llvm/bin bin)
  foreach(_tool IN ITEMS ${AMDGPU_ARCH_EXECUTABLE} ${NVPTX_ARCH_EXECUTABLE})
    execute_process(COMMAND ${_tool}
                    OUTPUT_VARIABLE _out
                    RESULT_VARIABLE _res
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    ERROR_QUIET)
    if(_res EQUAL 0 AND _out)
      string(REPLACE "\n" ";" _out "${_out}")
      list(APPEND _archs ${_out})
    endif()
  endforeach()
  list(FILTER _archs EXCLUDE REGEX "^gfx000$")   # CPU agent reported by ROCm
  list(REMOVE_DUPLICATES _archs)
  set(${out_var} "${_archs}" PARENT_SCOPE)
endfunction()

if(NOT OFFLOAD_ARCH)
  _t8_amd_detect_offload_arch(_offload_arch)
  if(_offload_arch)
    message(STATUS "Detected GPU offload architecture(s): ${_offload_arch}")
  else()
    message(STATUS "No GPU detected, building OpenMP without offloading")
  endif()
elseif(OFFLOAD_ARCH STREQUAL "NONE")
  set(_offload_arch "")
else()
  set(_offload_arch "${OFFLOAD_ARCH}")
endif()

if(_offload_arch AND NOT CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  message(WARNING "--offload-arch needs a Clang-based compiler (clang++/amdclang++), "
                  "but the compiler is ${CMAKE_CXX_COMPILER_ID}. Disabling offloading.")
  set(_offload_arch "")
endif()

string(REPLACE ";" "," T8_AMD_OFFLOAD_ARCH_FLAG "${_offload_arch}")

function(t8_amd_enable_openmp target)
  target_link_libraries(${target} PRIVATE OpenMP::OpenMP_CXX)
  if(T8_AMD_OFFLOAD_ARCH_FLAG)
    target_compile_options(${target} PRIVATE --offload-arch=${T8_AMD_OFFLOAD_ARCH_FLAG})
    target_link_options(${target} PRIVATE --offload-arch=${T8_AMD_OFFLOAD_ARCH_FLAG})
  endif()
endfunction()
