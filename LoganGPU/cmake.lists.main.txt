cmake_minimum_required(VERSION 3.14 FATAL_ERROR)

project(logan LANGUAGES C CXX CUDA)

set(CMAKE_CXX_STANDARD 14 CACHE STRING "")
set(CMAKE_CXX_EXTENSIONS OFF CACHE BOOL "")
set(CMAKE_CUDA_STANDARD 14 CACHE STRING "")
set(CMAKE_CUDA_EXTENSIONS OFF CACHE BOOL "")

find_package(OpenMP REQUIRED)

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.12)
  cmake_policy(SET CMP0074 NEW) 
endif()

if (${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.18)
  if(NOT CMAKE_CUDA_ARCHITECTURES)

    if (CUDA_VERSION VERSION_GREATER_EQUAL 11.0)
      set(CMAKE_CUDA_ARCHITECTURES 60 70 80)
    elseif (CUDA_VERSION VERSION_GREATER_EQUAL 10.0)
      set(CMAKE_CUDA_ARCHITECTURES 60 70)
    else()
      set(CMAKE_CUDA_ARCHITECTURES 50 60 70)
    endif()
    
    message(STATUS "Setting CMAKE_CUDA_ARCHITECTURES=${CMAKE_CUDA_ARCHITECTURES} -- override with specific target if this build fails")
  else()
    message(STATUS "Using specified CMAKE_CUDA_ARCHITECTURES=${CMAKE_CUDA_ARCHITECTURES}")
  endif()
endif()

file(GLOB sources ${PROJECT_SOURCE_DIR}/main.cpp)
file(GLOB headers ${PROJECT_SOURCE_DIR}/RunLoganAligner.hpp interface.hpp logan.hpp score.hpp seed.hpp functions.hpp functions_static.hpp)

set_source_files_properties(${sources} PROPERTIES
    LANGUAGE CUDA
    LINKER_LANGUAGE CUDA)

add_executable(logan ${sources} ${headers})
target_compile_options(logan PUBLIC
    $<$<COMPILE_LANGUAGE:CUDA>:-gencode=arch=compute_52,code=sm_52>
    $<$<COMPILE_LANGUAGE:CUDA>:-gencode=arch=compute_60,code=sm_60>
    $<$<COMPILE_LANGUAGE:CUDA>:-gencode=arch=compute_61,code=sm_61>
    $<$<COMPILE_LANGUAGE:CUDA>:-gencode=arch=compute_70,code=sm_70>
    $<$<COMPILE_LANGUAGE:CUDA>:-gencode=arch=compute_70,code=compute_70>
)

target_link_libraries(logan PUBLIC OpenMP::OpenMP_CXX)
target_include_directories(logan PUBLIC ${PROJECT_SOURCE_DIR})

get_property(LANGUAGES GLOBAL PROPERTY ENABLED_LANGUAGES)
if(CMAKE_CUDA_COMPILER AND "CUDA" IN_LIST LANGUAGES)
    target_compile_options(logan PUBLIC $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler -fopenmp -Wno-deprecated-gpu-targets>)
endif()
