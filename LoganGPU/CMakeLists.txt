cmake_minimum_required(VERSION 3.14)

# if(NOT ENABLE_CUDA)
#     message(FATAL_ERROR "Trying to build LOGAN but CUDA is not enabled")
# endif()

set(CMAKE_CUDA_STANDARD 14 CACHE STRING "")
set(CMAKE_CUDA_STANDARD_REQUIRED ON CACHE BOOL "")
set(CMAKE_CUDA_EXTENSIONS OFF CACHE BOOL "")

find_package(OpenMP REQUIRED)

if(CMAKE_VERSION VERSION_GREATER_EQUAL 3.12)
    # use the CUDA_ROOT variable 
    cmake_policy(SET CMP0074 NEW) 
endif()

enable_language(CUDA)
if (${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.18)
    if(NOT CMAKE_CUDA_ARCHITECTURES)

        if (CUDA_VERSION VERSION_GREATER_EQUAL 11.0)
        set(CMAKE_CUDA_ARCHITECTURES 50 60 70 80)
        elseif (CUDA_VERSION VERSION_GREATER_EQUAL 10.0)
        set(CMAKE_CUDA_ARCHITECTURES 50 60 70)
        else()
        set(CMAKE_CUDA_ARCHITECTURES 50 60 70)
        endif()
        
        message(STATUS "CMAKE_CUDA_ARCHITECTURES=${CMAKE_CUDA_ARCHITECTURES} -- override with specific target if this build failed.")
    else()
        message(STATUS "Using specified CMAKE_CUDA_ARCHITECTURES=${CMAKE_CUDA_ARCHITECTURES}")
    endif()
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CUDA_SEPARABLE_COMPILATION ON)
set(CMAKE_CUDA_PTX_COMPILATION ON)

if (${CMAKE_VERSION} VERSION_LESS 3.18)
    find_package(CUDA REQUIRED)
    # CUDA_SELECT_NVCC_ARCH_FLAGS(ARCH_FLAGS_LIST)
    if(ARCH_FLAGS_LIST)
        STRING(REPLACE ";" " " ARCH_FLAGS "${ARCH_FLAGS_LIST}")
    else()
        foreach(_CMAKE_CUDA_COMPUTE_CAPABILITY 50 60 70)
            set(ARCH_FLAGS "${ARCH_FLAGS} -gencode arch=compute_${_CMAKE_CUDA_COMPUTE_CAPABILITY},code=sm_${_CMAKE_CUDA_COMPUTE_CAPABILITY}")
        endforeach()
        if(NOT CMAKE_CUDA_COMPUTE_CAPABILITY)
            set(CMAKE_CUDA_COMPUTE_CAPABILITY 70 CACHE STRING "")
        endif()
        message(WARNING "Could not auto-detect the GPU arch flags building for all gpu compute capabilities 3 - 7 and PTX for ${CMAKE_CUDA_COMPUTE_CAPABILITY}")
        set(ARCH_FLAGS "${ARCH_FLAGS} -gencode arch=compute_${CMAKE_CUDA_COMPUTE_CAPABILITY},code=sm_${CMAKE_CUDA_COMPUTE_CAPABILITY}")
        set(ARCH_FLAGS "${ARCH_FLAGS} -gencode arch=compute_${CMAKE_CUDA_COMPUTE_CAPABILITY},code=compute_${CMAKE_CUDA_COMPUTE_CAPABILITY}")
    endif()
endif()

message(STATUS "ARCH_FLAGS=${ARCH_FLAGS}")

set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --default-stream per-thread -Xcompiler=-Wall,${OpenMP_CXX_FLAGS} ${ARCH_FLAGS}")
message(STATUS "CMAKE_CUDA_FLAGS=${CMAKE_CUDA_FLAGS}")

add_library(LOGAN_LIBRARY_obj OBJECT RunLoganAligner.cpp)
if (${CMAKE_VERSION} VERSION_GREATER_EQUAL 3.13)
    target_link_libraries(LOGAN_LIBRARY_obj INTERFACE OpenMP::OpenMP_CXX)
endif()

set_source_files_properties(RunLoganAligner.cpp PROPERTIES LANGUAGE CUDA LINKER_LANGUAGE CUDA)
option(LOGAN_SHARED "LOGAN shared library" OFF)
option(LOGAN_STATIC "LOGAN static library" ON)

if (LOGAN_SHARED)
    add_library(LOGAN_LIBRARY_shared SHARED $<TARGET_OBJECTS:LOGAN_LIBRARY_obj>)
    set_property(TARGET LOGAN_LIBRARY_shared PROPERTY CUDA_RESOLVE_DEVICE_SYMBOLS ON)
    target_link_libraries(LOGAN_LIBRARY_shared INTERFACE OpenMP::OpenMP_CXX) 
    install(TARGETS LOGAN_LIBRARY_shared LIBRARY DESTINATION lib)
endif()

if (LOGAN_STATIC)
    add_library(LOGAN_LIBRARY_static STATIC $<TARGET_OBJECTS:LOGAN_LIBRARY_obj>)
    target_link_libraries(LOGAN_LIBRARY_static INTERFACE OpenMP::OpenMP_CXX)
    set_property(TARGET LOGAN_LIBRARY_static PROPERTY CUDA_RESOLVE_DEVICE_SYMBOLS ON)
    install(TARGETS LOGAN_LIBRARY_static ARCHIVE DESTINATION lib)
endif()
