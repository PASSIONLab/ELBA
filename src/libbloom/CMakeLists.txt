cmake_minimum_required(VERSION 3.14)

message("Building libbloom")
add_subdirectory(murmur2)
include_directories(murmur2)

set(HIPMER_BLOOM64 ON)

if(HIPMER_BLOOM64)
add_library(LIBBLOOM-NOUPC OBJECT bloom64.c)
add_library(LIBBLOOM OBJECT bloom64.c)
else()
add_library(LIBBLOOM-NOUPC OBJECT bloom.c)
add_library(LIBBLOOM OBJECT bloom.c)
endif()

