cmake_minimum_required(VERSION 3.11)
project(lbl_pisa)
set(CMAKE_CXX_STANDARD 14)

include_directories(
        "/usr/local/include"
        "~/sali/git/bitbucket/combinatorial-blas-2.0/CombBLAS/_install/include"
        "~/sali/software/seqan"
        "~/sali/software/seqan/seqan")

link_directories("~/sali/git/bitbucket/combinatorial-blas-2.0/CombBLAS/_install/lib")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-DSEQAN_ENABLE_DEBUG=0")

find_package(MPI REQUIRED)

add_executable(pisa
        src/Alphabet.cpp
        include/Alphabet.hpp
        src/Constants.cpp
        include/Constants.hpp
        include/cxxopts.hpp
        src/DistributedAligner.cpp
        include/DistributedAligner.hpp
        src/DistributedFastaData.cpp
        include/DistributedFastaData.h
        src/FastaData.cpp
        include/FastaData.hpp
        include/Kmer.hpp
        src/main.cpp
        src/ParallelFastaReader.cpp
        include/ParallelFastaReader.hpp
        src/ParallelOps.cpp
        include/ParallelOps.hpp
        src/TraceUtils.cpp
        include/TraceUtils.hpp
        include/Types.hpp
        include/Utils.hpp)

target_link_libraries(pisa CombBLAS GraphGenlib Usortlib)
if (MPI_FOUND)
    include_directories(${MPI_INCLUDE_PATH})
    target_link_libraries(pisa ${MPI_LIBRARIES})
endif()