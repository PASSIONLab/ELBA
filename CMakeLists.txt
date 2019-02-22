cmake_minimum_required(VERSION 3.11)
project(lbl_dal)


include_directories(
        "/usr/local/include"
        "~/sali/git/bitbucket/combinatorial-blas-2.0/CombBLAS/_install/include"
        "~/sali/software/seqan"
        "~/sali/software/seqan/seqan"
)

link_directories("~/sali/git/bitbucket/combinatorial-blas-2.0/CombBLAS/_install/lib")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

find_package(MPI REQUIRED)

include_directories(.)

add_executable(lbl_dal
        Alphabet.cpp
        Alphabet.hpp
        Constants.cpp
        Constants.hpp
        cxxopts.hpp
        TraceUtils.cpp
        TraceUtils.hpp
        DistributedFastaData.cpp
        DistributedFastaData.h
        FastaData.cpp
        FastaData.hpp
        Kmer.hpp
        main.cpp
        ParallelFastaReader.cpp
        ParallelFastaReader.hpp
        ParallelOps.cpp
        ParallelOps.hpp
        Types.hpp
        Utils.hpp DistributedAligner.cpp DistributedAligner.hpp)

target_link_libraries(lbl_dal CombBLAS GraphGenlib Usortlib)
if (MPI_FOUND)
    include_directories(${MPI_INCLUDE_PATH})
    target_link_libraries(lbl_dal ${MPI_LIBRARIES})
endif()
