cmake_minimum_required(VERSION 3.11)
project(dibella)
set(CMAKE_CXX_STANDARD 14)

add_subdirectory($ENV{BLOOM_HOME})

include_directories(
	$ENV{SEQAN_HOME}/include 
	$ENV{COMBBLAS_HOME}/CombBLAS/install/include)
	# possibly add path to boost library (not necessary on Cori, just type 'module load boost')
link_directories(
    $ENV{COMBBLAS_HOME}/CombBLAS/install/lib 
    $ENV{BLOOM_HOME}/build)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

if(DISCOVER_LIBC)
  find_library(libc c)
else()
  set(libc "")
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-std=c++11 -fpermissive -w -O3 -DSEQAN_ENABLE_DEBUG=0 -mavx2")

find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

find_package(ZLIB)

add_definitions(-DMAX_KMER_SIZE=32)
add_definitions(-DERR_THRESHOLD=2)
add_definitions(-DMAX_NUM_READS=8)
add_definitions(-DDELTACHERNOFF=0.1)

add_executable(dibella
        src/Alphabet.cpp
        include/Alphabet.hpp
        src/Constants.cpp
        include/Constants.hpp
        include/cxxopts.hpp
        include/MPIUtils.hpp
        include/Defines.hpp
        include/VectorMap.hpp
        src/DistributedPairwiseRunner.cpp
        include/DistributedPairwiseRunner.hpp
        src/DistributedFastaData.cpp
        include/DistributedFastaData.hpp
        src/FastaData.cpp
        include/FastaData.hpp
        src/HashFuncs.c
        include/HashFuncs.h
        src/main.cpp
        src/ParallelFastaReader.cpp
        include/ParallelFastaReader.hpp
        src/ParallelOps.cpp
        include/ParallelOps.hpp
        src/TraceUtils.cpp
        include/TraceUtils.hpp
        include/Types.hpp
        include/Utils.hpp
        include/pw/PairwiseFunction.hpp
        src/pw/PairwiseFunction.cpp
        include/pw/SeedExtendXdrop.hpp
        src/pw/SeedExtendXdrop.cpp
        include/pw/OverlapFinder.hpp
        src/pw/OverlapFinder.cpp
        include/AlignmentInfo.hpp
        include/ScoreMat.hpp
        src/ScoreMat.cpp
        src/MemoryChk.c
        include/MemoryChk.h
        src/Buffer.c
        include/Buffer.h
        include/pw/BandedAligner.hpp
        include/pw/FullAligner.hpp
        src/pw/BandedAligner.cpp
        src/pw/FullAligner.cpp
        include/kmer/CommonKmers.hpp
        include/kmer/KmerIntersectSR.hpp
        include/kmer/MatrixEntry.hpp
        include/kmer/KmerOps.hpp
        src/kmer/KmerOps.cpp
        src/kmer/Kmer.cpp
        include/kmer/Kmer.hpp
        include/TransitiveReductionSR.hpp
	)

target_link_libraries(dibella PUBLIC ZLIB::ZLIB CombBLAS GraphGenlib Usortlib MURMURHASH LIBBLOOM-NOUPC)

if (MPI_FOUND)
    include_directories(${MPI_INCLUDE_PATH})
    target_link_libraries(dibella PUBLIC ${MPI_LIBRARIES})
endif()

if(OPENMP_FOUND OR OpenMP_CXX_FOUND)
    # Set THREADED if OpenMP is found
    target_compile_definitions(dibella PUBLIC THREADED)
    if(TARGET OpenMP::OpenMP_CXX)
        target_link_libraries(dibella PUBLIC OpenMP::OpenMP_CXX)
    else()
        target_compile_options(dibella PUBLIC "${OpenMP_CXX_FLAGS}")
        target_link_libraries(dibella PUBLIC "${OpenMP_CXX_FLAGS}")
    endif()
endif()
