cmake_minimum_required(VERSION 3.11)
project(elba LANGUAGES C CXX CUDA)
set(CMAKE_CXX_STANDARD 14)

add_subdirectory($ENV{BLOOM_HOME})

include_directories(
	  $ENV{SEQAN_HOME}/include 
      $ENV{COMBBLAS_HOME}/CombBLAS/install/include
      $ENV{COMBBLAS_HOME}/LoganGPU
	)
# possibly add path to boost library (not necessary on Cori, just type 'module load boost')
add_subdirectory($ENV{COMBBLAS_HOME}/LoganGPU)
link_directories(
        $ENV{COMBBLAS_HOME}/CombBLAS/install/lib 
        $ENV{BLOOM_HOME}/build
	    $ENV{COMBBLAS_HOME}/LoganGPU/build
  )
     
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

if(DISCOVER_LIBC)
  find_library(libc c)
else()
  set(libc "")
endif()

# remove '-mavx2' when running on IBM Power 9 (Summit)
set(CMAKE_CXX_FLAGS_RELEASE "-std=c++14 -w -O3 -DSEQAN_ENABLE_DEBUG=0 -mavx2")

find_package(MPI REQUIRED)
find_package(OpenMP REQUIRED)

find_package(ZLIB)

add_executable(elba
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
        include/pw/FullAligner.hpp
        src/pw/FullAligner.cpp
        include/kmer/CommonKmers.hpp
        include/kmer/KmerIntersectSR.hpp
        include/kmer/MatrixEntry.hpp
        include/kmer/KmerOps.hpp
        src/kmer/KmerOps.cpp
        src/kmer/Kmer.cpp
        include/kmer/Kmer.hpp
        include/TransitiveReduction.hpp
        include/CC.h
        include/ContigGeneration.hpp
        include/pw/GPULoganAligner.hpp
        src/pw/GPULoganAligner.cpp
	)

# LOWER_KMER_FREQ and UPPER_KMER_FREQ based on E. coli CCS 29x (need to change for any given dataset based on error rate and depth)
set(LOWER_KMER_FREQ 20
    CACHE STRING "Lower Bound for the K-mer Frequency")

set(UPPER_KMER_FREQ 30
    CACHE STRING "Upper Bound for the K-mer Frequency")
    
set(MAX_KMER_SIZE 32
    CACHE STRING "Longest K-mer Length")

# DELTACHERNOFF can be 0.7-1 for CCS and 0.1-0.2 for CLR
set(DELTACHERNOFF 0.9
    CACHE STRING "Delta Chernoff's bound")
    
set(FUZZ 1000
    CACHE STRING "Fuzz Value for Transitive Transitive Reduction")

set_property(TARGET elba PROPERTY CUDA_RESOLVE_DEVICE_SYMBOLS OFF)
set(ELBA_LINK_LIBRARY ${ELBA_LINK_LIBRARY} LOGAN_LIBRARY_static CombBLAS LIBBLOOM-NOUPC MURMURHASH)
target_link_libraries(elba PUBLIC ${ELBA_LINK_LIBRARY})

target_compile_definitions(
  elba
  PRIVATE
    -DMAX_KMER_SIZE=${MAX_KMER_SIZE}
    -DDELTACHERNOFF=${DELTACHERNOFF}
    -DLOWER_KMER_FREQ=${LOWER_KMER_FREQ}
    -DUPPER_KMER_FREQ=${UPPER_KMER_FREQ}
    -DFUZZ=${FUZZ}
)

target_link_libraries(elba PUBLIC ZLIB::ZLIB)

if (MPI_FOUND)
    include_directories(${MPI_INCLUDE_PATH})
    target_link_libraries(elba PUBLIC ${MPI_LIBRARIES})
endif()

find_library(CUDART_LIBRARY cudart ${CMAKE_CUDA_IMPLICIT_LINK_DIRECTORIES})
if (CUDA_FOUND)
	message(STATUS "Cuda found...")
endif()
message(STATUS "${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}")
message(STATUS "${CUDART_LIBRARY}")
target_include_directories(elba PRIVATE ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES})
target_link_libraries(elba PUBLIC ${CUDART_LIBRARY})

if(OPENMP_FOUND OR OpenMP_CXX_FOUND)
    # Set THREADED if OpenMP is found
    target_compile_definitions(elba PUBLIC THREADED)
    if(TARGET OpenMP::OpenMP_CXX)
        target_link_libraries(elba PUBLIC OpenMP::OpenMP_CXX)
    else()
        target_compile_options(elba PUBLIC "${OpenMP_CXX_FLAGS}")
        target_link_libraries(elba PUBLIC "${OpenMP_CXX_FLAGS}")
    endif()
endif()
