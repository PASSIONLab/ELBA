cmake_minimum_required(VERSION 3.11)
project(lbl_pisa)
set(CMAKE_CXX_STANDARD 14)

#IF(CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "DEBUG")
#    SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address -fno-elide-constructors")
#ENDIF(CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "DEBUG")

# Download and unpack googletest at configure time
configure_file(CMakeLists.txt.in googletest-download/CMakeLists.txt)
execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
    message(FATAL_ERROR "CMake step for googletest failed: ${result}")
endif()
execute_process(COMMAND ${CMAKE_COMMAND} --build .
        RESULT_VARIABLE result
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/googletest-download )
if(result)
    message(FATAL_ERROR "Build step for googletest failed: ${result}")
endif()

# Prevent overriding the parent project's compiler/linker
# settings on Windows
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

# Add googletest directly to our build. This defines
# the gtest and gtest_main targets.
add_subdirectory(${CMAKE_CURRENT_BINARY_DIR}/googletest-src
        ${CMAKE_CURRENT_BINARY_DIR}/googletest-build
        EXCLUDE_FROM_ALL)

# The gtest/gtest_main targets carry header search path
# dependencies automatically when using CMake 2.8.11 or
# later. Otherwise we have to add them here ourselves.
if (CMAKE_VERSION VERSION_LESS 2.8.11)
    include_directories("${gtest_SOURCE_DIR}/include")
endif()

include_directories(
        "/usr/local/include"
        "~/sali/git/bitbucket/combinatorial-blas-2.0/CombBLAS/_install/include"
#        "~/sali/software/seqan"
#        "~/sali/software/seqan/seqan")
        ~/sali/software/seqan-v2.4.0/include/
        /usr/local/include/boost/)
#
link_directories("~/sali/git/bitbucket/combinatorial-blas-2.0/CombBLAS/_install/lib")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS_RELEASE "-DSEQAN_ENABLE_DEBUG=0")
#set(CMAKE_CXX_FLAGS_DEBUG "-DDEBUG=1")

find_package(MPI REQUIRED)

add_executable(pisa
        src/Alphabet.cpp
        include/Alphabet.hpp
        src/Constants.cpp
        include/Constants.hpp
        include/cxxopts.hpp
        src/DistributedPairwiseRunner.cpp
        include/DistributedPairwiseRunner.hpp
        src/DistributedFastaData.cpp
        include/DistributedFastaData.hpp
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
        include/pw/BandedAligner.hpp
        include/pw/FullAligner.hpp
        src/pw/BandedAligner.cpp
        src/pw/FullAligner.cpp
        include/kmer/CommonKmers.hpp
        include/kmer/KmerIntersectSR.hpp
        include/kmer/MatrixEntry.hpp
        include/kmer/KmerOps.hpp
        src/kmer/KmerOps.cpp
        include/kmer/Kmer.hpp
        include/NearestKmers2.hpp
        src/NearestKmers2.cpp
        include/kmer/SubKmerIntersectSR.hpp
        )

add_executable(nkmers
        src/NearestKmers.cpp
        src/Alphabet.cpp include/Alphabet.hpp
        src/ScoreMat.cpp include/ScoreMat.hpp)

add_executable(nkmers2
        src/NearestKmers2.cpp include/NearestKmers2.hpp
        src/Alphabet.cpp include/Alphabet.hpp
        src/ScoreMat.cpp include/ScoreMat.hpp)

add_executable(seqanvecs src/seqantest/seqanvec.cpp)
add_executable(seqanalign src/seqantest/alignertest.cpp)

add_executable(kunigen
        TestKmerUniverseGeneration.cpp
        include/kmer/KmerOps.hpp
        include/Alphabet.hpp
        include/ParallelOps.hpp
        src/ParallelOps.cpp
        src/Alphabet.cpp)

add_executable(simpletests SimpleTests.cpp include/MinMaxHeap.hpp)
target_link_libraries(simpletests gtest_main)

add_executable(mcb MatrixCreationBug.cpp)

add_executable(ipack IntPack.cpp)

target_link_libraries(pisa CombBLAS GraphGenlib Usortlib)
target_link_libraries(mcb CombBLAS GraphGenlib Usortlib)
target_link_libraries(kunigen CombBLAS GraphGenlib Usortlib)
if (MPI_FOUND)
    include_directories(${MPI_INCLUDE_PATH})
    target_link_libraries(pisa ${MPI_LIBRARIES})
    target_link_libraries(mcb ${MPI_LIBRARIES})
    target_link_libraries(kunigen ${MPI_LIBRARIES})
endif()