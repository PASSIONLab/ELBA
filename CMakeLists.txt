cmake_minimum_required(VERSION 3.10)
project(lbl_dal)
include_directories("/usr/local/include" "~/sali/git/bitbucket/combinatorial-blas-2.0/CombBLAS/_install/include")
link_directories("~/sali/git/bitbucket/combinatorial-blas-2.0/CombBLAS/_install/lib")
find_package(MPI REQUIRED)
#find_package(Boost 1.68 COMPONENTS program_options REQUIRED)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -DDEBUG")

add_executable(lbl_dal main.cpp Constants.hpp ParallelOps.cpp ParallelOps.hpp Constants.cpp ParallelFastaReader.cpp ParallelFastaReader.hpp FastaData.cpp FastaData.hpp Alphabet.hpp Alphabet.cpp Kmer.hpp Utils.hpp cxxopts.hpp Types.hpp)

add_executable(comb_blas ComBlasTest.cpp)
add_executable(seg_test SegTest.cpp)

target_link_libraries(lbl_dal CombBLAS GraphGenlib Usortlib)
target_link_libraries(comb_blas CombBLAS GraphGenlib Usortlib)
target_link_libraries(seg_test CombBLAS)

if (MPI_FOUND)
    include_directories(${MPI_INCLUDE_PATH})
    target_link_libraries(lbl_dal ${MPI_LIBRARIES})
    target_link_libraries(comb_blas ${MPI_LIBRARIES})
    target_link_libraries(seg_test ${MPI_LIBRARIES})
endif()


