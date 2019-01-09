cmake_minimum_required(VERSION 3.10)
project(lbl_dal)
include_directories("/usr/local/include" "/Users/esaliya/sali/git/bitbucket/combinatorial-blas-2.0/CombBLAS/_install/include")
link_directories("/Users/esaliya/sali/git/bitbucket/combinatorial-blas-2.0/CombBLAS/_install/lib")
find_package(MPI REQUIRED)
#find_package(Boost 1.68 COMPONENTS program_options REQUIRED)

set(CMAKE_CXX_STANDARD 14)

add_executable(lbl_dal main.cpp constants.h ParallelOps.cpp ParallelOps.h constants.cpp ParallelFastaReader.cpp ParallelFastaReader.hpp FastaData.cpp FastaData.hpp)

add_executable(comb_blas ComBlasTest.cpp)
target_link_libraries(comb_blas CombBLAS GraphGenlib Usortlib)


#if(Boost_FOUND)
#    include_directories(${Boost_INCLUDE_DIR})
#    target_link_libraries(lbl_dal ${Boost_LIBRARIES})
#endif()

if (MPI_FOUND)
    include_directories(${MPI_INCLUDE_PATH})
    target_link_libraries(lbl_dal ${MPI_LIBRARIES})
    target_link_libraries(comb_blas ${MPI_LIBRARIES})
endif()


