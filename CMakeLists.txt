cmake_minimum_required(VERSION 3.18)

project(search_algo)
message(STATUS "Start Cmake...")

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# the .h files
include_directories(/usr/lib/llvm-14/lib/clang/14.0.0/include)
add_library(search_algo_interface INTERFACE)
target_include_directories(search_algo_interface INTERFACE include)

# debug mode
add_compile_options(-g)

# add openMP
find_package(MPI REQUIRED)
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# project files
set(PROJECT_SOURCES ${PROJECT_SOURCE_DIR}/src)
FILE(GLOB SEARCH_ALGO_SOURCE_FILES ${PROJECT_SOURCES}/*.c*)

# create the bin
foreach(X ${SEARCH_ALGO_SOURCE_FILES})
    get_filename_component(raw_name ${X} NAME)
    string(REGEX REPLACE "\\.[^.]*$" "" cooked_name ${raw_name})
    add_executable(${cooked_name} ${X})
    target_link_libraries(${cooked_name} PRIVATE search_algo_interface)
    target_link_libraries(${cooked_name} PRIVATE OpenMP::OpenMP_CXX)
    target_link_libraries(${cooked_name} PRIVATE MPI::MPI_CXX)
endforeach()