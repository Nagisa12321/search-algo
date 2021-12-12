cmake_minimum_required(VERSION 3.10)

project(search_algo)
message(STATUS "Start Cmake...")

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# the .h files
add_library(search_algo_interface INTERFACE)
target_include_directories(search_algo_interface INTERFACE include)

# debug mode
add_compile_options(-g)

# add openMP
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    # set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# project files
set(PROJECT_SOURCES ${PROJECT_SOURCE_DIR}/src)
FILE(GLOB SEARCH_ALGO_SOURCE_FILES ${PROJECT_SOURCES}/*.c*)

# create the bin
foreach(X ${SEARCH_ALGO_SOURCE_FILES})
    get_filename_component(raw_name ${X} NAME)
    string(REGEX REPLACE "\\.[^.]*$" "" cooked_name ${raw_name})

    message(STATUS "make and link the test file \"${cooked_name}\"")
    add_executable(${cooked_name} ${X})
    target_link_libraries(${cooked_name} PRIVATE search_algo_interface)
endforeach()