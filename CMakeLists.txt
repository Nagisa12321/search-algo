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

# project files
set(PROJECT_SOURCES ${PROJECT_SOURCE_DIR}/src)
FILE(GLOB SEARCH_ALGO_SOURCE_FILES ${PROJECT_SOURCES}/*.cpp)

# create the bin
foreach(X ${SEARCH_ALGO_SOURCE_FILES})
    get_filename_component(raw_name ${X} NAME)
    string(REGEX REPLACE "\\.[^.]*$" "" cooked_name ${raw_name})

    message(STATUS "make and link the test file \"${cooked_name}\"")
    add_executable(${cooked_name} ${X})
    target_link_libraries(${cooked_name} PRIVATE search_algo_interface)
endforeach()