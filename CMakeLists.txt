cmake_minimum_required (VERSION 3.20...3.31)
project (nano CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
 
list (APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_SOURCE_DIR}/seqan3/cmake")
 
find_package (seqan3 3.0 REQUIRED)

add_subdirectory(source)