cmake_minimum_required(VERSION 3.13)

project( 
    solver VERSION 0.1
    DESCRIPTION "This library contains a solver which uses the Newton-Raphson method to solve systems of equations with different numerical types."
    LANGUAGES CXX
    )

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

################################################################################
# Dependencies
################################################################################

# GMP, MPFR, MPC
find_library(gmp "gmp")
foreach(lib gmp mpfr)
    find_library("${lib}" "${lib}")
    if ( ${${lib}} STREQUAL "${lib}-NOTFOUND" ) 
        message( FATAL_ERROR 
"Library ${lib} not found. One way to specify the path to libraries is by adding -DCMAKE_LIBRARY_PATH= <list-of-library-paths> to the cmake call, for example: 
cmake .. -DCMAKE_LIBRARY_PATH=\"/home/user/libs\"
If several paths are specified, they should be separated by ; .")
    endif()
endforeach()

# Boost
set(Boost_USE_STATIC_LIBS ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)

find_package(Boost 1.69 
    COMPONENTS program_options
    REQUIRED
    )

if ( Boost_FOUND ) 
    include_directories(${Boost_INCLUDE_DIRS})
endif()

#------------------------------------------------------------------------------
# Executable
add_executable(tf-ricpad src/main.cpp)
include_directories(src/include)
target_link_libraries( 
    tf-ricpad PUBLIC
    gmp
    mpfr
    Boost::boost
    Boost::program_options
    )
