cmake_minimum_required(VERSION 3.6)
project(mpn C)

# Get the current git working branch
execute_process(
        COMMAND git log -1 --pretty='%h'
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE COMMIT
        OUTPUT_STRIP_TRAILING_WHITESPACE)

add_definitions("-DGITCOMMIT=\"${COMMIT}\"")

set(CMAKE_C_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_C_FLAGS "-Wall -Wextra -Wshadow -Wpedantic")
set(CMAKE_C_FLAGS_DEBUG "-g")
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG")

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIR})

find_package(Curses REQUIRED)
include_directories(${CURSES_INCLUDE_DIR})

add_executable(mpn mpn.c mpn.h amatrix.c amatrix.h auxx.c auxx.h pmatrix.c pmatrix.h reader.c reader.h mc.c mc.h libprogressbar/progressbar.c libprogressbar/progressbar.h)

target_link_libraries(mpn ${GSL_LIBRARIES})
target_link_libraries(mpn ${CURSES_LIBRARIES})
target_link_libraries(mpn m)
