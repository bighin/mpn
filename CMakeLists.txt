cmake_minimum_required(VERSION 3.6)
project(mpn C CXX)

#
# Register own CMake extensions:
#
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

# Get the current git working branch
execute_process(
        COMMAND git log -1 --pretty='%h'
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE COMMIT
        OUTPUT_STRIP_TRAILING_WHITESPACE)

add_definitions("-DGITCOMMIT=\"${COMMIT}\"")

set(CMAKE_CXX_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wshadow -Wpedantic")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIR})

find_package(Curses REQUIRED)
include_directories(${CURSES_INCLUDE_DIR})

set(ALPS_FORCE_NO_COMPILER_CHECK "1")
find_package(ALPSCore REQUIRED)
include_directories(${ALPSCore_INCLUDE_DIR})

add_executable(mpn mpn.c mpn.h amatrix.c amatrix.h auxx.c auxx.h pmatrix.c pmatrix.h loaderis.c loaderis.h mc.c mc.h libprogressbar/progressbar.c libprogressbar/progressbar.h inih/ini.c inih/ini.h main.c config.c config.h multiplicity.c multiplicity.h cache.c cache.h permutations.c permutations.h weight.c weight.h weight2.c weight2.h sampling.cpp sampling.h)

target_link_libraries(mpn ${GSL_LIBRARIES})
target_link_libraries(mpn ${CURSES_LIBRARIES})
target_link_libraries(mpn ${ALPSCore_LIBRARIES})
target_link_libraries(mpn m)
