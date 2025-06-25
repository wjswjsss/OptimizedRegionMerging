cmake_minimum_required(VERSION 3.10)
project(region_merging LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# — PCL & Boost as before —
find_package(PCL REQUIRED COMPONENTS common io filters segmentation surface visualization)
find_package(Boost REQUIRED)

# — LAStools LASlib (headers & static lib live under LASlib sub-dirs) —
find_path(LAS_INCLUDE_DIR
  NAMES lasreader.hpp laswriter.hpp
  PATHS /usr/local/include/LASlib
)
find_library(LAS_LIBRARY
  NAMES LASlib
  PATHS /usr/local/lib/LASlib
)

message(STATUS "LAS_INCLUDE_DIR = ${LAS_INCLUDE_DIR}")
message(STATUS "LAS_LIBRARY     = ${LAS_LIBRARY}")

if (NOT LAS_INCLUDE_DIR OR NOT LAS_LIBRARY)
  message(FATAL_ERROR
    "Could not find LAStools LASlib!\n"
    " Headers expected in /usr/local/include/LASlib\n"
    " Library expected in   /usr/local/lib/LASlib\n"
    "You can override with -DLAS_INCLUDE_DIR=… and -DLAS_LIBRARY=…"
  )
endif()

# — Include paths —
include_directories(
  ${PROJECT_SOURCE_DIR}/include
  ${PCL_INCLUDE_DIRS}
  ${Boost_INCLUDE_DIRS}
  ${LAS_INCLUDE_DIR}
)

# — Collect sources & build executable —
file(GLOB SRCS "${PROJECT_SOURCE_DIR}/src/*.cpp")
add_executable(${PROJECT_NAME} ${SRCS})

# — Link against PCL, Boost and LASlib —
target_link_libraries(${PROJECT_NAME}
  ${PCL_LIBRARIES}
  ${Boost_LIBRARIES}
  ${LAS_LIBRARY}
)
