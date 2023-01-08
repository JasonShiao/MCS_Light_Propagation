if (POLICY CMP0048)
  cmake_policy(SET CMP0048 NEW)
endif (POLICY CMP0048)

project(MCSLightPropagation VERSION 1.0.0)

cmake_minimum_required(VERSION 3.15)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# fetch latest argparse
include(FetchContent)
FetchContent_Declare(
    argparse
    GIT_REPOSITORY https://github.com/p-ranav/argparse.git
)
FetchContent_MakeAvailable(argparse)


add_executable(${PROJECT_NAME})

target_sources(${PROJECT_NAME} PRIVATE
  src/main.cpp
  src/AbsorbOnly.cpp
  src/MediaStructure.cpp
  src/Mismatch_Bound.cpp
  src/Rand.cpp
  src/Sim_Refl_Trans.cpp
  src/Store_grid_info.cpp
)

target_include_directories(${PROJECT_NAME} PRIVATE src)

target_link_libraries(${PROJECT_NAME} argparse)

#set(ARGPARSE_INSTALL OFF CACHE BOOL "Enable install")
#add_subdirectory(third_party/argparse)
#
#target_link_libraries(${PROJECT_NAME} argparse)

