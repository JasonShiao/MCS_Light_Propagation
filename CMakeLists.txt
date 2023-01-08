project(MCSLightPropagation)

cmake_minimum_required(VERSION 3.15)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

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

target_include_directories(${PROJECT_NAME} PRIVATE
  src
)
