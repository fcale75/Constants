cmake_minimum_required(VERSION 3.20)

project(constants_flint LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_path(FLINT_INCLUDE_DIR
  NAMES flint/flint.h
  HINTS
    /opt/homebrew/opt/flint/include
    /opt/homebrew/include
    /usr/local/include
    /usr/include
)

find_library(FLINT_LIBRARY
  NAMES flint
  HINTS
    /opt/homebrew/opt/flint/lib
    /opt/homebrew/lib
    /usr/local/lib
    /usr/lib
)

find_library(GMP_LIBRARY
  NAMES gmp
  HINTS
    /opt/homebrew/lib
    /usr/local/lib
    /usr/lib
)

find_library(MPFR_LIBRARY
  NAMES mpfr
  HINTS
    /opt/homebrew/lib
    /usr/local/lib
    /usr/lib
)

if(NOT FLINT_INCLUDE_DIR OR NOT FLINT_LIBRARY)
  message(FATAL_ERROR "FLINT not found. Install with: brew install flint")
endif()

add_executable(constants_objective
  src/constants_objective.cpp
)

target_include_directories(constants_objective PRIVATE ${FLINT_INCLUDE_DIR})
target_link_libraries(constants_objective PRIVATE ${FLINT_LIBRARY})

if(GMP_LIBRARY)
  target_link_libraries(constants_objective PRIVATE ${GMP_LIBRARY})
endif()

if(MPFR_LIBRARY)
  target_link_libraries(constants_objective PRIVATE ${MPFR_LIBRARY})
endif()
