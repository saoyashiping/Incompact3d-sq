# - Find the real 2decomp-fft (decomp2d) installation without auto-download.

if(NOT DEFINED DECOMP2D_ROOT AND DEFINED ENV{DECOMP2D_ROOT})
  set(DECOMP2D_ROOT "$ENV{DECOMP2D_ROOT}")
endif()

set(_decomp2d_hints)
if(DEFINED DECOMP2D_ROOT AND NOT DECOMP2D_ROOT STREQUAL "")
  list(APPEND _decomp2d_hints "${DECOMP2D_ROOT}")
endif()
list(APPEND _decomp2d_hints ${CMAKE_PREFIX_PATH})

find_package(decomp2d CONFIG QUIET HINTS ${_decomp2d_hints})

if(NOT decomp2d_FOUND)
  find_path(DECOMP2D_INCLUDE_DIR
    NAMES decomp_2d_constants.mod
    HINTS ${_decomp2d_hints}
    PATH_SUFFIXES include include/decomp2d include/2decomp-fft)

  find_library(DECOMP2D_LIBRARY
    NAMES decomp2d 2decomp_fft decomp2d_fft
    HINTS ${_decomp2d_hints}
    PATH_SUFFIXES lib lib64)

  if(DECOMP2D_INCLUDE_DIR AND DECOMP2D_LIBRARY)
    add_library(decomp2d::decomp2d UNKNOWN IMPORTED)
    set_target_properties(decomp2d::decomp2d PROPERTIES
      IMPORTED_LOCATION "${DECOMP2D_LIBRARY}"
      INTERFACE_INCLUDE_DIRECTORIES "${DECOMP2D_INCLUDE_DIR}")

    if(NOT TARGET decomp2d)
      add_library(decomp2d INTERFACE)
      target_link_libraries(decomp2d INTERFACE decomp2d::decomp2d)
      target_include_directories(decomp2d INTERFACE "${DECOMP2D_INCLUDE_DIR}")
    endif()

    set(decomp2d_FOUND TRUE)
    set(2decomp_INCLUDE_DIR "${DECOMP2D_INCLUDE_DIR}")
  endif()
endif()

if(decomp2d_FOUND)
  message(STATUS "2decomp-fft FOUND (real external install)")
else()
  message(FATAL_ERROR
    "DECOMP2D not found. Please install 2decomp-fft once and configure with "
    "-DCMAKE_PREFIX_PATH=/home/sq/opt/2decomp-fft or set DECOMP2D_ROOT.")
endif()
