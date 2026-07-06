cmake_minimum_required(VERSION 3.14)
project(vargrid VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Default to an optimised build: the CG solver is ~10x slower unoptimised.
if (NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
  message(STATUS "No build type specified, defaulting to Release")
endif()

include(GNUInstallDirs)

# --- Dependencies ---

find_package(PkgConfig REQUIRED)
find_package(Threads REQUIRED)

# BOM libraries
find_package(bom-util REQUIRED)
find_package(bom-core REQUIRED)

# HDF5 (try config first, then module)
find_package(hdf5 COMPONENTS C HL CONFIG)
if (hdf5_FOUND)
  list(APPEND DEP_LIBS hdf5-shared)
  message(STATUS "Found HDF5 via hdf5-config.cmake")
else()
  find_package(HDF5 MODULE COMPONENTS C HL)
  if (HDF5_FOUND)
    list(APPEND DEP_INCLUDES ${HDF5_INCLUDE_DIRS})
    list(APPEND DEP_LIBS ${HDF5_LIBRARIES})
    message(STATUS "Found HDF5 via FindHDF5.cmake")
  endif()
endif()

# NetCDF (try config first, then pkg-config)
find_package(netCDF CONFIG)
if (netCDF_FOUND)
  list(APPEND DEP_LIBS netCDF::netcdf)
  message(STATUS "Found netCDF via netCDFConfig.cmake")
else()
  pkg_check_modules(netCDF netcdf)
  if (netCDF_FOUND)
    list(APPEND DEP_INCLUDES ${netCDF_INCLUDE_DIRS})
    list(APPEND DEP_LIBS ${netCDF_LIBRARIES})
    message(STATUS "Found netCDF via pkg-config")
  endif()
endif()

pkg_check_modules(PROJ REQUIRED IMPORTED_TARGET proj)
pkg_check_modules(UDUNITS2 QUIET udunits2)

# --- Sources ---

set(SOURCES
  src/main.cc
  src/reader.cc
  src/writer.cc
  src/grid.cc
  src/gridding/cappi.cc
  src/gridding/leroi.cc
  src/gridding/observation_operator.cc
  src/gridding/variational.cc
  src/post/pipeline.cc
  src/post/steiner.cc
)

# --- Executable ---

add_executable(vargrid ${SOURCES})

target_compile_definitions(vargrid PRIVATE VARGRID_VERSION="${PROJECT_VERSION}")

target_include_directories(vargrid PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${DEP_INCLUDES}
)

target_link_libraries(vargrid PRIVATE
  bom-util::bom-util
  bom-core::bom-core
  ${DEP_LIBS}
  PkgConfig::PROJ
  Threads::Threads
)

if (UDUNITS2_FOUND)
  target_link_libraries(vargrid PRIVATE PkgConfig::UDUNITS2)
endif()

# --- Tests (dependency-free, off by default) ---

option(VARGRID_BUILD_TESTS "Build unit tests" OFF)
if (VARGRID_BUILD_TESTS)
  enable_testing()
  add_executable(test_leroi_core tests/test_leroi_core.cc)
  target_include_directories(test_leroi_core PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
  add_test(NAME leroi_core COMMAND test_leroi_core)
endif()

# --- Install ---

install(TARGETS vargrid RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR})
