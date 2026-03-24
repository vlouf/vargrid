cmake_minimum_required(VERSION 3.14)
project(vargrid VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(GNUInstallDirs)

# --- Find BOM libraries ---
# These are expected at CMAKE_INSTALL_PREFIX (e.g. /opt/bom).
# Pass -DCMAKE_INSTALL_PREFIX=/opt/bom when invoking cmake.
find_package(PkgConfig)
find_package(bom-util REQUIRED)
list(APPEND DEPENDENCY_LIBRARIES bom-util::bom-util)

find_package(bom-core REQUIRED)
list(APPEND DEPENDENCY_LIBRARIES bom-core::bom-core)

# fails link on el8 without this
find_package(Threads)
list(APPEND DEPENDENCY_LIBRARIES Threads::Threads)
# --- System dependencies ---

find_package(hdf5 COMPONENTS C HL CONFIG)
if (hdf5_FOUND)  
  list(APPEND DEPENDENCY_INCLUDES_PRIVATE ${HDF5_INCLUDE_DIRS})
  list(APPEND DEPENDENCY_LIBRARIES hdf5-shared)
  message("Found HDF5 via hdf5-config.cmake")
else()
  find_package(HDF5 MODULE COMPONENTS C HL)
  if (HDF5_FOUND)    
    list(APPEND DEPENDENCY_INCLUDES_PRIVATE ${HDF5_INCLUDE_DIRS})
    list(APPEND DEPENDENCY_LIBRARIES ${HDF5_LIBRARIES})
    message("Found HDF5 via FindHDF5.cmake module")  
  endif()
endif()

find_package(netCDF CONFIG)
if(netCDF_FOUND)  
  list(APPEND DEPENDENCY_INCLUDES_PRIVATE ${netCDF_INCLUDE_DIRS})
  list(APPEND DEPENDENCY_LIBRARIES netCDF::netcdf)
  message("Found netCDF via netCDFConfig.cmake")
else()
  pkg_check_modules(netCDF netcdf)
  if(netCDF_FOUND)    
    list(APPEND DEPENDENCY_INCLUDES_PRIVATE ${netCDF_INCLUDE_DIRS})
    list(APPEND DEPENDENCY_LIBRARIES ${netCDF_LIBRARIES})
    message("Found netcdf via pkg-config")  
  endif()
endif()

pkg_check_modules(PROJ REQUIRED IMPORTED_TARGET proj)
pkg_check_modules(UDUNITS2 QUIET udunits2)

# --- Source files ---
# Expected directory layout:
#
#   vargrid/
#   ├── CMakeLists.txt          (this file)
#   └── src/
#       ├── pch.h               (from layered-flow)
#       ├── array_operations.h  (from layered-flow)
#       ├── array_operations.cc (from layered-flow)
#       ├── cappi.h             (from layered-flow)
#       ├── cappi.cc            (from layered-flow)
#       ├── io.h                (from layered-flow)
#       ├── io.cc               (from layered-flow)
#       ├── metadata.h          (from layered-flow)
#       ├── metadata.cc         (from layered-flow)
#       └── vargrid/
#           ├── main.cc
#           ├── config.h
#           ├── cost_function.h
#           ├── observation_operator.h
#           ├── observation_operator.cc
#           ├── solver.h
#           ├── vargrid.h
#           └── vargrid.cc

set(LAYERED_FLOW_SOURCES
  src/io.cc
  src/cappi.cc
  src/metadata.cc
  src/array_operations.cc
)

set(VARGRID_SOURCES
  src/vargrid/main.cc
  src/vargrid/observation_operator.cc
  src/vargrid/vargrid.cc
)

# --- Executable ---
add_executable(vargrid
  ${LAYERED_FLOW_SOURCES}
  ${VARGRID_SOURCES}
)

target_include_directories(vargrid PRIVATE
  ${BOM_INCLUDE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${CMAKE_CURRENT_SOURCE_DIR}/src/vargrid
)

target_link_libraries(vargrid PRIVATE
  ${BOM_CORE_LIB}
  ${BOM_UTIL_LIB}
  ${DEPENDENCY_LIBRARIES}
  PkgConfig::PROJ
)

if(UDUNITS2_FOUND)
  target_link_libraries(vargrid PRIVATE PkgConfig::UDUNITS2)
endif()

find_package(Threads REQUIRED)
target_link_libraries(vargrid PRIVATE Threads::Threads)

# --- Install ---
install(TARGETS vargrid RUNTIME DESTINATION bin)