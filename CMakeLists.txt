cmake_minimum_required(VERSION 3.20)
project(PSLP LANGUAGES C CXX) # The library is written in C, but we also have some C++ test files

set(PSLP_VERSION_MAJOR 0)
set(PSLP_VERSION_MINOR 0)
set(PSLP_VERSION_PATCH 4)
set(PSLP_VERSION "${PSLP_VERSION_MAJOR}.${PSLP_VERSION_MINOR}.${PSLP_VERSION_PATCH}")
add_compile_definitions(PSLP_VERSION="${PSLP_VERSION}")

message(STATUS "Configuring PSLP (version ${PSLP_VERSION})")

include(GNUInstallDirs)

#Enable export of all symbols on Windows for shared libraries
if(BUILD_SHARED_LIBS)
  set(CMAKE_WINDOWS_EXPORT_ALL_SYMBOLS ON)
endif()

#Build position - independent code
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

#Disable C compiler extensions
set(CMAKE_C_EXTENSIONS OFF)

#This is needed for clock_gettime on Linux without compiler extensions
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
    add_compile_definitions(_POSIX_C_SOURCE=200809L)
endif()

#Append custom CMake modules
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

#Option for shared / static libraries
option(BUILD_SHARED_LIBS "Build libraries as shared as opposed to static" ON)

#Enable RPATH support
include(AddInstallRPATHSupport)
add_install_rpath_support(
  BIN_DIRS "${CMAKE_INSTALL_FULL_BINDIR}"
  LIB_DIRS "${CMAKE_INSTALL_FULL_LIBDIR}"
  INSTALL_NAME_DIR "${CMAKE_INSTALL_FULL_LIBDIR}"
  USE_LINK_PATH
)

#Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Add uninstall target
include(AddUninstallTarget)

#Specify the C standard
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

#Collect source files
file(GLOB CORE_SOURCES CONFIGURE_DEPENDS
    src/core/*.c
    src/explorers/*.c
)

# distinguish between public and private headers
set(PUBLIC_HEADERS
  include/PSLP/PSLP_API.h
 
)

file(GLOB PRIVATE_HEADERS CONFIGURE_DEPENDS
    include/core/*.h
    include/data_structures/*.h
    include/explorers/*.h
)


# Handle debugger conditionally
option(ENABLE_DEBUGGER "Enable debugger functionality" OFF)

if(NOT (ENABLE_DEBUGGER OR CMAKE_BUILD_TYPE STREQUAL "Debug"))
    list(REMOVE_ITEM CORE_SOURCES 
        "${CMAKE_CURRENT_SOURCE_DIR}/src/core/Debugger.c"
    )
    message(STATUS "Debugger disabled")
else()
    message(STATUS "Debugger enabled")
endif()

# Create PSLP library
add_library(PSLP ${CORE_SOURCES})
target_include_directories(PSLP
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/PSLP>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/PSLP>
  PRIVATE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/core>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/data_structures>
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include/explorers>
)

# Common compile options
if(MSVC)
  target_compile_options(PSLP PRIVATE /W4)
  set(CMAKE_DEBUG_POSTFIX "d")
  set_target_properties(PSLP PROPERTIES DEBUG_POSTFIX "d")
  target_compile_options(PSLP PRIVATE
    $<$<CONFIG:Debug>:/Od>
    $<$<CONFIG:Release>:/O2 /DNDEBUG>
    $<$<CONFIG:RelWithDebInfo>:/O2 /DNDEBUG>
    $<$<CONFIG:MinSizeRel>:/O1 /DNDEBUG>
  )
else()
  target_compile_options(PSLP PRIVATE -Wall -Wextra -Wpedantic -Wno-sign-compare)
  target_compile_options(PSLP PRIVATE
    $<$<CONFIG:Debug>:-g -O0>
    $<$<CONFIG:Release>:-O3 -DNDEBUG>
    $<$<CONFIG:RelWithDebInfo>:-O2 -g -DNDEBUG>
    $<$<CONFIG:MinSizeRel>:-Os -DNDEBUG>
  )
endif()

if(ENABLE_DEBUGGER OR CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_compile_definitions(PSLP PRIVATE DEBUGGER_ENABLED)
else()
    target_compile_definitions(PSLP PRIVATE NDEBUG)
endif()


# Link pthread library or compatibility layer
if(NOT WIN32)
    find_package(Threads REQUIRED)
    target_link_libraries(PSLP PRIVATE Threads::Threads)
endif()

# Testing
option(PSLP_BUILD_TESTING "Build PSLP tests" OFF)
if(PSLP_BUILD_TESTING)
  include(CTest)
  enable_testing()
  message(STATUS "PSLP testing enabled (this should be disabled for production builds!!!)")

  add_executable(all_tests tests/all_tests.c tests/test_macros.c)

  # Link PSLP library
  target_link_libraries(all_tests PRIVATE PSLP)
  target_compile_definitions(PSLP PRIVATE TESTING)

    # Add include directories explicitly for the test target
  target_include_directories(all_tests PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/include/core
      ${CMAKE_CURRENT_SOURCE_DIR}/include/data_structures
      ${CMAKE_CURRENT_SOURCE_DIR}/include/explorers
      ${CMAKE_CURRENT_SOURCE_DIR}/include/PSLP
  )

  target_compile_definitions(all_tests PRIVATE TESTING)

  add_test(NAME all_tests COMMAND all_tests)
else()
  message(STATUS "Testing disabled")
endif()

# Link math library on UNIX
if(UNIX)
  target_link_libraries(PSLP PRIVATE m)
endif()

option(COVERAGE "Enable coverage reporting" OFF)

if(COVERAGE)
    message(STATUS "Building with coverage flags")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -g --coverage")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} --coverage")
endif()

add_executable(test_cpp EXCLUDE_FROM_ALL tests/test_cpp.cpp)
target_link_libraries(test_cpp PRIVATE PSLP)


# Set output directories to match MSVC/Visual Studio defaults
foreach(OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES})
  string(TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG_UPPER)
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG_UPPER} "${CMAKE_BINARY_DIR}/${OUTPUTCONFIG}")
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG_UPPER} "${CMAKE_BINARY_DIR}/${OUTPUTCONFIG}")
  set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG_UPPER} "${CMAKE_BINARY_DIR}/${OUTPUTCONFIG}")
endforeach()

# Also set for single-config generators
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")

# Debug output for library location
message(STATUS "PSLP library will be output to: ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
message(STATUS "PSLP archive will be output to: ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}")
if(MSVC)
  message(STATUS "Expected Debug import library: ${CMAKE_BINARY_DIR}/lib/PSLPd.lib")
endif()

# Print all library and archive output directories for all configurations
foreach(OUTPUTCONFIG ${CMAKE_CONFIGURATION_TYPES})
  string(TOUPPER ${OUTPUTCONFIG} OUTPUTCONFIG_UPPER)
  message(STATUS "[DEBUG] CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG_UPPER}: ${CMAKE_RUNTIME_OUTPUT_DIRECTORY_${OUTPUTCONFIG_UPPER}}")
  message(STATUS "[DEBUG] CMAKE_LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG_UPPER}: ${CMAKE_LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG_UPPER}}")
  message(STATUS "[DEBUG] CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG_UPPER}: ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY_${OUTPUTCONFIG_UPPER}}")
  message(STATUS "[DEBUG] Expected import lib: ${CMAKE_LIBRARY_OUTPUT_DIRECTORY_${OUTPUTCONFIG_UPPER}}/PSLPd.lib")
endforeach()

# Print single-config generator output directories
message(STATUS "[DEBUG] CMAKE_RUNTIME_OUTPUT_DIRECTORY: ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
message(STATUS "[DEBUG] CMAKE_LIBRARY_OUTPUT_DIRECTORY: ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
message(STATUS "[DEBUG] CMAKE_ARCHIVE_OUTPUT_DIRECTORY: ${CMAKE_ARCHIVE_OUTPUT_DIRECTORY}")
if(MSVC)
  message(STATUS "[DEBUG] Expected Debug import library: ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}/PSLPd.lib")
endif()

# Installation and export
install(TARGETS PSLP
  EXPORT PSLPTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/PSLP
)

add_custom_command(TARGET PSLP POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E echo "[DEBUG] Listing contents of $<TARGET_FILE_DIR:PSLP>"
  COMMAND ${CMAKE_COMMAND} -E env ls -l $<TARGET_FILE_DIR:PSLP>
)

install(FILES ${PUBLIC_HEADERS} DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/PSLP)

# Export the targets file
install(EXPORT PSLPTargets
  FILE PSLPTargets.cmake
  NAMESPACE PSLP::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/PSLP
)

# Create and install Config file
include(CMakePackageConfigHelpers)
configure_package_config_file(
  ${CMAKE_CURRENT_SOURCE_DIR}/cmake/PSLPConfig.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/PSLPConfig.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/PSLP
)

write_basic_package_version_file(
  ${CMAKE_CURRENT_BINARY_DIR}/PSLPConfigVersion.cmake
  VERSION ${PSLP_VERSION}
  COMPATIBILITY AnyNewerVersion
)

install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/PSLPConfig.cmake
  ${CMAKE_CURRENT_BINARY_DIR}/PSLPConfigVersion.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/PSLP
)