cmake_minimum_required(VERSION 3.16)
project(zspace VERSION 1.0.0 LANGUAGES C)

set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)
set(CMAKE_C_EXTENSIONS OFF)

# Build types
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Options
option(BUILD_SHARED_LIBS "Build shared library libzspace" ON)
option(BUILD_STATIC_LIBS "Build static library libzspace_static" ON)
option(BUILD_CLI "Build zspace_card CLI" ON)
option(BUILD_TESTS "Build tests and verification harnesses" ON)
option(ENABLE_SANITIZERS "Enable AddressSanitizer/UBsan (debug only)" OFF)

# Compiler warnings
if(CMAKE_C_COMPILER_ID MATCHES "GNU|Clang")
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3 -march=native")
    endif()
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra -Wpedantic -Wshadow -Wpointer-arith -Wcast-qual -Wcast-align -Wstrict-prototypes -Wmissing-prototypes")
    find_package(OpenMP)
    if(OpenMP_C_FOUND)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        link_libraries(${OpenMP_C_LIBRARIES})
        message(STATUS "OpenMP enabled: ${OpenMP_C_FLAGS}")
    endif()
    if(ENABLE_SANITIZERS AND CMAKE_BUILD_TYPE STREQUAL "Debug")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address,undefined -fno-omit-frame-pointer")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address,undefined")
    endif()
endif()

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/generated)

# Find cfitsio
find_library(CFITSIO_LIBRARY NAMES cfitsio PATHS /usr/lib /usr/local/lib /opt/homebrew/lib $ENV{CONDA_PREFIX}/lib)
find_path(CFITSIO_INCLUDE_DIR NAMES fitsio.h PATHS /usr/include /usr/local/include /opt/homebrew/include $ENV{CONDA_PREFIX}/include)
if(CFITSIO_LIBRARY AND CFITSIO_INCLUDE_DIR)
    set(HAVE_CFITSIO 1)
    message(STATUS "cfitsio found: ${CFITSIO_INCLUDE_DIR}, ${CFITSIO_LIBRARY}")
else()
    set(HAVE_CFITSIO 0)
    message(WARNING "cfitsio not found - FITS ingestion will be disabled. Install cfitsio-dev (Linux) or cfitsio (macOS/conda).")
endif()

# Source files - single library: BLS + validator + ingestion + all kernels
file(GLOB GENERATED_BLS_KERNELS "${CMAKE_CURRENT_SOURCE_DIR}/generated/purce_src_bls_*.c")
file(GLOB GENERATED_VALIDATOR_KERNELS "${CMAKE_CURRENT_SOURCE_DIR}/generated/purce_src_*.c")
list(FILTER GENERATED_VALIDATOR_KERNELS EXCLUDE REGEX "purce_src_bls_")

set(ZSPACE_SOURCES
    src/zspace_bls.c
    src/zspace_eph.c
    src/zspace_audit.c
    src/zspace_core.c
    src/zspace_ingestion.c
    src/zspace_batch.c
    ${GENERATED_BLS_KERNELS}
    ${GENERATED_VALIDATOR_KERNELS}
)

# Generated header
set(GENERATED_HEADER "${CMAKE_CURRENT_SOURCE_DIR}/generated/purce_src.h")

# Static library - BLS only (working kernels)
if(BUILD_STATIC_LIBS)
    add_library(zspace_static STATIC ${ZSPACE_SOURCES})
    target_include_directories(zspace_static PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR}/generated)
    target_link_libraries(zspace_static m)
    set_target_properties(zspace_static PROPERTIES OUTPUT_NAME zspace)
endif()

# Shared library - BLS only
if(BUILD_SHARED_LIBS)
    add_library(zspace SHARED ${ZSPACE_SOURCES})
    target_include_directories(zspace PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR}/generated)
    target_link_libraries(zspace m)
    set_target_properties(zspace PROPERTIES OUTPUT_NAME zspace SOVERSION 1 VERSION 1.0.0)
endif()

# CLI links against BLS library
if(BUILD_CLI)
    add_executable(zspace_card src/zspace_card.c)
    target_link_libraries(zspace_card zspace)
    target_link_libraries(zspace_card m)
endif()

# Tests
if(BUILD_TESTS)
    enable_testing()
    add_executable(verify_kernels tests/verify_kernels.c)
    target_link_libraries(verify_kernels zspace)
    target_link_libraries(verify_kernels m)
    add_test(NAME verify_kernels COMMAND verify_kernels)
    add_executable(bench_bls tests/bench_bls.c)
    target_link_libraries(bench_bls zspace)
    target_link_libraries(bench_bls m)
    if(OpenMP_C_FOUND)
        target_link_libraries(bench_bls ${OpenMP_C_LIBRARIES})
    endif()
endif()

# Install
install(TARGETS zspace zspace_static
    ARCHIVE DESTINATION lib
    LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin
)
install(TARGETS zspace_card RUNTIME DESTINATION bin)
install(FILES src/zspace_core.h DESTINATION include)
install(FILES src/zspace_bls.h DESTINATION include)
install(FILES ${CMAKE_CURRENT_SOURCE_DIR}/generated/purce_src.h DESTINATION include)