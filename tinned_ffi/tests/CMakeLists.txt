cmake_minimum_required(VERSION 3.20)
project(tinned_c_tests C)

# --- Options users can override -------------------------
set(WORKSPACE_DIR "${CMAKE_CURRENT_LIST_DIR}/../.." CACHE PATH "Path to Cargo workspace root")
set(TINNED_CRATE_DIR "${CMAKE_CURRENT_LIST_DIR}/.." CACHE PATH "Path to Rust FFI crate")
set(CARGO_PROFILE "release" CACHE STRING "cargo profile: debug or release")
set(TINNED_FEATURES "c-headers" CACHE STRING "Cargo features for building/generating headers")
# --------------------------------------------------------

# Target dir
if(DEFINED ENV{CARGO_TARGET_DIR})
  set(TARGET_DIR "$ENV{CARGO_TARGET_DIR}")
else()
  set(TARGET_DIR "${WORKSPACE_DIR}/target")
endif()

# Profile handling
if (CARGO_PROFILE STREQUAL "release")
  set(CARGO_BUILD_ARGS --release)
  set(PROFILE_DIR "release")
else()
  set(CARGO_BUILD_ARGS "")     # debug build has no flag
  set(PROFILE_DIR "debug")
endif()

# Generate header via your test-based generator
set(GEN_HEADER "${TINNED_CRATE_DIR}/include/tinned.h")
add_custom_command(
  OUTPUT "${GEN_HEADER}"
  COMMAND ${CMAKE_COMMAND} -E env
          RUSTFLAGS=$ENV{RUSTFLAGS}
          cargo test
          --manifest-path "${TINNED_CRATE_DIR}/Cargo.toml"
          --features "${TINNED_FEATURES}"
          -- --exact header_gen::generate_c_header --nocapture
  BYPRODUCTS "${GEN_HEADER}"
  WORKING_DIRECTORY "${TINNED_CRATE_DIR}"
  COMMENT "Generating tinned.h via safer-ffi header builder"
  VERBATIM
)

add_custom_target(generate_header ALL DEPENDS "${GEN_HEADER}")

# Build the Rust FFI library
add_custom_target(rust_lib ALL
  COMMAND ${CMAKE_COMMAND} -E env
          RUSTFLAGS=$ENV{RUSTFLAGS}
          cargo build
          --manifest-path "${TINNED_CRATE_DIR}/Cargo.toml"
          --features "${TINNED_FEATURES}"
          ${CARGO_BUILD_ARGS}
  WORKING_DIRECTORY "${TINNED_CRATE_DIR}"
  COMMENT "Building Rust FFI library"
  VERBATIM
)

# Pick up the built cdylib from the workspace target dir
if (WIN32)
  set(RUST_RUNTIME "${TARGET_DIR}/${PROFILE_DIR}/tinned_ffi.dll")
  set(RUST_IMP_LIB "${TARGET_DIR}/${PROFILE_DIR}/tinned_ffi.lib")
  add_library(tinned_ffi SHARED IMPORTED GLOBAL)
  set_target_properties(tinned_ffi PROPERTIES
    IMPORTED_LOCATION "${RUST_RUNTIME}"
    IMPORTED_IMPLIB   "${RUST_IMP_LIB}"
  )
elseif(APPLE)
  set(RUST_RUNTIME "${TARGET_DIR}/${PROFILE_DIR}/libtinned_ffi.dylib")
  add_library(tinned_ffi SHARED IMPORTED GLOBAL)
  set_target_properties(tinned_ffi PROPERTIES
    IMPORTED_LOCATION "${RUST_RUNTIME}"
  )
else()
  set(RUST_RUNTIME "${TARGET_DIR}/${PROFILE_DIR}/libtinned_ffi.so")
  add_library(tinned_ffi SHARED IMPORTED GLOBAL)
  set_target_properties(tinned_ffi PROPERTIES
    IMPORTED_LOCATION "${RUST_RUNTIME}"
  )
endif()

add_dependencies(tinned_ffi rust_lib)

# C tests
add_executable(test_api
  test_api.c
  test_symbol.c
  # add more test_*.c files here
)

# Include the generated header & test includes
target_include_directories(test_api PRIVATE
  "${TINNED_CRATE_DIR}/include"
  "${CMAKE_CURRENT_LIST_DIR}"
)

target_link_libraries(test_api PRIVATE tinned_ffi)
add_dependencies(test_api generate_header)

# Make sure the Rust runtime library is found at run-time
if (APPLE)
  set_target_properties(test_api PROPERTIES
    BUILD_RPATH "@loader_path"
    INSTALL_RPATH "@loader_path"
  )
elseif(UNIX AND NOT APPLE)
  set_target_properties(test_api PROPERTIES
    BUILD_RPATH "\$ORIGIN"
    INSTALL_RPATH "\$ORIGIN"
  )
endif()

# Also copy the Rust runtime next to the test executable (robust across platforms)
add_custom_command(TARGET test_api POST_BUILD
  COMMAND ${CMAKE_COMMAND} -E copy_if_different
          "${RUST_RUNTIME}"
          "$<TARGET_FILE_DIR:test_api>"
  COMMENT "Copying Rust runtime next to test_api"
)