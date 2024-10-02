cmake_minimum_required(VERSION 3.0.0)

project("Subdivision Surfaces")

# Maybe stop from CMAKEing in the wrong place
if (CMAKE_BINARY_DIR STREQUAL CMAKE_SOURCE_DIR)
    message(FATAL_ERROR "Source and build directories cannot be the same. Go use the /build directory.")
endif()

### Configure output locations
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

### Compiler options
set( CMAKE_EXPORT_COMPILE_COMMANDS 1 ) # Emit a compile flags file to support completion engines 
set(CMAKE_CXX_STANDARD 17)

if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using Clang (linux or apple) or GCC
  message("Using clang/gcc compiler flags")
  add_definitions(-DUSING_GNU)
  SET(BASE_CXX_FLAGS "-std=c++17 -Wall -Wextra -g3")
  SET(DISABLED_WARNINGS " -Wno-unused-parameter -Wno-unused-variable -Wno-unused-function -Wno-deprecated-declarations -Wno-missing-braces")
  SET(TRACE_INCLUDES " -H -Wno-error=unused-command-line-argument")

  if ("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    message("Setting clang-specific options")
    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -ferror-limit=5 -fcolor-diagnostics")
    SET(CMAKE_CXX_FLAGS_DEBUG          "-fsanitize=address -fno-limit-debug-info")
  elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    SET(BASE_CXX_FLAGS "${BASE_CXX_FLAGS} -fmax-errors=5")
    message("Setting gcc-specific options")
    SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} -Wno-maybe-uninitialized -Wno-format-zero-length -Wno-unused-but-set-parameter -Wno-unused-but-set-variable")
  endif()

  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
  add_compile_options(-Wno-sign-compare -Wno-unknown-pragmas -Wno-unused-local-typedefs -Wno-catch-value)
  #SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TRACE_INCLUDES}") # uncomment if you need to track down where something is getting included from
  SET(CMAKE_CXX_FLAGS_DEBUG          "${CMAKE_CXX_FLAGS_DEBUG} -g3")
  SET(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
  SET(CMAKE_CXX_FLAGS_RELEASE        "-march=native -O3 -DNDEBUG")
  SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "MSVC")
  # using Visual Studio C++
  message("Using Visual Studio compiler flags")
  add_definitions(-DUSING_MSVC)
  
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  set(CMAKE_CXX_EXTENSIONS OFF)
  set(BASE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
  set(BASE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP") # parallel build
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4267\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4244\"")  # ignore conversion to smaller type (fires more aggressively than the gcc version, which is annoying)
  SET(DISABLED_WARNINGS "${DISABLED_WARNINGS} /wd\"4305\"")  # ignore truncation on initialization
  SET(CMAKE_CXX_FLAGS "${BASE_CXX_FLAGS} ${DISABLED_WARNINGS}")
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /MD")
  set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} /MDd")
  add_compile_options(/wd4018 /wd4068 /wd4505 /wd5030)
  add_definitions(/D "_CRT_SECURE_NO_WARNINGS")
  add_definitions(-DNOMINMAX)
  add_definitions(-D_USE_MATH_DEFINES)
  set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)

else()
  # unrecognized
  message( FATAL_ERROR "Unrecognized compiler [${CMAKE_CXX_COMPILER_ID}]" )
endif()

# find_package(Eigen3 3.3 REQUIRED NO_MODULE)



# Add polyscope
message("\n\n == CMAKE recursively building Polyscope\n")
add_subdirectory(../deps/polyscope deps/polyscope)


# Add libIGL
# With these lines commented out, we use libIGL in "header only" mode.  Uncomment to recurse and get the full package.
# (and see the libIGL build documentation)
message("\n\n == CMAKE recursively building libIGL\n")
add_subdirectory(../deps/libigl deps/libigl)

# Create an executable
file(GLOB SOURCES  "src/SphereGeneration.cpp" "src/LoopSubdivision.cpp" "src/HalfedgeDS.cpp" "src/HalfedgeBuilder.cpp" "src/main.cpp")

add_executable(
    td_7_executable
    ${SOURCES}
)

# Link settings
target_link_libraries(td_7_executable polyscope igl::core )

add_library(td_7_lib
  src/HalfedgeDS.cpp 
  src/HalfedgeDS.h
  src/HalfedgeBuilder.cpp 
  src/HalfedgeBuilder.h
  src/LoopSubdivision.cpp
  src/LoopSubdivision.h
  src/SphereGeneration.cpp
  src/SphereGeneration.h
)



target_link_libraries(td_7_lib 
  igl::core
)

#### Do the unit testing

message("\n\n == CMAKE including the googletest\n")
add_subdirectory(../deps/googletest deps/googletest)

add_executable(unit_tests tests/unittests.cpp)
target_link_libraries(unit_tests td_7_lib gtest gtest_main)

add_test(
  NAME AllTests 
  COMMAND unit_tests 
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
)