# Minimum cmake version
cmake_minimum_required (VERSION 3.0.0)

# Name of project and that it is C++ only.
project (core CXX)

# ----------------------------------------------------------------------------
# Dependencies
# ----------------------------------------------------------------------------

# Search for zlib as a dependency for SeqAn.
find_package (ZLIB)

# Load the SeqAn module and fail if not found.
find_package (SeqAn REQUIRED)

# ----------------------------------------------------------------------------
# Build Setup
# ----------------------------------------------------------------------------

# Add include directories.
include_directories (${SEQAN_INCLUDE_DIRS})

# Add definitions set by find_package (SeqAn).
add_definitions (${SEQAN_DEFINITIONS})

# Add CXX flags found by find_package (SeqAn).
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${SEQAN_CXX_FLAGS}")



find_package( Threads )
link_libraries( ${CMAKE_THREAD_LIBS_INIT} )



# Add executable and link against SeqAn dependencies.
add_executable (core core.cpp)
target_link_libraries (core ${SEQAN_LIBRARIES})
