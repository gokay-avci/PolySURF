cmake_minimum_required(VERSION 3.10)
project(MOFid_Builder)

# Use C++14
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 1. Find OpenBabel
find_package(PkgConfig REQUIRED)
pkg_check_modules(OPENBABEL REQUIRED openbabel-3)

include_directories(${OPENBABEL_INCLUDE_DIRS})
link_directories(${OPENBABEL_LIBRARY_DIRS})

# 2. Define placeholder paths
add_definitions(-DLOCAL_OB_DATADIR="/opt/homebrew/share/openbabel/3.1.0")
add_definitions(-DLOCAL_OB_LIBDIR="/opt/homebrew/lib/openbabel/3.1.0")

# 3. Generate config header
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/config_sbu.h 
"#define MOFID_VERSION_MAJOR 1
#define MOFID_VERSION_MINOR 0
")
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# 4. Define the Source Files
# We grab all .cpp files, then REMOVE the ones that are standalone tools
file(GLOB ALL_SOURCES "*.cpp")

# These files contain their own main() function and must NOT be linked into sbu
list(REMOVE_ITEM ALL_SOURCES 
    "${CMAKE_CURRENT_SOURCE_DIR}/sbu.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/tsfm_smiles.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/searchdb.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/compare.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/sobgrep.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/remove_solvent.cpp" 
)

# 5. Compile sbu
# sbu needs its own source + the library sources (ALL_SOURCES)
add_executable(sbu sbu.cpp ${ALL_SOURCES})
target_link_libraries(sbu ${OPENBABEL_LIBRARIES})

# 6. Compile tsfm_smiles (Separate Binary)
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tsfm_smiles.cpp")
    add_executable(tsfm_smiles tsfm_smiles.cpp)
    target_link_libraries(tsfm_smiles ${OPENBABEL_LIBRARIES})
endif()