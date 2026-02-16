MOFid-Rust
A robust Rust implementation of the MOFid generator. This tool wraps OpenBabel (C++) and Systre (Java) to generate unique identifiers for Metal-Organic Frameworks.



./target/release/mofid_rust --help


The total binary size is > 100 Mb on this case 


if binary size is of concern the following commands can be included within cargo.toml

[profile.release]
strip = true           # Automatically strip symbols from the binary.
opt-level = "z"        # Optimize for size.
lto = true             # Enable Link Time Optimization (removes dead code).
codegen-units = 1      # Maximize size reduction optimization (slower build time).
panic = "abort"        # Removes unwinding logic (makes crashes slightly less pretty but smaller).




1. Prerequisites
Before compiling, ensure you have the following installed on your system (commands below are for macOS via Homebrew):
Rust: curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
CMake: brew install cmake
Pkg-Config: brew install pkg-config
OpenBabel 3: brew install open-babel
Java (OpenJDK): brew install openjdk
Git: brew install git
2. Setup Java (macOS specific)
Systre (the topology calculator) requires Java. On macOS with Homebrew, you must manually link the OpenJDK installation so the system can find it.
code
Bash
# 1. Install Java
brew install openjdk

# 2. Link it to the system VM folder (Required!)
sudo ln -sfn /opt/homebrew/opt/openjdk/libexec/openjdk.jdk /Library/Java/JavaVirtualMachines/openjdk.jdk

# 3. Verify
java -version

# Should output: openjdk version "25.x.x" (or similar)
3. Compile the C++ Dependencies (sbu & tsfm_smiles)
The Rust code relies on two external C++ binaries. You must compile these from the original MOFid source code using a corrected configuration.
A. Download Source
Clone the original repository:
code
Bash
git clone https://github.com/snurr-group/mofid.git
cd mofid/src
B. Configure CMakeLists.txt
The original CMake file is outdated. Create a new CMakeLists.txt in the src folder with the content below.
Critical: Check your OpenBabel version path (ls /opt/homebrew/share/openbabel/) and update the 3.1.0 in the file below if necessary.
code
Cmake
cmake_minimum_required(VERSION 3.10)
project(MOFid_Builder)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 1. Find OpenBabel
find_package(PkgConfig REQUIRED)
pkg_check_modules(OPENBABEL REQUIRED openbabel-3)
include_directories(${OPENBABEL_INCLUDE_DIRS})
link_directories(${OPENBABEL_LIBRARY_DIRS})

# 2. Hardcode Paths (Fixes "Unable to find plugins" error)
# Ensure these match your actual Homebrew paths!
add_definitions(-DLOCAL_OB_DATADIR="/opt/homebrew/share/openbabel/3.1.0")
add_definitions(-DLOCAL_OB_LIBDIR="/opt/homebrew/lib/openbabel/3.1.0")

# 3. Generate Config Header
file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/config_sbu.h 
"#define MOFID_VERSION_MAJOR 1
#define MOFID_VERSION_MINOR 0
")
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# 4. Compile sbu (Excluding conflicting mains)
file(GLOB ALL_SOURCES "*.cpp")
list(REMOVE_ITEM ALL_SOURCES 
    "${CMAKE_CURRENT_SOURCE_DIR}/sbu.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/tsfm_smiles.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/searchdb.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/compare.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/sobgrep.cpp"
    "${CMAKE_CURRENT_SOURCE_DIR}/remove_solvent.cpp" 
)
add_executable(sbu sbu.cpp ${ALL_SOURCES})
target_link_libraries(sbu ${OPENBABEL_LIBRARIES})

# 5. Compile tsfm_smiles
if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/tsfm_smiles.cpp")
    add_executable(tsfm_smiles tsfm_smiles.cpp)
    target_link_libraries(tsfm_smiles ${OPENBABEL_LIBRARIES})
endif()
C. Build Binaries
code
Bash
mkdir build
cd build
cmake ..
make
# Result: 'sbu' and 'tsfm_smiles' binaries are created in this folder.
4. Project Installation & Structure
Create your Rust project directory.
Create bin and resources folders.
Copy the compiled binaries and the Systre resources into place.
Your directory must look exactly like this:
code
Text
mofid_rust/
├── Cargo.toml
├── src/
│   └── main.rs
├── bin/
│   ├── sbu              <-- Copied from C++ build step above
│   └── tsfm_smiles      <-- Copied from C++ build step above
└── resources/
    ├── Systre-experimental-20.8.0.jar  <-- From original mofid/resources/
    └── RCSRnets.arc                    <-- From original mofid/resources/
OpenBabel Data Link:
The Rust wrapper specifically looks for a data folder in the working directory. Create a symbolic link to your system's OpenBabel data:
code
Bash
# Run this in the root of 'mofid_rust'
ln -s /opt/homebrew/share/openbabel/3.1.0 ./data
Build Rust Project:
code
Bash
cargo build --release
5. Usage
Run the compiled Rust binary from the project root.
Generate MOFid (.cif input)
Calculates the ID, SMILES, and Topology.
code
Bash
./target/release/mofid_rust run-mofid 01_MIL-88.cif
Use --json for structured output.
Export Tables
Converts a file containing a list of MOFids into analysis-ready TSV tables.
code
Bash
./target/release/mofid_rust export-tables --input list_of_ids.txt --output-dir tables/
Remove Metals
Filters a TSV file to remove entries containing metal atoms (uses OpenBabel detection).
code
Bash
./target/release/mofid_rust remove-metals --input tables/smiles.tsv
6. Troubleshooting

Error: MOFid-v1.TIMEOUT
Cause: Java is not installed, or the path to Systre...jar is wrong.
Fix: Run java -version. If it fails, see "Setup Java" section above.


Error: Unable to find OpenBabel plugins
Cause: The C++ sbu binary was compiled with generic paths.
Fix: Update the CMakeLists.txt to point to /opt/homebrew/lib/openbabel/3.1.0 specifically and recompile sbu.


Error: dyld: Library not loaded
Cause: System cannot find shared libraries.
Fix: export DYLD_LIBRARY_PATH=/opt/homebrew/lib:$DYLD_LIBRARY_PATH

