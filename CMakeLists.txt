cmake_minimum_required (VERSION 3.14)
project (cpp_stub VERSION 0.0.1)
enable_testing()

# Set common project paths
set(CPP_SRC_PATH "src/cpp")

# Set some Abaqus integration variables
set(UMAT "umat")

# Make the code position independent
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Configure version
configure_file("${CPP_SRC_PATH}/config.h.in" "${PROJECT_BINARY_DIR}/config.h")

# Find eigen
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH} $ENV{HOME}/.local $ENV{HOME}/local/ $ENV{HOME})
find_package (Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})

# Find bash (Required for abaqus integration tests)
find_program(BASH_PROGRAM bash)
if(BASH_PROGRAM)
    message(STATUS "Found bash: ${BASH_PROGRAM}")
else()
    message(STATUS "bash not found. Abaqus integration tests will be skipped")
endif()

# Find Abaqus (Required for abaqus integration tests)
find_program(ABAQUS_PROGRAM 
             NAMES abaqus abq2020 abq2019hf4 abq2019 abq2018hf11 abq2018 abq2017hf2 abq2017
             HINTS /apps/abaqus/Commands
            )
if(ABAQUS_PROGRAM)
    message(STATUS "Found Abaqus: ${ABAQUS_PROGRAM}")
else()
    message(STATUS "Abaqus not found. Abaqus integration tests will be skipped")
endif()

# Add the cmake folder to locate the FindSphinx module
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

# Find related, but external, projects
# TODO: Decide how to handle master/dev build differences. bash env variables for GIT TAG?
include(FetchContent)
FetchContent_Declare(
    error_tools
    GIT_REPOSITORY ssh://git@xcp-stash.lanl.gov:7999/mm/error_tools.git
    GIT_TAG        origin/dev
)
FetchContent_Declare(
    vector_tools
    GIT_REPOSITORY ssh://git@xcp-stash.lanl.gov:7999/mm/vector_tools.git
    GIT_TAG        origin/dev
)
FetchContent_Declare(
    constitutive_tools
    GIT_REPOSITORY ssh://git@xcp-stash.lanl.gov:7999/mm/constitutive_tools.git
    GIT_TAG        origin/dev
)
FetchContent_Declare(
    stress_tools
    GIT_REPOSITORY ssh://git@xcp-stash.lanl.gov:7999/mm/stress_tools.git
    GIT_TAG        origin/dev
)
FetchContent_Declare(
    solver_tools
    GIT_REPOSITORY ssh://git@xcp-stash.lanl.gov:7999/mm/solver_tools.git
    GIT_TAG        origin/dev
)
FetchContent_MakeAvailable(error_tools vector_tools constitutive_tools stress_tools solver_tools)

# Add project source directories
add_subdirectory ("${CPP_SRC_PATH}")

# Only add tests and documentation for current project builds. Protects downstream project builds.
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    # Find Boost. Required for tests
    find_package(Boost 1.53.0 REQUIRED COMPONENTS unit_test_framework)
    # Add tests and docs
    add_subdirectory ("${CPP_SRC_PATH}/tests")
    add_subdirectory ("docs")
    # Add Abaqus integration tests is bash and abaqus are available
    if(BASH_PROGRAM AND ABAQUS_PROGRAM)
        add_test(NAME test_abaqus
                 COMMAND ${BASH_PROGRAM} test_abaqus.sh ${ABAQUS_PROGRAM} single_element_c3d8 
                     ${CMAKE_BINARY_DIR}/${CPP_SRC_PATH}/${UMAT}.o
                 WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/src/abaqus/tests
                )
    endif()
endif()

# Set the c++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "-Wall -ansi -pedantic -O3 -fmax-errors=5 -ggdb")
