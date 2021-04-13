# Project
add_library (${PROJECT_NAME} SHARED "${PROJECT_NAME}.cpp" "${PROJECT_NAME}.h")
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 11)
target_compile_options(${PROJECT_NAME} PUBLIC)
target_include_directories(${PROJECT_NAME} PUBLIC .
                           "${vector_tools_SOURCE_DIR}/${CPP_SRC_PATH}"
                           "${error_tools_SOURCE_DIR}/${CPP_SRC_PATH}"
                           "${constitutive_SOURCE_DIR}/${CPP_SRC_PATH}"
                           "${stress_tools_SOURCE_DIR}/${CPP_SRC_PATH}"
                           "${solver_tools_SOURCE_DIR}/${CPP_SRC_PATH}"
                           "${abaqus_tools_SOURCE_DIR}/${CPP_SRC_PATH}")
target_link_libraries(${PROJECT_NAME} Eigen3::Eigen
                                      error_tools
                                      constitutive_tools
                                      stress_tools
                                      solver_tools)

# Abaqus UMAT interface
add_library (${UMAT} SHARED "${UMAT}.cpp" "${UMAT}.h")
set_target_properties(${UMAT} PROPERTIES CXX_STANDARD 11
                                      PREFIX ""
                                      SUFFIX ".o")
target_compile_options(${UMAT} PUBLIC)
target_include_directories(${UMAT} PUBLIC .
                           "${vector_tools_SOURCE_DIR}/${CPP_SRC_PATH}")
target_link_libraries(${UMAT} PUBLIC ${PROJECT_NAME})

# Conditionally add abaqus integration test(s) only for current project builds. Protects downstream project builds.
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    if(BASH_PROGRAM AND ABAQUS_PROGRAM)
        # Copy abaqus source directory to build directory
        add_custom_command(TARGET ${UMAT} POST_BUILD
                           COMMAND ${CMAKE_COMMAND} -E copy_directory
                               "${PROJECT_SOURCE_DIR}/${ABAQUS_SRC_PATH}" "${PROJECT_BINARY_DIR}/${ABAQUS_SRC_PATH}"
                          )
        # Add abaqus umat integration tests that run abaqus against input files
        add_test(NAME test_abaqus_umat_integration
                 COMMAND ${BASH_PROGRAM} test_abaqus.sh ${ABAQUS_PROGRAM} single_element_c3d8 $<TARGET_FILE:${UMAT}>
                 WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/${ABAQUS_SRC_PATH}/tests
                )
    endif()
endif()
