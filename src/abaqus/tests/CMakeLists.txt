# Explicitly set list of input files for integration testing
set(ABAQUS_INPUT_FILES
    "../single_element_c3d8.inp"
   )
# Create variable for hardcoded integration shell script name
set(ABAQUS_INTEGRATION_TEST "test_abaqus.sh")
# Create integration test target and associated dependencies
set(ABAQUS_TEST_NAME "test_abaqus_integration")
add_custom_target(${ABAQUS_TEST_NAME} ALL
                  DEPENDS ${UMAT})
# Copy abaqus source directory to build directory
add_custom_command(TARGET ${ABAQUS_TEST_NAME} PRE_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                       "${PROJECT_SOURCE_DIR}/${ABAQUS_SRC_PATH}" "${PROJECT_BINARY_DIR}/${ABAQUS_SRC_PATH}"
                   DEPENDS ${ABAQUS_INPUT_FILES} ${ABAQUS_INTEGRATION_TEST}
                  )
# Add abaqus cpp_stub_umat integration tests that run abaqus against input files
foreach(INPUT_FILE ${ABAQUS_INPUT_FILES})
    add_test(NAME ${ABAQUS_TEST_NAME}[${INPUT_FILE}]
             COMMAND ${BASH_PROGRAM} ${ABAQUS_INTEGRATION_TEST} ${ABAQUS_PROGRAM} ${INPUT_FILE} $<TARGET_FILE:${UMAT}>
             WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/${ABAQUS_SRC_PATH}/tests
            )
endforeach()
