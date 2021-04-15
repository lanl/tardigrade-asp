set(ABAQUS_TEST_NAME "test_abaqus_integration")
add_custom_target(${ABAQUS_TEST_NAME} ALL
                  DEPENDS ${UMAT})
# Copy abaqus source directory to build directory
add_custom_command(TARGET ${ABAQUS_TEST_NAME} PRE_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                       "${PROJECT_SOURCE_DIR}/${ABAQUS_SRC_PATH}" "${PROJECT_BINARY_DIR}/${ABAQUS_SRC_PATH}"
                  )
# Add abaqus umat integration tests that run abaqus against input files
add_test(NAME test_abaqus_umat_integration
         COMMAND ${BASH_PROGRAM} test_abaqus.sh ${ABAQUS_PROGRAM} single_element_c3d8 $<TARGET_FILE:${UMAT}>
         WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/${ABAQUS_SRC_PATH}/tests
        )
