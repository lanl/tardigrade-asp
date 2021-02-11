set(TEST_NAME "test_${PROJECT_NAME}")
add_executable (${TEST_NAME} "${TEST_NAME}.cpp")
add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
set_target_properties(${TEST_NAME} PROPERTIES CXX_STANDARD 11)
target_link_libraries(${TEST_NAME} PUBLIC ${PROJECT_NAME}
                      error_tools
                      constitutive_tools
                      stress_tools
                      solver_tools)

# TODO: Move the abaqus tools header and tests to a stand alone project
set(TEST_ABAQUS_TOOLS "test_${ABAQUS_TOOLS}")
add_executable (${TEST_ABAQUS_TOOLS} "${TEST_ABAQUS_TOOLS}.cpp")
add_test(NAME ${TEST_ABAQUS_TOOLS} COMMAND ${TEST_ABAQUS_TOOLS})
set_target_properties(${TEST_ABAQUS_TOOLS} PROPERTIES CXX_STANDARD 11)
target_include_directories(${TEST_ABAQUS_TOOLS} PUBLIC .
                           "${vector_tools_SOURCE_DIR}/${CPP_SRC_PATH}")
