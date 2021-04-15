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
