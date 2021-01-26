# Project
add_library (${PROJECT_NAME} "${PROJECT_NAME}.cpp" "${PROJECT_NAME}.h")
target_include_directories(${PROJECT_NAME} PUBLIC .)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 11)
target_compile_options(${PROJECT_NAME} PUBLIC)

# Abaqus UMAT interface
add_library (umat "umat.cpp" "umat.h")
set_target_properties(umat PROPERTIES CXX_STANDARD 11)
target_compile_options(umat PUBLIC)

target_link_libraries(${PROJECT_NAME} Eigen3::Eigen 
                                      error_tools 
                                      constitutive_tools
                                      stress_tools
                                      solver_tools)
target_link_libraries(umat ${PROJECT_NAME})

target_include_directories(${PROJECT_NAME} PUBLIC
                           "${vector_tools_SOURCE_DIR}/src/cpp"
                           "${error_tools_SOURCE_DIR}/src/cpp"
                           "${constitutive_SOURCE_DIR}/src/cpp"
                           "${stress_tools_SOURCE_DIR}/src/cpp"
                           "${solver_tools_SOURCE_DIR}/src/cpp")
