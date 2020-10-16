add_library (${PROJECT_NAME} "${PROJECT_NAME}.cpp" "${PROJECT_NAME}.h")

target_include_directories(${PROJECT_NAME} PUBLIC .)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 11)
target_compile_options(${PROJECT_NAME} PUBLIC
                       -Wall -ansi -pedantic -O3 -fmax-errors=5 -ggdb)

target_link_libraries(${PROJECT_NAME} Eigen3::Eigen 
                                      error_tools 
                                      constitutive_tools
                                      stress_tools
                                      solver_tools)

target_include_directories(${PROJECT_NAME} PUBLIC
                           "${vector_tools_SOURCE_DIR}/src/cpp"
                           "${error_tools_SOURCE_DIR}/src/cpp"
                           "${constitutive_SOURCE_DIR}/src/cpp"
                           "${stress_tools_SOURCE_DIR}/src/cpp"
                           "${solver_tools_SOURCE_DIR}/src/cpp")
