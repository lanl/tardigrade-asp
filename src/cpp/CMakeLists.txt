# Project
add_library (${PROJECT_NAME} SHARED "${PROJECT_NAME}.cpp" "${PROJECT_NAME}.h")
target_include_directories(${PROJECT_NAME} PUBLIC .)
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 11)
target_compile_options(${PROJECT_NAME} PUBLIC)

# Abaqus UMAT interface
add_library (umat SHARED "umat.cpp" "umat.h")
target_include_directories(umat PUBLIC .)
set_target_properties(umat PROPERTIES CXX_STANDARD 11
                                      PREFIX ""
                                      SUFFIX ".o")
target_compile_options(umat PUBLIC)


target_include_directories(${PROJECT_NAME} PUBLIC)
