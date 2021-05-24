# Project
add_library (${PROJECT_NAME} SHARED "${PROJECT_NAME}.cpp" "${PROJECT_NAME}.h")
set_target_properties(${PROJECT_NAME} PROPERTIES CXX_STANDARD 11
                                                 PUBLIC_HEADER ${PROJECT_NAME}.h)
target_link_libraries(${PROJECT_NAME} error_tools)
target_compile_options(${PROJECT_NAME} PUBLIC)

# Abaqus UMAT interface
add_library(${UMAT} SHARED "${UMAT}.cpp" "${UMAT}.h")
set_target_properties(${UMAT} PROPERTIES CXX_STANDARD 11
                              PUBLIC_HEADER ${UMAT}.h)
                              PREFIX ""
                              SUFFIX ".o")
target_compile_options(${UMAT} PUBLIC)
target_link_libraries(${UMAT} PUBLIC ${PROJECT_NAME})

# Local builds of upstream projects require local include paths
if(NOT cmake_build_type_lower STREQUAL "release")
    target_include_directories(${PROJECT_NAME} PUBLIC
                               "${vector_tools_SOURCE_DIR}/${CPP_SRC_PATH}"
                               "${error_tools_SOURCE_DIR}/${CPP_SRC_PATH}"
                               "${abaqus_tools_SOURCE_DIR}/${CPP_SRC_PATH}")
    target_include_directories(${UMAT} PUBLIC
                               "${vector_tools_SOURCE_DIR}/${CPP_SRC_PATH}")
endif()
