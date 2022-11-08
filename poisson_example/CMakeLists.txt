cmake_minimum_required(VERSION 3.5)
project(poisson_example)
cmake_policy(SET CMP0004 NEW)

find_package(DOLFIN REQUIRED)
include(${DOLFIN_USE_FILE})

set(SRC_DIR src/)
set(SOURCES
    ${SRC_DIR}main.cpp
    ${SRC_DIR}TDF_file.hpp
    ${SRC_DIR}TDF_file.cpp
    ${SRC_DIR}forms.h
)

add_executable(${PROJECT_NAME} ${SOURCES})
target_link_libraries(${PROJECT_NAME} dolfin)
