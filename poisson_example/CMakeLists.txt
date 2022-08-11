cmake_minimum_required(VERSION 3.5)
project(poisson_example)
cmake_policy(SET CMP0004 NEW)

find_package(DOLFIN REQUIRED)
include(${DOLFIN_USE_FILE})

set(SRC_DIR src/)
set(SOURCES
    ${SRC_DIR}main.cpp
    ${SRC_DIR}TLV_file.hpp
    ${SRC_DIR}TLV_file.cpp
    ${SRC_DIR}forms.h
)

add_executable(${PROJECT_NAME} ${SOURCES})
target_link_libraries(${PROJECT_NAME} dolfin)


#target_include_directories(${PROJECT_NAME}
    #PUBLIC
    #$<BUILD_INTERFACE:/usr/lib/petsc/include>
    #$<BUILD_INTERFACE:/usr/lib/slepc/include>
#)
