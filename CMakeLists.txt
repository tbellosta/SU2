cmake_minimum_required(VERSION 3.10)

set(PROJECT_NAME SU2_dev)
project(${PROJECT_NAME})

# -------------------------------------------------------------------------------------- Compiler options
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Debug)
endif()

message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
set(CMAKE_CXX_FLAGS_DEBUG "-g -fsanitize=address")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3")

# -------------------------------------------------------------------------------------- Setup paths
set(ROOT_DIR ${CMAKE_CURRENT_SOURCE_DIR}/SU2_CFD)
set(BIN_DIR ${ROOT_DIR}/bin)
set(INC_DIR ${ROOT_DIR}/include)
set(SRC_DIR ${ROOT_DIR}/src)
set(TEST_DIR ${ROOT_DIR}/test)
set(COMMON_DIR ${CMAKE_CURRENT_SOURCE_DIR}/Common)

# -------------------------------------------------------------------------------------- Source files
file(GLOB SOURCE ${SRC_DIR}/*
        ${INC_DIR}/*
        ${COMMON_DIR}/*)

include_directories(
        ${SRC_DIR}
        ${INC_DIR}
        ${COMMON_DIR})


# -------------------------------------------------------------------------------------- MorphoGen
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR})
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${BIN_DIR})
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${BIN_DIR})

add_executable(${PROJECT_NAME} ${SOURCE})
