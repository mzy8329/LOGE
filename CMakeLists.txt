cmake_minimum_required(VERSION 3.10)
project(LOGE)

set(SRC_DIR src)
set(INCLUDE_DIR include)

file(GLOB SOURCES "${SRC_DIR}/*.cpp")

find_package(OpenCV REQUIRED)

include_directories(
    ${INCLUDE_DIR}
    ${EIGEN3_INCLUDE_DIR}
    ${OpenCV_INCLUDE_DIRS}
    )
add_executable(LOGE ${SOURCES})
target_link_libraries(LOGE ${OpenCV_LIBS})