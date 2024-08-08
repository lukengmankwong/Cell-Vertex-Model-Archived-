cmake_minimum_required(VERSION 3.30)
project(cellvertexmodel)


include_directories(${PROJECT_SOURCE_DIR}/inc)

set(SOURCES
    src/main.cpp
    src/global.cpp
    src/functions.cpp
    src/vertex.cpp
    src/edge.cpp
    src/cell.cpp
)

add_executable(${PROJECT_NAME} ${SOURCES})

find_package(CGAL REQUIRED)
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()
cmake_policy(SET CMP0167 NEW)

target_link_libraries(${PROJECT_NAME} CGAL::CGAL)
