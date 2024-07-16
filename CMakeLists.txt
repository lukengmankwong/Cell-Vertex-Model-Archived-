cmake_minimum_required(VERSION 3.30)
project(cellvertexmodel)


find_package(CGAL REQUIRED)

add_executable(${PROJECT_NAME} functions.cpp vertex.cpp edge.cpp cell.cpp main.cpp)

target_link_libraries(${PROJECT_NAME} CGAL::CGAL)
