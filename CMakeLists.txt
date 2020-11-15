cmake_minimum_required(VERSION 3.0.0)
project(NQueens VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 11)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address,undefined")

include(CTest)
enable_testing()

add_executable(NQueens main.cpp printBoard.cpp)

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
