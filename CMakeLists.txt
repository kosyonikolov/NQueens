cmake_minimum_required(VERSION 3.0.0)
project(NQueens VERSION 0.1.0)

set(CMAKE_CXX_STANDARD 11)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address,undefined")

include(CTest)
enable_testing()

option(DRAW_BOARD "Draw and save the solution board" OFF)

if (DRAW_BOARD)
    find_package(OpenCV REQUIRED)
    include_directories(${OpenCV_INCLUDE_DIRS})
    add_definitions(-DDRAW_SOLUTION_BOARD)
endif()

add_executable(NQueens main.cpp printBoard.cpp drawBoard.cpp)
if (DRAW_BOARD)
    target_link_libraries(NQueens ${OpenCV_LIBS})
endif()

set(CPACK_PROJECT_NAME ${PROJECT_NAME})
set(CPACK_PROJECT_VERSION ${PROJECT_VERSION})
include(CPack)
