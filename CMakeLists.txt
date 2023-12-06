cmake_minimum_required(VERSION 3.16)
project(deque)

SET(CMAKE_INSTALL_RPATH "${PROJECT_SOURCE_DIR}/bin")
SET(TASK_NAME Deque)

add_compile_options(-pedantic -Werror -Wextra -fsanitize=address -fsanitize=undefined -std=c++20)

add_link_options(-pedantic -Werror -Wextra -fsanitize=address -fsanitize=undefined -std=c++20)

find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

enable_testing()
add_executable(${TASK_NAME} tests.cpp)

add_executable(deque_pt2_stress_test stress_test.cpp)

add_test(${TASK_NAME} ${Testing_SOURCE_DIR}/bin/testing)

target_link_libraries(${TASK_NAME} Threads::Threads ${GTEST_LIBRARIES} ${GMOCK_BOTH_LIBRARIES})

install (
        TARGETS ${TASK_NAME}
        RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/bin
)