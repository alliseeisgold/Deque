cmake_minimum_required(VERSION 3.16)
project(deque)

SET(CMAKE_INSTALL_RPATH "${PROJECT_SOURCE_DIR}/bin")

add_compile_options(-pedantic -Werror -Wextra -fsanitize=address -fsanitize=undefined -std=c++20)

add_link_options(-pedantic -Werror -Wextra -fsanitize=address -fsanitize=undefined -std=c++20)

find_package(GTest REQUIRED)
include_directories(${GTEST_INCLUDE_DIRS})

enable_testing()
add_executable(main_tests tests.cpp)

add_executable(stress_tests stress_test.cpp)

add_test(main_tests ${Testing_SOURCE_DIR}/bin/testing)

target_link_libraries(main_tests Threads::Threads ${GTEST_LIBRARIES} ${GMOCK_BOTH_LIBRARIES})

install (
        TARGETS ${TASK_NAME}
        RUNTIME DESTINATION ${PROJECT_SOURCE_DIR}/bin
)