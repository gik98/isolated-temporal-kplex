cmake_minimum_required(VERSION 3.10)
project(ThesisPlex)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

if (CMAKE_BUILD_TYPE EQUAL "Release")
  add_definitions(-DSPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_INFO)
else()
  add_definitions(-DSPDLOG_ACTIVE_LEVEL=SPDLOG_LEVEL_DEBUG)
endif()

find_package(spdlog REQUIRED)
find_package(OpenMP)

add_compile_options(-Wall -Wextra)

# Google ASAN
set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fsanitize=address -Wno-unknown-pragmas -Wno-unused-parameter")
set (CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fsanitize=address")

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -march=native -O2")

find_package(Boost)

include_directories(./inc)
include_directories(${Boost_INCLUDE_DIRS})

file(GLOB SOURCES "src/lib/*.cpp")

# Unit testing
#option(PACKAGE_TESTS "Build the tests" ON)
#if(PACKAGE_TESTS)
#    enable_testing()
#    include(GoogleTest)
#    add_subdirectory(test)
#endif()

add_executable(plex ${SOURCES} src/main.cpp)
target_link_libraries(plex PRIVATE spdlog::spdlog)
target_link_libraries(plex PUBLIC OpenMP::OpenMP_CXX)

file(COPY datasets DESTINATION .)
