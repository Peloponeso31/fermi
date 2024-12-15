# Minimum CMake version required
cmake_minimum_required(VERSION 3.16)

# Project name and version
project(fermi VERSION 1.0)

# Set the C++ standard (e.g., C++17)
set(CMAKE_CXX_STANDARD 23)
set(CMAKE_CXX_STANDARD_REQUIRED True)
message(STATUS "C Standard: ${CMAKE_CXX_STANDARD}")

# Specify the main executable and its source files
add_executable(${PROJECT_NAME}
        main.cpp
)

# Find and link Boost libraries
find_package(Boost REQUIRED COMPONENTS program_options)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    target_link_libraries(${PROJECT_NAME} PRIVATE Boost::program_options)
endif()

# Include other library folders
# Add the path to additional libraries (change "lib" to your actual folder path)
# Example: If you have "lib/mylib" with a CMakeLists.txt, you can use add_subdirectory
#add_subdirectory(lib/mylib)

# If you just need the headers without a CMakeLists.txt in the library folder
# Include the directory to find header files
include_directories(
        "math-parser-benchmark-project"
        "tabulate/include"
)

# Link libraries (if needed)
# Example: If "mylib" generates a library target, link it to your project
# target_link_libraries(${PROJECT_NAME} PRIVATE mylib)
