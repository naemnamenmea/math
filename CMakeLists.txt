cmake_minimum_required(VERSION 3.27)
project(dummy)

include(FetchContent)

FetchContent_Declare(
    misc
    GIT_REPOSITORY https://github.com/naemnamenmea/misc.git
    GIT_TAG        main
)
FetchContent_MakeAvailable(misc)

if (MSVC)
    # warning level 4
    add_compile_options(/W4 /WX)
else()
    # additional warnings
    add_compile_options(-Wall -Wextra -Wpedantic -Werror)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_SUPPRESS_REGENERATION OFF)

file(GLOB SOURCE_FILES
    public/*.hpp
    src/*.cpp
)

add_library(math STATIC ${SOURCE_FILES})

target_include_directories(math PUBLIC public)

target_link_libraries(math PUBLIC misc)
