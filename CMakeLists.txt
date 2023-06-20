cmake_minimum_required(VERSION 3.27)
project(dummy)

include(FetchContent)

FetchContent_Declare(
    misc
    GIT_REPOSITORY https://github.com/naemnamenmea/misc.git
    GIT_TAG        main
)
FetchContent_MakeAvailable(misc)

file(GLOB SOURCE_FILES
    public/*.hpp
    src/*.cpp
)

add_library(math STATIC ${SOURCE_FILES})

target_include_directories(math PUBLIC public)

target_link_libraries(math PUBLIC misc)
