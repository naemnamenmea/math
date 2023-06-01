file(GLOB SOURCE_FILES
    public/*.hpp
    src/*.cpp
)

add_library(math STATIC ${SOURCE_FILES})

target_include_directories(math PUBLIC public)

target_link_libraries(math PUBLIC misc)
