include(FetchContent)

set(FETCHCONTENT_QUIET FALSE)

FetchContent_Declare(
    tomlpp
    URL https://github.com/marzer/tomlplusplus/archive/refs/tags/v3.3.0.tar.gz
    DOWNLOAD_EXTRACT_TIMESTAMP TRUE
    )

FetchContent_GetProperties(tomlpp)

if(NOT tomlpp_POPULATED)
    FetchContent_Populate(tomlpp)
    add_subdirectory(${tomlpp_SOURCE_DIR} ${tomlpp_BINARY_DIR})
endif()

set(TOMLPP_INCLUDE_DIR ${tomlpp_SOURCE_DIR}/include
    CACHE PATH 
    "include directory of toml++/toml.h")