# ==========================================
# Compile the following libraries:
# - cdisort: via cmake instead of makefile
# ==========================================

add_subdirectory(cdisort)

set(CDISORT_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}
    CACHE PATH 
    "cdisort include directory")
