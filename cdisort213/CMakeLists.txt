# Install cdisort submodule
#
#   library: libcdisort_debug.o
#   library: libcdisort_release.o
#
# Define the following symbols
#
#   pydisort::cdisort
#
# Normal usage would be:
#
#   target_link_libraries(pydisort::cdisort)

string(TOLOWER ${CMAKE_BUILD_TYPE} buildl)
string(TOUPPER ${CMAKE_BUILD_TYPE} buildu)

set(namel cdisort)
string(TOUPPER ${namel} nameu)

file(GLOB src_files *.c)

add_library(${namel}_${buildl}
    OBJECT
    ${src_files}
    )

set_target_properties(${namel}_${buildl}
    PROPERTIES
    COMPILE_FLAGS "${CMAKE_C_FLAGS_${buildu}}"
    )

target_link_libraries(${namel}_${buildl} m memorypool::memorypool)

add_library(pydisort::cdisort ALIAS ${namel}_${buildl})



if (CUDAToolkit_FOUND)
  file(GLOB cu_src_files cu_src/*.cu)

  add_library(${namel}_cuda_${buildl}
      OBJECT
      ${cu_src_files}
      locate.cu
      )

  set_target_properties(${namel}_cuda_${buildl}
      PROPERTIES
      CUDA_SEPARABLE_COMPILATION ON
      )

  target_link_libraries(${namel}_cuda_${buildl}
      PUBLIC
      $<IF:$<BOOL:${CUDAToolkit_FOUND}>,CUDA::cudart,>
      $<IF:$<BOOL:${CUDAToolkit_FOUND}>,CUDA::cudart_static,>
      $<IF:$<BOOL:${CUDAToolkit_FOUND}>,${TORCH_CUDA_LIBRARY},>
      $<IF:$<BOOL:${CUDAToolkit_FOUND}>,${C10_CUDA_LIBRARY},>
    m memorypool::memorypool_cu
  )

  target_include_directories(${namel}_cuda_${buildl}
      PUBLIC
      ${CMAKE_CURRENT_SOURCE_DIR}
  )

  add_library(pydisort::cdisort_cu ALIAS ${namel}_cuda_${buildl})
endif()
