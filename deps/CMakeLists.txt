
include(FetchContent)


# Add the CJSON library
list(APPEND CMAKE_MESSAGE_CONTEXT "cjson")
FetchContent_Declare(cjson_lib
  GIT_REPOSITORY https://github.com/DaveGamble/cJSON.git
  GIT_TAG d348621ca93571343a56862df7de4ff3bc9b5667 # version 1.7.15
)
set(ENABLE_CJSON_TEST OFF CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(cjson_lib)
set(CJSON_HEADER_FILE ${cjson_lib_SOURCE_DIR}/cJSON.h)
configure_file(${CJSON_HEADER_FILE} ${CMAKE_BINARY_DIR}/include/cjson/cJSON.h)
target_include_directories(cjson
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
)
add_library(CJSON::CJSON ALIAS cjson)
list(POP_BACK CMAKE_MESSAGE_CONTEXT)

# Add OSQP 
list(APPEND CMAKE_MESSAGE_CONTEXT "osqp")
FetchContent_Declare(qdldl_lib
    GIT_REPOSITORY https://github.com/osqp/qdldl.git
    GIT_TAG 7d16b70a10a152682204d745d814b6eb63dc5cd2  # version 0.1.5
)
FetchContent_MakeAvailable(qdldl_lib)
add_library(QDLDL::QDLDL ALIAS qdldlobject)
list(POP_BACK CMAKE_MESSAGE_CONTEXT)