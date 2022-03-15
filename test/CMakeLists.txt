# function add_ulqr_test(name)
#
# Adds a new executable test called <name>_test.
# Assumes thes source code is in a file called <name>_test.
function (add_project_test name)
  set(TEST_NAME ${name}_test)
  add_executable(${TEST_NAME}
    ${TEST_NAME}.c

    test_utils.h
    test_utils.c
  )
  target_link_libraries(${TEST_NAME}
    PUBLIC
    csc
    PRIVATE 
    simpletest 
    m
  )
  add_test(NAME ${TEST_NAME} COMMAND ${TEST_NAME})
  target_compile_definitions(${TEST_NAME}
    PRIVATE
    JSONDATA="${CMAKE_CURRENT_LIST_DIR}/2-1-32_banded.json"
  )
endfunction()

# Build the testing library
add_library(simpletest
  simpletest/simpletest.h
  simpletest/simpletest.c
)
target_include_directories(simpletest
  PRIVATE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
)

# Add tests
add_project_test(csc)