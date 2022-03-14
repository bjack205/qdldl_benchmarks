cmake_minimum_required(VERSION 3.15)
list(APPEND CMAKE_MESSAGE_CONTEXT ldlbench)
project(ldlbench VERSION 0.1 LANGUAGES C)

include(CMakePrintHelpers)

##############################
# Options
##############################

##############################
# Dependencies 
##############################
add_subdirectory(deps)

##############################
# BUILD
##############################
add_subdirectory(src)

##############################
# TESTING 
##############################
enable_testing()
add_subdirectory(test)