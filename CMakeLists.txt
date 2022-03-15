cmake_minimum_required(VERSION 3.15)
list(APPEND CMAKE_MESSAGE_CONTEXT ldlbench)
project(ldlbench VERSION 0.1 LANGUAGES C)

include(CMakePrintHelpers)

##############################
# Options
##############################
SET(PARDISO_DIR "/home/brian/Software" CACHE STRING "Directory containing pardiso library.")
SET(PARDISO_LIB "pardiso600-GNU800-X86-64" CACHE STRING "Name of pardiso library, without suffix.")

##############################
# Dependencies 
##############################
add_subdirectory(deps)

##############################
# BUILD
##############################
# Compile options
add_compile_options(-Wall -Wextra -pedantic -Werror)

# Build
add_subdirectory(src)

##############################
# TESTING 
##############################
enable_testing()
add_subdirectory(test)