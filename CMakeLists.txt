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

# find_pardiso()
# # Find LAPACK
# find_package(LAPACK REQUIRED)
# cmake_print_variables(LAPACK_FOUND)

# # Find BLAS
# find_package(BLAS REQUIRED)
# cmake_print_variables(BLAS_FOUND)

# # Find OpenMP
# find_package(OpenMP REQUIRED)
# cmake_print_variables(OpenMP_FOUND)

# # Find Pardiso
# SET(FIND_LIBPARDISO LIBPARDISO-NOTFOUND)
# find_library(FIND_LIBPARDISO ${PARDISO_LIB} 
#     HINTS ${PARDISO_DIR} 
#     PATH_SUFFIXES "pardiso-6.0"
# )
# cmake_print_variables(FIND_LIBPARDISO)
# add_library(PARDISO SHARED IMPORTED)
# set_target_properties(PARDISO PROPERTIES
#     IMPORTED_LOCATION ${FIND_LIBPARDISO}
# )

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