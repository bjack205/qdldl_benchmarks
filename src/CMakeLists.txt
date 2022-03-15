add_library(csc
    csc.h
    csc.c
)
target_link_libraries(csc
    PRIVATE
    CJSON::CJSON
)
target_include_directories(csc
    PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
)

add_library(solvers
    qdldl_solver.h
    qdldl_solver.c

    pardiso_solver.h
    pardiso_solver.c

    utils.h
    utils.c
)
target_link_libraries(solvers
    PUBLIC
    QDLDL::QDLDL
    csc
    PARDISO::PARDISO
)
target_include_directories(solvers
    PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
)

add_executable(solve_qdldl
    solve_qdldl.c
)
target_link_libraries(solve_qdldl
    PRIVATE
    solvers
)

add_executable(solve_pardiso
    solve_pardiso.c
)
target_link_libraries(solve_pardiso
    PRIVATE
    solvers
)

add_executable(pardiso_example
    pardiso_sym.c
)
target_link_libraries(pardiso_example
    PRIVATE
    PARDISO::PARDISO
)