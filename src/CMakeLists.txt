add_library(solvers
    csc.h
    csc.c

    kkt.h
    kkt.c

    qdldl_solver.h
    qdldl_solver.c

    utils.h
    utils.c
)
target_link_libraries(solvers
    PUBLIC
    QDLDL::QDLDL

    PRIVATE
    CJSON::CJSON
)
target_include_directories(solvers
    PUBLIC
    ${CMAKE_CURRENT_LIST_DIR}
)

# QDLDL Executable 
add_executable(solve_qdldl
    solve_qdldl.c
)
target_link_libraries(solve_qdldl
    PRIVATE
    solvers
)

if (PARDISO_FOUND)
    message(STATUS "Building Pardiso targets")

    add_library(pardiso_solver
        pardiso_solver.h
        pardiso_solver.c
    )
    target_link_libraries(pardiso_solver
        PUBLIC
        PARDISO::PARDISO
    )

    # PARDISO Executable
    add_executable(solve_pardiso
        solve_pardiso.c
    )
    target_link_libraries(solve_pardiso
        PRIVATE
        solvers
        pardiso_solver
    )

    # Sample from PARDISO website
    add_executable(pardiso_example
        pardiso_sym.c
    )
    target_link_libraries(pardiso_example
        PRIVATE
        PARDISO::PARDISO
    )
endif()
