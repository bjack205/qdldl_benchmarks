add_library(csc
    csc.h
    csc.c
)
target_link_libraries(csc
    PRIVATE
    CJSON::CJSON
)