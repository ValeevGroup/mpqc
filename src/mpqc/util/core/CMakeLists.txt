set(sources
        backtrace.cpp
        backtrace.h
        exenv.cpp
        exenv.h
        exception.cpp
        exception.h
        formio.cpp
        formio.h
        )

add_mpqc_library(util_core sources sources "MADworld" "mpqc/util/core")
