#-------------------------------------------------------------
# This script defines a variable JOBS to be used for parallel
# compilation in install and build cmake scripts.
#-------------------------------------------------------------

# Use this module to get the number of processors/cores available
include(ProcessorCount)

# Create variable NO_JOBS_DEFINED_PREVIOUSLY for later use
if(NOT DEFINED JOBS)
    set(NO_JOBS_DEFINED_PREVIOUSLY TRUE)
else()
    set(NO_JOBS_DEFINED_PREVIOUSLY FALSE)
endif()

# Get the maximum number of available processors/cores
ProcessorCount(MAX_JOBS)

# Ensure some parallel build for the external dependencies
if(NOT DEFINED JOBS OR JOBS LESS 1)
    if(MAX_JOBS LESS 4)
        set(JOBS 1)
    else()
        # Remove two cores/processors to prevent freezing of the system
        math(EXPR JOBS "${MAX_JOBS} - 2")
    endif()
endif()

# Issue a message to indicate how many parallel jobs are used.
if(NO_JOBS_DEFINED_PREVIOUSLY)
    message(STATUS "Using ${JOBS} parallel jobs for the compilation. This system permits a total of ${MAX_JOBS} parallel jobs.")
    message(STATUS "To change the default number of parallel jobs, use 'cmake -DJOBS=<number> -P <script-name>'.")
endif()
