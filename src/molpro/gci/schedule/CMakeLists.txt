LibraryManager_Append(${PROJECT_NAME}
        SOURCES SharedCounterGA.cpp SharedCounterMPI3.cpp
        PUBLIC_HEADER SharedCounter.h SharedCounterBase.h SharedCounterGA.h SharedCounterMPI3.h
        )

if (GCI_SHAREDCOUNTER_MPI3)
    message(STATUS "Building gci with MPI3 SharedCounter")
    target_compile_definitions(${PROJECT_NAME} PUBLIC GCI_SCHEDULE_SharedCounterMPI3)
elseif (GCI_SHAREDCOUNTER_GA)
    message(STATUS "Building gci with GA SharedCounter")
    target_compile_definitions(${PROJECT_NAME} PUBLIC GCI_SCHEDULE_SharedCounterGA)
else ()
    message(FATAL_ERROR "SharedCounter not specified, turn on one of the options: GCI_SHAREDCOUNTER_MPI3, GCI_SHAREDCOUNTER_GA")
endif ()
