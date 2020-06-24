LibraryManager_Append(${PROJECT_NAME}
        SOURCES SharedCounterGA.cpp SharedCounterMPI3.cpp
        PUBLIC_HEADER SharedCounterBase.h SharedCounterGA.h SharedCounterMPI3.h
        )

target_compile_definitions(${PROJECT_NAME} PRIVATE GCI_SCHEDULE_SHAREDCOUNTER_TYPE=SharedCounterBasic)