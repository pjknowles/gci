LibraryManager_Append(${PROJECT_NAME}
        SOURCES SharedCounterGA.cpp SharedCounterMPI3.cpp
        PUBLIC_HEADER SharedCounter.h SharedCounterBase.h SharedCounterGA.h SharedCounterMPI3.h
        )

#TODO add option to use GA shared counter. This should not be the default just because GA was found.
#target_compile_definitions(${PROJECT_NAME} PUBLIC GCI_SCHEDULE_SharedCounterGA)
target_compile_definitions(${PROJECT_NAME} PUBLIC GCI_SCHEDULE_SharedCounterMPI3)
