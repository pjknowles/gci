configure_file(${PROJECT_SOURCE_DIR}/gci.fcidump ${PROJECT_BINARY_DIR}/gci.fcidump COPYONLY)
set(TESTS he2.fcidump he2n.fcidump hedz.fcidump hevdz.fcidump nevdz.fcidump)
foreach (fil ${TESTS})
    get_filename_component(nam ${fil} NAME_WE)
    configure_file(${fil} ${CMAKE_CURRENT_BINARY_DIR}/${file} COPYONLY)
    add_test(NAME ${nam}
            COMMAND ${MPIEXEC_EXECUTABLE} ${MPIOPTIONS} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS}
            ${PROJECT_BINARY_DIR}/gci.exe ${MPIEXEC_POSTFLAGS} ${fil}
            )
endforeach ()

add_subdirectory(unit_test)
