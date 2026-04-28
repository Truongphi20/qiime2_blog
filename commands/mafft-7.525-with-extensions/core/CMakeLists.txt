
set(COMMON_SOURCES
    mtxutl.c io.c defs.c mltaln9.c Galign11.c Lalign11.c genalign11.c
)

set(FFT_SOURCES
    fft.c fftFunctions.c
)

set(ALIGN_SOURCES
    tddis.c constants.c partSalignmm.c Lalignmm.c rna.c Salignmm.c Dalignmm.c
    Falign.c Falign_localhom.c SAalignmm.c MSalignmm.c
)

macro(add_mafft_executable name)
    add_executable(${name} ${ARGN})
    target_link_libraries(${name} PRIVATE ${MATH_LIBRARY} ${THREADS_LIB})
    install(TARGETS ${name} DESTINATION ${MAFFT_LIBEXECDIR})
endmacro()

add_mafft_executable(dvtditr ${COMMON_SOURCES} ${FFT_SOURCES} ${ALIGN_SOURCES} nj.c tditeration.c dvtditr.c treeOperation.c addfunctions.c)
add_mafft_executable(dndfast7 ${COMMON_SOURCES} tddis.c constants.c dndfast7.c)
add_mafft_executable(dndblast ${COMMON_SOURCES} tddis.c constants.c dndblast.c)
add_mafft_executable(sextet5 ${COMMON_SOURCES} tddis.c constants.c sextet5.c)
add_mafft_executable(mafft-distance ${COMMON_SOURCES} tddis.c constants.c mafft-distance.c)
add_mafft_executable(pairlocalalign ${COMMON_SOURCES} ${FFT_SOURCES} ${ALIGN_SOURCES} MSalign11.c suboptalign11.c pairlocalalignmain.c pairlocalalign.c)
add_mafft_executable(multi2hat3s ${COMMON_SOURCES} ${FFT_SOURCES} tddis.c constants.c multi2hat3s.c)
add_mafft_executable(pairash ${COMMON_SOURCES} ${FFT_SOURCES} ${ALIGN_SOURCES} MSalign11.c suboptalign11.c pairash.c)
add_mafft_executable(addsingle ${COMMON_SOURCES} ${FFT_SOURCES} ${ALIGN_SOURCES} addsingle.c addfunctions.c)
add_mafft_executable(maffttext2hex maffttext2hex.c)
add_mafft_executable(hex2maffttext hex2maffttext.c)
add_mafft_executable(splittbfast ${COMMON_SOURCES} ${FFT_SOURCES} ${ALIGN_SOURCES} splittbfast.c)
add_mafft_executable(disttbfast ${COMMON_SOURCES} ${FFT_SOURCES} ${ALIGN_SOURCES} disttbfast.c addfunctions.c)
add_mafft_executable(tbfast ${COMMON_SOURCES} ${FFT_SOURCES} ${ALIGN_SOURCES} tbfast.c addfunctions.c pairlocalalign.c MSalign11.c)
add_mafft_executable(nodepair ${COMMON_SOURCES} ${FFT_SOURCES} ${ALIGN_SOURCES} nodepair.c addfunctions.c pairlocalalign.c MSalign11.c)
add_mafft_executable(mafft-profile ${COMMON_SOURCES} ${FFT_SOURCES} ${ALIGN_SOURCES} mafft-profile.c)
add_mafft_executable(f2cl ${COMMON_SOURCES} constants.c f2cl.c)
add_mafft_executable(mccaskillwrap ${COMMON_SOURCES} constants.c mccaskillwrap.c)
add_mafft_executable(contrafoldwrap ${COMMON_SOURCES} constants.c contrafoldwrap.c)
add_mafft_executable(countlen ${COMMON_SOURCES} countlen.c)
add_mafft_executable(seq2regtable io.c mtxutl.c defs.c seq2regtable.c)
add_mafft_executable(regtable2seq ${COMMON_SOURCES} regtable2seq.c)
add_mafft_executable(score ${COMMON_SOURCES} constants.c score.c)
add_mafft_executable(getlag ${COMMON_SOURCES} ${FFT_SOURCES} ${ALIGN_SOURCES} getlag.c)
add_mafft_executable(dndpre ${COMMON_SOURCES} constants.c dndpre.c)
add_mafft_executable(setcore ${COMMON_SOURCES} ${FFT_SOURCES} ${ALIGN_SOURCES} setcore.c)
add_mafft_executable(filter ${COMMON_SOURCES} filter.c)
add_mafft_executable(replaceu ${COMMON_SOURCES} replaceu.c)
add_mafft_executable(restoreu ${COMMON_SOURCES} restoreu.c)
add_mafft_executable(setdirection ${COMMON_SOURCES} setdirection.c)
add_mafft_executable(makedirectionlist ${COMMON_SOURCES} ${FFT_SOURCES} ${ALIGN_SOURCES} makedirectionlist.c addfunctions.c)
add_mafft_executable(version version.c)

# dash_client
if(ENABLE_DASH_CLIENT)
    find_program(GO_EXECUTABLE go)
    if(GO_EXECUTABLE)
        add_custom_command(
            OUTPUT dash_client
            COMMAND ${CMAKE_COMMAND} -E env CGO_ENABLED=0 ${GO_EXECUTABLE} build --ldflags "-extldflags '-static'" -o dash_client ${CMAKE_CURRENT_SOURCE_DIR}/dash_client.go
            DEPENDS dash_client.go
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
        )
        add_custom_target(dash_client_target ALL DEPENDS dash_client)
        install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/dash_client DESTINATION ${MAFFT_LIBEXECDIR})
    else()
        message(WARNING "Go executable not found. dash_client will not be built.")
    endif()
endif()

# Scripts
set(LIBDIR "${CMAKE_INSTALL_FULL_LIBEXECDIR}/mafft")
set(BINDIR "${CMAKE_INSTALL_FULL_BINDIR}")

macro(configure_mafft_script template_name output_name placeholder value)
    file(READ ${template_name} TEMPLATE_CONTENT)
    string(REPLACE "${placeholder}" "${value}" SCRIPT_CONTENT "${TEMPLATE_CONTENT}")
    file(WRITE ${CMAKE_CURRENT_BINARY_DIR}/${output_name} "${SCRIPT_CONTENT}")
    install(PROGRAMS ${CMAKE_CURRENT_BINARY_DIR}/${output_name} DESTINATION ${CMAKE_INSTALL_BINDIR})
endmacro()

configure_mafft_script(mafft.tmpl mafft _LIBDIR ${LIBDIR})
configure_mafft_script(mafft-homologs.tmpl mafft-homologs.rb _BINDIR ${BINDIR})
configure_mafft_script(mafft-sparsecore.tmpl mafft-sparsecore.rb _BINDIR ${BINDIR})

# Perl scripts
install(PROGRAMS mafftash_premafft.tmpl DESTINATION ${MAFFT_LIBEXECDIR} RENAME mafftash_premafft.pl)
install(PROGRAMS seekquencer_premafft.tmpl DESTINATION ${MAFFT_LIBEXECDIR} RENAME seekquencer_premafft.pl)

# Symlinks
set(MAFFT_LINKS
    linsi ginsi einsi fftns fftnsi nwns nwnsi xinsi qinsi
    mafft-linsi mafft-ginsi mafft-einsi mafft-fftns mafft-fftnsi mafft-nwns mafft-nwnsi mafft-xinsi mafft-qinsi
)

foreach(link ${MAFFT_LINKS})
    install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink mafft ${link} WORKING_DIRECTORY \$ENV{DESTDIR}${CMAKE_INSTALL_FULL_BINDIR})")
endforeach()

install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ../${MAFFT_LIBEXECDIR}/mafft-profile mafft-profile WORKING_DIRECTORY \$ENV{DESTDIR}${CMAKE_INSTALL_FULL_BINDIR})")
install(CODE "execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink ../${MAFFT_LIBEXECDIR}/mafft-distance mafft-distance WORKING_DIRECTORY \$ENV{DESTDIR}${CMAKE_INSTALL_FULL_BINDIR})")

# Man pages
install(FILES mafft.1 mafft-homologs.1 DESTINATION ${CMAKE_INSTALL_MANDIR}/man1)
