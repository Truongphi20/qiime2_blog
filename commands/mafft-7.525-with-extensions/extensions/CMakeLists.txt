
set(MXSCARNA_SRC_DIR mxscarna_src)
set(PROBCONS_DIR ${MXSCARNA_SRC_DIR}/probconsRNA)
set(VIENNA_DIR ${MXSCARNA_SRC_DIR}/vienna)

add_executable(mxscarnamod
    ${MXSCARNA_SRC_DIR}/Main.cc
    ${MXSCARNA_SRC_DIR}/McCaskill.cpp
    ${MXSCARNA_SRC_DIR}/vienna/energy_param.cpp
    ${MXSCARNA_SRC_DIR}/seq2scs.cpp
    ${MXSCARNA_SRC_DIR}/Globaldp.cpp
    ${MXSCARNA_SRC_DIR}/postProcessings.cpp
    ${MXSCARNA_SRC_DIR}/AlifoldMEA.cpp
)

target_include_directories(mxscarnamod PRIVATE
    ${MXSCARNA_SRC_DIR}
    ${PROBCONS_DIR}
    ${VIENNA_DIR}
)

target_compile_definitions(mxscarnamod PRIVATE
    NumInsertStates=1
    VERSION="2.0"
    NDEBUG
)

# C++98 is requested in the original Makefile
set_target_properties(mxscarnamod PROPERTIES
    CXX_STANDARD 98
    CXX_STANDARD_REQUIRED ON
)

target_link_libraries(mxscarnamod PRIVATE ${MATH_LIBRARY})

install(TARGETS mxscarnamod DESTINATION ${MAFFT_LIBEXECDIR})
