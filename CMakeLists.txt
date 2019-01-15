# Minimum required Version
cmake_minimum_required(VERSION 2.6 FATAL_ERROR)

# Projects name
project(spin_hamiltonian)

option(DEBUG "Enable debugging" OFF)

enable_language(Fortran)

set(SOURCE_SPIN ${SOURCE_SPIN}  source/main.F90
                                source/mod_tools.F90
                                source/mod_types.F90
                                source/mod_parameters.F90
                                source/mod_constants.F90
                                source/mod_spins.F90
                                source/mod_hamiltonian.F90
                                source/eigensolver.F90
                                )

# Execulable names
set(SPIN_EXE sh.exe)
set(CMAKE_Fortran_COMPILER ifort)

add_executable(${SPIN_EXE} ${SOURCE_SPIN})

set(CMAKE_BINARY_DIR ${CMAKE_SOURCE_DIR}/bin)
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR})
set(LIBRARY_OUTPUT_PATH ${CMAKE_BINARY_DIR})

execute_process(COMMAND git describe --abbrev=4 --dirty --always --tags OUTPUT_VARIABLE GIT_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)

target_link_libraries(${SPIN_EXE} -mkl)
target_compile_options(${SPIN_EXE} PRIVATE -no-wrap-margin
                                            -fpp
                                            -DVERSION="${GIT_VERSION}"
                                            -mkl
                                            -xHost)


if(DEBUG)
  target_compile_options(${SPIN_EXE} PRIVATE  -C
                                              -CB
                                              -check uninit
                                              -debug all
                                              -warn all
                                              -ftrapuv
                                              -traceback
                                              -fpp
                                              -DVERSION="${GIT_VERSION}"
                                              -g
                                              -O0)
else()
  target_compile_options(${SPIN_EXE} PRIVATE -O3)
endif(DEBUG)
