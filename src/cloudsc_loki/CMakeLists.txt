# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Define this dwarf variant as an ECBuild feature
ecbuild_add_option( FEATURE CLOUDSC_LOKI
    DESCRIPTION "Use Loki source-to-source transformations with CLOUDSC " DEFAULT OFF
    CONDITION Serialbox_FOUND OR HDF5_FOUND
)

# Define the CLAW-based Loki dwarf variants as an ECBuild feature
ecbuild_add_option( FEATURE CLOUDSC_LOKI_CLAW
    DESCRIPTION "Build optimized CPU and GPU version of CLOUDSC derived from Loki+CLAW" DEFAULT OFF
    CONDITION HAVE_CLOUDSC_LOKI
)

if( HAVE_CLOUDSC_LOKI )

    ####################################################
    ##  Define various pre-processing modes via Loki  ##
    ####################################################

    set( COMMON_MODULE "${CMAKE_CURRENT_SOURCE_DIR}/../common/module" )
    set( COMMON_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/../common/include" )
    set( XMOD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/xmod" )

    set( LOKI_FRONTEND "fp" CACHE STRING "Frontend parser for Loki transforms" )

    # OFP frontend cannot deal with statement functions, so we toggle them here
    set( CLOUDSC_DEFINE_STMT_FUNC "" )
    if(NOT "${LOKI_FRONTEND}" STREQUAL "ofp")
    set( CLOUDSC_DEFINE_STMT_FUNC CLOUDSC_STMT_FUNC )
    endif()

    ####################################################
    ##  Idempotence mode:                             ##
    ##   * Internal "do-nothing" mode for Loki debug  ##
    ####################################################

    loki_transform(
        COMMAND convert
        OUTPUT
            loki-idem/cloudsc.idem.F90
            loki-idem/cloudsc_driver_loki_mod.idem.F90
        BUILDDIR ${CMAKE_CURRENT_BINARY_DIR}/loki-idem
        DEPENDS
            cloudsc.F90
            cloudsc_driver_loki_mod.F90
            ${_OMNI_DEPENDENCIES}
        MODE idem
        CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
        CPP
        DEFINITIONS
            ${CLOUDSC_DEFINE_STMT_FUNC}
        FRONTEND ${LOKI_FRONTEND}
        SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}
            ${COMMON_MODULE}
        INCLUDES
            ${COMMON_INCLUDE}
        XMOD
            ${_TARGET_XMOD_DIR}
            ${XMOD_DIR}
    )

    ecbuild_add_executable( TARGET dwarf-cloudsc-loki-idem
        SOURCES
            dwarf_cloudsc.F90
            loki-idem/cloudsc_driver_loki_mod.idem.F90
            loki-idem/cloudsc.idem.F90
        LIBS
            cloudsc-common-lib
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )

    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-idem-serial
        COMMAND bin/dwarf-cloudsc-loki-idem
        ARGS 1 100 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )
    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-idem-omp
        COMMAND bin/dwarf-cloudsc-loki-idem
        ARGS 4 100 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 4
        CONDITION HAVE_OMP
    )
    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-idem-mpi
        COMMAND bin/dwarf-cloudsc-loki-idem
        ARGS 1 100 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        MPI 2
        OMP 1
        CONDITION HAVE_MPI
    )
    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-idem-mpi-omp
        COMMAND bin/dwarf-cloudsc-loki-idem
        ARGS 4 100 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        MPI 2
        OMP 4
        CONDITION HAVE_OMP AND HAVE_MPI
    )

    ############################################################
    ##  Idempotence mode with pool-allocator for temporaries: ##
    ##   * Internal "do-nothing" mode for Loki debug          ##
    ############################################################

    loki_transform(
        COMMAND convert
        OUTPUT
            loki-idem-stack/cloudsc.idem_stack.F90
            loki-idem-stack/cloudsc_driver_loki_mod.idem_stack.F90
        BUILDDIR ${CMAKE_CURRENT_BINARY_DIR}/loki-idem-stack
        DEPENDS
            cloudsc.F90
            cloudsc_driver_loki_mod.F90
            ${_OMNI_DEPENDENCIES}
        MODE idem-stack
        CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
        CPP
        DEFINITIONS
            ${CLOUDSC_DEFINE_STMT_FUNC}
        FRONTEND ${LOKI_FRONTEND}
        SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}
            ${COMMON_MODULE}
        INCLUDES
            ${COMMON_INCLUDE}
        XMOD
            ${_TARGET_XMOD_DIR}
            ${XMOD_DIR}
    )

    ecbuild_add_executable( TARGET dwarf-cloudsc-loki-idem-stack
        SOURCES
            dwarf_cloudsc.F90
            loki-idem-stack/cloudsc_driver_loki_mod.idem_stack.F90
            loki-idem-stack/cloudsc.idem_stack.F90
        LIBS
            cloudsc-common-lib
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )

    if( CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
        target_compile_options( dwarf-cloudsc-loki-idem-stack PRIVATE "-fcray-pointer" )
    elseif( CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC" OR CMAKE_Fortran_COMPILER_ID MATCHES "PGI" )
        target_compile_options( dwarf-cloudsc-loki-idem-stack PRIVATE "-Mcray=pointer" )
    endif()

    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-idem-stack-serial
        COMMAND bin/dwarf-cloudsc-loki-idem-stack
        ARGS 1 100 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )
    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-idem-stack-omp
        COMMAND bin/dwarf-cloudsc-loki-idem-stack
        ARGS 4 100 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 4
        CONDITION HAVE_OMP
    )
    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-idem-stack-mpi
        COMMAND bin/dwarf-cloudsc-loki-idem-stack
        ARGS 1 100 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        MPI 2
        OMP 1
        CONDITION HAVE_MPI
    )
    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-idem-stack-mpi-omp
        COMMAND bin/dwarf-cloudsc-loki-idem-stack
        ARGS 4 100 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        MPI 2
        OMP 4
        CONDITION HAVE_OMP AND HAVE_MPI
    )


    ####################################################
    ##  SCA mode (Single Column Abstraction):         ##
    ##   * Extract de-vectorized SCA format code      ##
    ####################################################

    loki_transform(
        COMMAND convert
        OUTPUT
            loki-sca/cloudsc.sca.F90
            loki-sca/cloudsc_driver_loki_mod.sca.F90
        BUILDDIR ${CMAKE_CURRENT_BINARY_DIR}/loki-sca
        DEPENDS
            cloudsc.F90
            cloudsc_driver_loki_mod.F90
            ${_OMNI_DEPENDENCIES}
        MODE sca
        CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
        CPP
        DEFINITIONS
            ${CLOUDSC_DEFINE_STMT_FUNC}
        FRONTEND ${LOKI_FRONTEND}
        SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}
            ${COMMON_MODULE}
        INCLUDES
            ${COMMON_INCLUDE}
        XMOD
            ${_TARGET_XMOD_DIR}
            ${XMOD_DIR}
    )

    ecbuild_add_executable( TARGET dwarf-cloudsc-loki-sca
        SOURCES
            dwarf_cloudsc.F90
            loki-sca/cloudsc_driver_loki_mod.sca.F90
            loki-sca/cloudsc.sca.F90
        LIBS
            cloudsc-common-lib
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )

    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-sca-serial
        COMMAND bin/dwarf-cloudsc-loki-sca
        ARGS 1 100 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )
    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-sca-mpi
        COMMAND bin/dwarf-cloudsc-loki-sca
        ARGS 1 100 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        MPI 2
        OMP 1
        CONDITION HAVE_MPI
    )

    ####################################################
    ##  CLAW-CPU mode:                                ##
    ##   * Generate SCA code with CLAW annotations    ##
    ##   * Process with CLAW (CPU layout and OpenMP)  ##
    ####################################################
    if( HAVE_CLOUDSC_LOKI_CLAW AND TARGET clawfc )

        loki_transform(
            COMMAND convert
            OUTPUT
                loki-claw-cpu/cloudsc.claw.F90
                loki-claw-cpu/cloudsc_driver_loki_mod.claw.F90
            BUILDDIR ${CMAKE_CURRENT_BINARY_DIR}/loki-claw-cpu
            DEPENDS
                cloudsc.F90
                cloudsc_driver_loki_mod.F90
                ${_OMNI_DEPENDENCIES}
            MODE claw
            CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
            CPP
            DEFINITIONS
                ${CLOUDSC_DEFINE_STMT_FUNC}
            FRONTEND ${LOKI_FRONTEND}
            SOURCES
                ${CMAKE_CURRENT_SOURCE_DIR}
                ${COMMON_MODULE}
            INCLUDES
                ${COMMON_INCLUDE}
            XMOD
                ${_TARGET_XMOD_DIR}
                ${XMOD_DIR}
        )

        # We purposefully suppress CLAWs insertion of OpenMP loops,
        # as they are already present in the outer driver.
        claw_compile(
            MODEL_CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/claw_cloudsc.config
            TARGET cpu DIRECTIVE none
            INCLUDES ${COMMON_INCLUDE}
            XMOD ${_TARGET_XMOD_DIR} ${XMOD_DIR}
            SOURCE loki-claw-cpu/cloudsc.claw.F90
            OUTPUT loki-claw-cpu/cloudsc.claw.cpu.F90
        )
        claw_compile(
            MODEL_CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/claw_cloudsc.config
            TARGET cpu DIRECTIVE none
            INCLUDES ${COMMON_INCLUDE}
            XMOD ${_TARGET_XMOD_DIR} ${XMOD_DIR}
            SOURCE loki-claw-cpu/cloudsc_driver_loki_mod.claw.F90
            OUTPUT loki-claw-cpu/cloudsc_driver_loki_mod.claw.cpu.F90
            DEPENDS loki-claw-cpu/cloudsc.claw.cpu.F90
        )

        ecbuild_add_executable( TARGET dwarf-cloudsc-loki-claw-cpu
            SOURCES
                dwarf_cloudsc.F90
                loki-claw-cpu/cloudsc_driver_loki_mod.claw.cpu.F90
                loki-claw-cpu/cloudsc.claw.cpu.F90
            LIBS
                cloudsc-common-lib
            DEFINITIONS ${CLOUDSC_DEFINITIONS}
        )

        ecbuild_add_test(
            TARGET dwarf-cloudsc-loki-claw-cpu-serial
            COMMAND bin/dwarf-cloudsc-loki-claw-cpu
            ARGS 1 100 16
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
            OMP 1
        )
        ecbuild_add_test(
            TARGET dwarf-cloudsc-loki-claw-cpu-omp
            COMMAND bin/dwarf-cloudsc-loki-claw-cpu
            ARGS 4 100 16
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
            OMP 4
            CONDITION HAVE_OMP
        )
        ecbuild_add_test(
            TARGET dwarf-cloudsc-loki-claw-cpu-mpi
            COMMAND bin/dwarf-cloudsc-loki-claw-cpu
            ARGS 1 100 16
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
            MPI 2
            OMP 1
            CONDITION HAVE_MPI
        )
        ecbuild_add_test(
            TARGET dwarf-cloudsc-loki-claw-cpu-mpi-omp
            COMMAND bin/dwarf-cloudsc-loki-claw-cpu
            ARGS 4 100 16
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
            MPI 2
            OMP 4
            CONDITION HAVE_OMP AND HAVE_MPI
        )
    endif()

    ####################################################
    ##  CLAW-GPU mode:                                ##
    ##   * Generate SCA code with CLAW annotations    ##
    ##   * Process with CLAW (GPU layout and OpenACC) ##
    ####################################################
    if( HAVE_CLOUDSC_LOKI_CLAW AND TARGET clawfc )

        # Uses Loki-frontend CPP to switch to statement function variant again,
        # but suppresses inlining of stmt funcs by omitting `--include`

        loki_transform(
            COMMAND convert
            OUTPUT
                loki-claw-gpu/cloudsc.claw.F90
                loki-claw-gpu/cloudsc_driver_loki_mod.claw.F90
            BUILDDIR ${CMAKE_CURRENT_BINARY_DIR}/loki-claw-gpu
            DEPENDS
                cloudsc.F90
                cloudsc_driver_loki_mod.F90
                ${_OMNI_DEPENDENCIES}
            MODE claw
            CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
            CPP
            DATA_OFFLOAD
            REMOVE_OPENMP
            DEFINITIONS
                CLOUDSC_GPU_TIMING
                ${CLOUDSC_DEFINE_STMT_FUNC}
            FRONTEND ${LOKI_FRONTEND}
            SOURCES
                ${CMAKE_CURRENT_SOURCE_DIR}
                ${COMMON_MODULE}
            INCLUDES
                ${COMMON_INCLUDE}
            XMOD
                ${_TARGET_XMOD_DIR}
                ${XMOD_DIR}
        )

        claw_compile(
            MODEL_CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/claw_cloudsc.config
            TARGET gpu DIRECTIVE openacc
            INCLUDES ${COMMON_INCLUDE}
            XMOD ${_TARGET_XMOD_DIR} ${XMOD_DIR}
            SOURCE loki-claw-gpu/cloudsc.claw.F90
            OUTPUT loki-claw-gpu/cloudsc.claw.gpu.F90
        )
        claw_compile(
            MODEL_CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/claw_cloudsc.config
            TARGET gpu DIRECTIVE openacc
            INCLUDES ${COMMON_INCLUDE}
            XMOD ${_TARGET_XMOD_DIR} ${XMOD_DIR}
            SOURCE loki-claw-gpu/cloudsc_driver_loki_mod.claw.F90
            OUTPUT loki-claw-gpu/cloudsc_driver_loki_mod.claw.gpu.F90
            DEPENDS loki-claw-gpu/cloudsc.claw.gpu.F90
        )

        ecbuild_add_executable( TARGET dwarf-cloudsc-loki-claw-gpu
            SOURCES
                dwarf_cloudsc.F90
                loki-claw-gpu/cloudsc_driver_loki_mod.claw.gpu.F90
                loki-claw-gpu/cloudsc.claw.gpu.F90
            LIBS
                cloudsc-common-lib
            DEFINITIONS ${CLOUDSC_DEFINITIONS}
        )

        ecbuild_add_test(
            TARGET dwarf-cloudsc-loki-claw-gpu-serial
            COMMAND bin/dwarf-cloudsc-loki-claw-gpu
            ARGS 1 1280 128
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
            OMP 1
        )
    endif()


    ####################################################
    ##  "Single Column Coalesced" (SCC) mode          ##
    ##   * Removes horizontal vector loops            ##
    ##   * Invokes compute kernel as `!$acc vector`   ##
    ####################################################

    loki_transform(
        COMMAND convert
        OUTPUT
            loki-scc/cloudsc.scc.F90
            loki-scc/cloudsc_driver_loki_mod.scc.F90
        BUILDDIR ${CMAKE_CURRENT_BINARY_DIR}/loki-scc
        DEPENDS
            cloudsc.F90
            cloudsc_driver_loki_mod.F90
            ${_OMNI_DEPENDENCIES}
        MODE scc
        CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
        CPP
        DATA_OFFLOAD
        REMOVE_OPENMP
        DEFINITIONS
            CLOUDSC_GPU_TIMING
            ${CLOUDSC_DEFINE_STMT_FUNC}
        FRONTEND ${LOKI_FRONTEND}
        SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}
            ${COMMON_MODULE}
        INCLUDES
            ${COMMON_INCLUDE}
        XMOD
            ${_TARGET_XMOD_DIR}
            ${XMOD_DIR}
    )

    ecbuild_add_executable( TARGET dwarf-cloudsc-loki-scc
        SOURCES
            dwarf_cloudsc.F90
            loki-scc/cloudsc_driver_loki_mod.scc.F90
            loki-scc/cloudsc.scc.F90
        LIBS
            cloudsc-common-lib
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )

    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-scc-serial
        COMMAND bin/dwarf-cloudsc-loki-scc
        ARGS 1 1280 128
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
        ENVIRONMENT "NVCOMPILER_ACC_CUDA_HEAPSIZE=128M"
    )

    ######################################################
    ##  "Single Column Coalesced" (SCC) mode with stack ##
    ##   * Removes horizontal vector loops              ##
    ##   * Invokes compute kernel as `!$acc vector`     ##
    ##   * Allocates temporaries using pool allocator   ##
    ######################################################

    loki_transform(
        COMMAND convert
        OUTPUT
            loki-scc-stack/cloudsc.scc_stack.F90
            loki-scc-stack/cloudsc_driver_loki_mod.scc_stack.F90
        BUILDDIR ${CMAKE_CURRENT_BINARY_DIR}/loki-scc-stack
        DEPENDS
            cloudsc.F90
            cloudsc_driver_loki_mod.F90
            ${_OMNI_DEPENDENCIES}
        MODE scc-stack
        CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
        CPP
        DATA_OFFLOAD
        REMOVE_OPENMP
        DEFINITIONS
            CLOUDSC_GPU_TIMING
            ${CLOUDSC_DEFINE_STMT_FUNC}
        FRONTEND ${LOKI_FRONTEND}
        SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}
            ${COMMON_MODULE}
        INCLUDES
            ${COMMON_INCLUDE}
        XMOD
            ${_TARGET_XMOD_DIR}
            ${XMOD_DIR}
    )

    ecbuild_add_executable( TARGET dwarf-cloudsc-loki-scc-stack
        SOURCES
            dwarf_cloudsc.F90
            loki-scc-stack/cloudsc_driver_loki_mod.scc_stack.F90
            loki-scc-stack/cloudsc.scc_stack.F90
        LIBS
            cloudsc-common-lib
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )

    if( CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
        target_compile_options( dwarf-cloudsc-loki-scc-stack PRIVATE "-fcray-pointer" )
    elseif( CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC" OR CMAKE_Fortran_COMPILER_ID MATCHES "PGI" )
        target_compile_options( dwarf-cloudsc-loki-scc-stack PRIVATE "-Mcray=pointer" )
    endif()

    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-scc-stack-serial
        COMMAND bin/dwarf-cloudsc-loki-scc-stack
        ARGS 1 1280 128
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )


    ####################################################
    ##  SCC-hoist mode                                ##
    ##   * SCC with vector loop hoisted               ##
    ##   * Kernel is "seq, but args are full blocks   ##
    ##   * Temporary arrays hoisted to driver         ##
    ####################################################

    loki_transform(
        COMMAND convert
        OUTPUT
            loki-scc-hoist/cloudsc.scc_hoist.F90
            loki-scc-hoist/cloudsc_driver_loki_mod.scc_hoist.F90
        BUILDDIR ${CMAKE_CURRENT_BINARY_DIR}/loki-scc-hoist
        DEPENDS
            cloudsc.F90
            cloudsc_driver_loki_mod.F90
            ${_OMNI_DEPENDENCIES}
        MODE scc-hoist
        CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
        CPP
        DATA_OFFLOAD
        REMOVE_OPENMP
        DEFINITIONS
            CLOUDSC_GPU_TIMING
            ${CLOUDSC_DEFINE_STMT_FUNC}
        FRONTEND ${LOKI_FRONTEND}
        SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}
            ${COMMON_MODULE}
        INCLUDES
            ${COMMON_INCLUDE}
        XMOD
            ${_TARGET_XMOD_DIR}
            ${XMOD_DIR}
    )

    ecbuild_add_executable( TARGET dwarf-cloudsc-loki-scc-hoist
        SOURCES
            dwarf_cloudsc.F90
            loki-scc-hoist/cloudsc_driver_loki_mod.scc_hoist.F90
            loki-scc-hoist/cloudsc.scc_hoist.F90
        LIBS
            cloudsc-common-lib
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )

    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-scc-hoist-serial
        COMMAND bin/dwarf-cloudsc-loki-scc-hoist
        ARGS 1 1280 128
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )


    ####################################################
    ##  SCC CUF (CUDA Fortran)                        ##
    ##   * SCC with CUDA Fortran (CUF)                ##
    ##   * allows for 3 types/kind of transformation  ##
    ##    * parametrised (dimensions) local arrays    ##
    ##    * host side hoisted local arrays            ##
    ##    * dynamic memory allocation for local       ##
    ##      arrays (not recommended)                  ##
    ####################################################


if( HAVE_CUDA )

    # scc-cuf-parametrise
    loki_transform_convert(
        MODE cuf-parametrise FRONTEND ${LOKI_FRONTEND} CPP
        CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_cuf_loki.config
        PATH ${CMAKE_CURRENT_SOURCE_DIR}
        HEADERS
            ${COMMON_MODULE}/yomcst.F90
            ${COMMON_MODULE}/yomphyder.F90
            ${COMMON_MODULE}/yoethf.F90
            ${COMMON_MODULE}/yoecldp.F90
        INCLUDES ${COMMON_INCLUDE}
        DEFINITIONS CLOUDSC_GPU_TIMING
        DATA_OFFLOAD REMOVE_OPENMP
        XMOD ${_TARGET_XMOD_DIR} ${XMOD_DIR}
        OUTPATH ${CMAKE_CURRENT_BINARY_DIR}/loki-scc-cuf-parametrise
        OUTPUT
            loki-scc-cuf-parametrise/cuf_cloudsc_driver_loki_mod.cuf_parametrise.F90
            loki-scc-cuf-parametrise/cuf_cloudsc.cuf_parametrise.F90
        DEPENDS cuf_cloudsc.F90 cuf_cloudsc_driver_loki_mod.F90 ${_OMNI_DEPENDENCIES}
    )

    set_source_files_properties(
        loki-scc-cuf-parametrise/cuf_cloudsc_driver_loki_mod.cuf_parametrise.F90
        loki-scc-cuf-parametrise/cuf_cloudsc.cuf_parametrise.F90
        PROPERTIES COMPILE_FLAGS "-cuda -gpu=maxregcount:128"
    )

    ecbuild_add_executable( TARGET dwarf-cloudsc-loki-scc-cuf-parametrise
        SOURCES
            dwarf_cloudsc.F90
            loki-scc-cuf-parametrise/cuf_cloudsc_driver_loki_mod.cuf_parametrise.F90
            loki-scc-cuf-parametrise/cuf_cloudsc.cuf_parametrise.F90
        LIBS
            cloudsc-common-lib
        DEFINITIONS ${CLOUDSC_DEFINITIONS} CLOUDSC_GPU_SCC_CUF
    )

    target_link_options(dwarf-cloudsc-loki-scc-cuf-parametrise PUBLIC "-cuda")

    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-scc-cuf-parametrise-serial
        COMMAND bin/dwarf-cloudsc-loki-scc-cuf-parametrise
        ARGS 1 1280 128
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )

    # scc-cuf-hoist
    loki_transform_convert(
        MODE cuf-hoist FRONTEND ${LOKI_FRONTEND} CPP
        CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_cuf_loki.config
        PATH ${CMAKE_CURRENT_SOURCE_DIR}
        HEADERS
            ${COMMON_MODULE}/yomcst.F90
            ${COMMON_MODULE}/yomphyder.F90
            ${COMMON_MODULE}/yoethf.F90
            ${COMMON_MODULE}/yoecldp.F90
        INCLUDES ${COMMON_INCLUDE}
        DEFINITIONS CLOUDSC_GPU_TIMING
        DATA_OFFLOAD REMOVE_OPENMP
        XMOD ${_TARGET_XMOD_DIR} ${XMOD_DIR}
        OUTPATH ${CMAKE_CURRENT_BINARY_DIR}/loki-scc-cuf-hoist
        OUTPUT
            loki-scc-cuf-hoist/cuf_cloudsc_driver_loki_mod.cuf_hoist.F90
            loki-scc-cuf-hoist/cuf_cloudsc.cuf_hoist.F90
    DEPENDS cuf_cloudsc.F90 cuf_cloudsc_driver_loki_mod.F90 ${_OMNI_DEPENDENCIES}
    )

    set_source_files_properties(
        loki-scc-cuf-hoist/cuf_cloudsc_driver_loki_mod.cuf_hoist.F90
        loki-scc-cuf-hoist/cuf_cloudsc.cuf_hoist.F90
        PROPERTIES COMPILE_FLAGS "-cuda -gpu=maxregcount:128"
    )

    ecbuild_add_executable( TARGET dwarf-cloudsc-loki-scc-cuf-hoist
        SOURCES
            dwarf_cloudsc.F90
            loki-scc-cuf-hoist/cuf_cloudsc_driver_loki_mod.cuf_hoist.F90
            loki-scc-cuf-hoist/cuf_cloudsc.cuf_hoist.F90
        LIBS
            cloudsc-common-lib
        DEFINITIONS ${CLOUDSC_DEFINITIONS} CLOUDSC_GPU_SCC_CUF
    )

    target_link_options(dwarf-cloudsc-loki-scc-cuf-hoist PUBLIC "-cuda")

    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-scc-cuf-hoist-serial
        COMMAND bin/dwarf-cloudsc-loki-scc-cuf-hoist
        ARGS 1 1280 128
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )

endif()

    #############################################################################
    ## C-transpilation mode for generating vectorized C host code (experimental!)
    #############################################################################

    if ( NOT HAVE_SINGLE_PRECISION )

        loki_transform(
            COMMAND convert
            MODE c FRONTEND ${LOKI_FRONTEND} CPP
            CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
            SOURCES ${CMAKE_CURRENT_SOURCE_DIR}
            HEADERS
                ${COMMON_MODULE}/parkind1.F90
                ${COMMON_MODULE}/yomphyder.F90
                ${COMMON_MODULE}/yomcst.F90
                ${COMMON_MODULE}/yoethf.F90
                ${COMMON_MODULE}/yoecldp.F90
                ${COMMON_MODULE}/fcttre_mod.F90
                ${COMMON_MODULE}/fccld_mod.F90
            INCLUDES ${COMMON_INCLUDE}
            XMOD ${_TARGET_XMOD_DIR} ${XMOD_DIR}
            BUILDDIR ${CMAKE_CURRENT_BINARY_DIR}/loki-c
            OUTPUT
                loki-c/cloudsc_driver_loki_mod.c.F90
                loki-c/cloudsc_fc.F90 loki-c/cloudsc_c.c
                loki-c/yoethf_fc.F90 loki-c/yomcst_fc.F90
                loki-c/yoecldp_fc.F90
            DEPENDS cloudsc.F90 cloudsc_driver_loki_mod.F90 ${_OMNI_DEPENDENCIES}
        )

        # Define the CLAW-CPU build target for this variant
        ecbuild_add_executable( TARGET dwarf-cloudsc-loki-c
            SOURCES
                dwarf_cloudsc.F90
                loki-c/cloudsc_driver_loki_mod.c.F90
                loki-c/cloudsc_fc.F90
                loki-c/cloudsc_c.c
                loki-c/yoethf_fc.F90
                loki-c/yomcst_fc.F90
                loki-c/yoecldp_fc.F90
            LIBS
                cloudsc-common-lib
            DEFINITIONS ${CLOUDSC_DEFINITIONS}
        )

        ecbuild_add_test(
            TARGET dwarf-cloudsc-loki-c-serial
            COMMAND bin/dwarf-cloudsc-loki-c
            ARGS 1 100 16
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
            OMP 1
        )
        ecbuild_add_test(
            TARGET dwarf-cloudsc-loki-c-omp
            COMMAND bin/dwarf-cloudsc-loki-c
            ARGS 4 100 16
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
            OMP 4
            CONDITION HAVE_OMP
        )

    endif()

    # Create symlink for the input data
    if( HAVE_SERIALBOX )
      execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${CMAKE_CURRENT_SOURCE_DIR}/../../data ${CMAKE_CURRENT_BINARY_DIR}/../../../data )
    endif()

    if( HAVE_HDF5 )
      execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${CMAKE_CURRENT_SOURCE_DIR}/../../config-files/input.h5 ${CMAKE_CURRENT_BINARY_DIR}/../../../input.h5 )
      execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${CMAKE_CURRENT_SOURCE_DIR}/../../config-files/reference.h5 ${CMAKE_CURRENT_BINARY_DIR}/../../../reference.h5 )
    endif()
endif()
