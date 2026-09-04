set(MAXWELL_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")

file(GLOB_RECURSE MAXWELL_CORE_SOURCES CONFIGURE_DEPENDS
    "${MAXWELL_SOURCE_DIR}/classes/grid/*.cu"
    "${MAXWELL_SOURCE_DIR}/classes/source/*.cu"
    "${MAXWELL_SOURCE_DIR}/classes/pml/*.cu"
    "${MAXWELL_SOURCE_DIR}/classes/write-output/*.cu"
    "${MAXWELL_SOURCE_DIR}/classes/validation/*.cu"
)
kestrel_configure_sources(${MAXWELL_CORE_SOURCES})

add_library(kestrel-maxwell-core STATIC ${MAXWELL_CORE_SOURCES})
target_include_directories(kestrel-maxwell-core PUBLIC "${MAXWELL_SOURCE_DIR}")
target_link_libraries(kestrel-maxwell-core PUBLIC xpu::xpu)

set(MAXWELL_OPT_FLAGS
    -O3
    -march=native
    -mtune=native
    -funroll-loops
    -ffast-math
    -flto
    -DNDEBUG
)

set(MAXWELL_WARN_FLAGS
    -Wall
    -Wextra
    -Wpedantic
    -Wshadow
    -Wconversion
)

if(IS_DEBUG)
    if(KESTREL_CUDA)
        target_compile_options(kestrel-maxwell-core PRIVATE
            -O0 -g
            -Xcompiler=-Wall,-Wextra,-Wshadow,-Wconversion,-fno-omit-frame-pointer
        )
    else()
        target_compile_options(kestrel-maxwell-core PRIVATE
            -O0 -g3
            ${MAXWELL_WARN_FLAGS}
            -fno-omit-frame-pointer
            -fsanitize=address,undefined
        )
        target_link_options(kestrel-maxwell-core PUBLIC -fsanitize=address,undefined)
    endif()
elseif(IS_TEST)
    if(KESTREL_CUDA)
        target_compile_options(kestrel-maxwell-core PRIVATE
            -O3
            -Xcompiler=-O3,-march=native,-mtune=native,-funroll-loops,-ffast-math,-DNDEBUG
            -Xcompiler=-Wall,-Wextra,-Wshadow,-Wconversion
        )
    else()
        target_compile_options(kestrel-maxwell-core PRIVATE
            ${MAXWELL_OPT_FLAGS}
            ${MAXWELL_WARN_FLAGS}
        )
        target_link_options(kestrel-maxwell-core PUBLIC -flto)
    endif()
elseif(IS_RELEASE)
    if(KESTREL_CUDA)
        target_compile_options(kestrel-maxwell-core PRIVATE
            -O3
            -Xcompiler=-O3,-march=native,-mtune=native,-funroll-loops,-ffast-math,-DNDEBUG
        )
    else()
        target_compile_options(kestrel-maxwell-core PRIVATE ${MAXWELL_OPT_FLAGS})
        target_link_options(kestrel-maxwell-core PUBLIC -flto)
    endif()
endif()

find_package(OpenMP REQUIRED COMPONENTS CXX)
target_link_libraries(kestrel-maxwell-core PUBLIC OpenMP::OpenMP_CXX)
target_compile_options(kestrel-maxwell-core PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:-Xcompiler=-fopenmp>
)

if(ipo_supported AND NOT IS_DEBUG AND NOT KESTREL_CUDA)
    set_property(TARGET kestrel-maxwell-core PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

# Maxwell example executable:

file(GLOB_RECURSE MAXWELL_SIMULATION_SOURCES CONFIGURE_DEPENDS
    "${MAXWELL_SOURCE_DIR}/classes/simulation/*.cu"
)
kestrel_configure_sources(
    "${MAXWELL_SOURCE_DIR}/main.cu"
    ${MAXWELL_SIMULATION_SOURCES}
)

add_executable(
    kestrel-maxwell
    "${MAXWELL_SOURCE_DIR}/main.cu"
    ${MAXWELL_SIMULATION_SOURCES}
)
target_link_libraries(kestrel-maxwell PRIVATE kestrel-maxwell-core)

if(ipo_supported AND NOT IS_DEBUG AND NOT KESTREL_CUDA)
    set_property(TARGET kestrel-maxwell PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()

configure_file(config.cfg config.cfg COPYONLY)

# Maxwell tests are kept intact until the generic Grid and Field contracts exist.

if(BUILD_TESTS)
    file(GLOB MAXWELL_TEST_SOURCES CONFIGURE_DEPENDS
        "${CMAKE_CURRENT_SOURCE_DIR}/tests/test-*.cu"
    )
    kestrel_configure_sources(${MAXWELL_TEST_SOURCES})

    add_executable(kestrel-maxwell-tests ${MAXWELL_TEST_SOURCES})
    target_include_directories(
        kestrel-maxwell-tests
        PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/tests"
    )
    target_link_libraries(kestrel-maxwell-tests PRIVATE kestrel-maxwell-core)

    if(IS_DEBUG)
        if(KESTREL_CUDA)
            target_compile_options(kestrel-maxwell-tests PRIVATE
                -O1 -g
                -Xcompiler=-Wall,-Wextra,-Wshadow,-Wconversion,-fno-omit-frame-pointer
            )
        else()
            target_compile_options(kestrel-maxwell-tests PRIVATE
                -O1 -g3
                ${MAXWELL_WARN_FLAGS}
                -fno-omit-frame-pointer
                -fsanitize=address,undefined
            )
            target_link_options(kestrel-maxwell-tests PRIVATE -fsanitize=address,undefined)
        endif()
    elseif(IS_TEST)
        if(KESTREL_CUDA)
            target_compile_options(kestrel-maxwell-tests PRIVATE
                -O2
                -Xcompiler=-O2,-march=native,-DNDEBUG
                -Xcompiler=-Wall,-Wextra,-Wshadow,-Wconversion
            )
        else()
            target_compile_options(kestrel-maxwell-tests PRIVATE
                -O2
                -march=native
                -DNDEBUG
                ${MAXWELL_WARN_FLAGS}
            )
        endif()
    endif()

    if(ipo_supported AND IS_TEST AND NOT KESTREL_CUDA)
        set_property(
            TARGET kestrel-maxwell-tests
            PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE
        )
    endif()

    add_test(NAME maxwell-aligned-soa      COMMAND kestrel-maxwell-tests AlignedSoA)
    add_test(NAME maxwell-grid             COMMAND kestrel-maxwell-tests Grid)
    add_test(NAME maxwell-source           COMMAND kestrel-maxwell-tests Source)
    add_test(NAME maxwell-pml              COMMAND kestrel-maxwell-tests PML)
    add_test(NAME maxwell-integration      COMMAND kestrel-maxwell-tests Integration)
    add_test(NAME maxwell-output           COMMAND kestrel-maxwell-tests Output)
    add_test(NAME maxwell-validation-class COMMAND kestrel-maxwell-tests ValidationClass)
    add_test(NAME maxwell-all-tests        COMMAND kestrel-maxwell-tests)

    if(IS_DEBUG)
        set_tests_properties(maxwell-aligned-soa      PROPERTIES TIMEOUT 60)
        set_tests_properties(maxwell-grid             PROPERTIES TIMEOUT 60)
        set_tests_properties(maxwell-source           PROPERTIES TIMEOUT 60)
        set_tests_properties(maxwell-output           PROPERTIES TIMEOUT 60)
        set_tests_properties(maxwell-pml              PROPERTIES TIMEOUT 300)
        set_tests_properties(maxwell-integration      PROPERTIES TIMEOUT 600)
        set_tests_properties(maxwell-validation-class PROPERTIES TIMEOUT 300)
        set_tests_properties(maxwell-all-tests        PROPERTIES TIMEOUT 600)
    else()
        set_tests_properties(maxwell-aligned-soa      PROPERTIES TIMEOUT 30)
        set_tests_properties(maxwell-grid             PROPERTIES TIMEOUT 30)
        set_tests_properties(maxwell-source           PROPERTIES TIMEOUT 30)
        set_tests_properties(maxwell-output           PROPERTIES TIMEOUT 30)
        set_tests_properties(maxwell-pml              PROPERTIES TIMEOUT 120)
        set_tests_properties(maxwell-integration      PROPERTIES TIMEOUT 300)
        set_tests_properties(maxwell-validation-class PROPERTIES TIMEOUT 120)
        set_tests_properties(maxwell-all-tests        PROPERTIES TIMEOUT 300)
    endif()
endif()

# Maxwell benchmarks:

set(MAXWELL_BENCHMARK_SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/benchmark/benchmark.cu")
kestrel_configure_sources("${MAXWELL_BENCHMARK_SOURCE}")

add_executable(kestrel-maxwell-benchmark "${MAXWELL_BENCHMARK_SOURCE}")
target_link_libraries(kestrel-maxwell-benchmark PRIVATE kestrel-maxwell-core)

if(ipo_supported AND NOT IS_DEBUG AND NOT KESTREL_CUDA)
    set_property(
        TARGET kestrel-maxwell-benchmark
        PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE
    )
endif()

add_library(kestrel-maxwell-core-serial STATIC ${MAXWELL_CORE_SOURCES})
target_include_directories(kestrel-maxwell-core-serial PUBLIC "${MAXWELL_SOURCE_DIR}")
target_link_libraries(kestrel-maxwell-core-serial PUBLIC xpu::xpu)

if(IS_DEBUG)
    if(KESTREL_CUDA)
        target_compile_options(kestrel-maxwell-core-serial PRIVATE
            -O0 -g -Xcompiler=-Wall,-Wextra
        )
    else()
        target_compile_options(kestrel-maxwell-core-serial PRIVATE
            -O0 -g -Wall -Wextra -Wpedantic
        )
    endif()
else()
    if(KESTREL_CUDA)
        target_compile_options(kestrel-maxwell-core-serial PRIVATE
            -O3 -Xcompiler=-O3,-march=native,-DNDEBUG
        )
    else()
        target_compile_options(kestrel-maxwell-core-serial PRIVATE
            -O3 -march=native -DNDEBUG -flto
        )
        target_link_options(kestrel-maxwell-core-serial PUBLIC -flto)
    endif()
endif()

if(ipo_supported AND NOT IS_DEBUG AND NOT KESTREL_CUDA)
    set_property(
        TARGET kestrel-maxwell-core-serial
        PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE
    )
endif()

add_executable(kestrel-maxwell-benchmark-serial "${MAXWELL_BENCHMARK_SOURCE}")
target_link_libraries(
    kestrel-maxwell-benchmark-serial
    PRIVATE kestrel-maxwell-core-serial
)

if(ipo_supported AND NOT IS_DEBUG AND NOT KESTREL_CUDA)
    set_property(
        TARGET kestrel-maxwell-benchmark-serial
        PROPERTY INTERPROCEDURAL_OPTIMIZATION TRUE
    )
endif()
