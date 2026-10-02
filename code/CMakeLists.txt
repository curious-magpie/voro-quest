cmake_minimum_required(VERSION 3.20)
project(rvd VERSION 0.1.0 LANGUAGES C CXX)

# RVD: Voronoi and restricted Voronoi diagrams you can look inside. The design
# is ../design_document/DESIGN.tex; this is milestone M0 of it -- the build,
# the host and its windows, a box and a cloud of seeds.

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR})

option(RVD_BUILD_TESTS "Build the standalone test binaries" ON)
option(RVD_BUILD_VIEWER "Build the viewer host (needs mesh-viewer)" ON)

# The flags every target of ours gets, in one place.
#
# -ffp-contract=off: no fusing of a*b + c into one FMA instruction. GCC does it
# by default wherever the hardware has FMA (every ARM machine, and x86 built
# for a recent CPU), and it changes the last bit of the result. Two things
# depend on that bit not moving: the seeds a random seed draws, which must be
# the same on every machine, and, from M1, the exact predicates, whose
# error-free products and sums are only error-free without fusing.
add_library(rvd_flags INTERFACE)
target_compile_options(rvd_flags INTERFACE -Wall -Wextra -ffp-contract=off)

# glm is header-only; prefer the config package, fall back to the bare header path.
find_package(glm QUIET)

# -------------------------------------------------------------------- mesh-dev
#
# The pieces of mesh-dev RVD reuses (design, section 4.2): the mesh readers,
# SurfaceMesh and TetMesh with their adjacency, the Morton order, and the test
# harness. Taken by path and compiled here, never copied, because copies drift.
# They need nothing but glm and two of mesh-dev's own small headers.
#
# RVD_MESH_DEV_DIR points at the checkout. The default is the one next to this
# project:
#
#   cmake -S . -B build -DRVD_MESH_DEV_DIR=$HOME/phd/my_software/mesh-dev
get_filename_component(RVD_DEFAULT_MESH_DEV_DIR
                       "${PROJECT_SOURCE_DIR}/../../mesh-dev" ABSOLUTE)
set(RVD_MESH_DEV_DIR "${RVD_DEFAULT_MESH_DEV_DIR}" CACHE PATH
    "The mesh-dev checkout whose io/ and mesh/ RVD builds against")
if(NOT EXISTS "${RVD_MESH_DEV_DIR}/src/mesh/tet_mesh.h")
    message(FATAL_ERROR "no mesh-dev in ${RVD_MESH_DEV_DIR}: point "
                        "RVD_MESH_DEV_DIR at a checkout of it")
endif()

add_library(rvd_mesh_dev STATIC
    ${RVD_MESH_DEV_DIR}/src/io/mesh_reader.cpp
    ${RVD_MESH_DEV_DIR}/src/io/msh_reader.cpp
    ${RVD_MESH_DEV_DIR}/src/io/obj_reader.cpp
    ${RVD_MESH_DEV_DIR}/src/mesh/morton_order.cpp
    ${RVD_MESH_DEV_DIR}/src/mesh/surface_mesh.cpp
    ${RVD_MESH_DEV_DIR}/src/mesh/tet_mesh.cpp
)
target_include_directories(rvd_mesh_dev PUBLIC ${RVD_MESH_DEV_DIR}/src)
target_link_libraries(rvd_mesh_dev PRIVATE rvd_flags)
if(glm_FOUND)
    target_link_libraries(rvd_mesh_dev PUBLIC glm::glm)
endif()

# --------------------------------------------------------------------- the core
#
# Everything that does not touch a window or the GPU. The app links it, and so
# does every test; a machine with no display builds and tests all of it.
add_library(rvd_core STATIC
    src/domain/box.cpp
    src/points/samplers.cpp
    src/session/rvd_session.cpp
)
target_include_directories(rvd_core PUBLIC src)
target_link_libraries(rvd_core PUBLIC rvd_mesh_dev)
target_link_libraries(rvd_core PRIVATE rvd_flags)

# ------------------------------------------------------------------ the viewer
#
# mesh-viewer, as mesh-dev takes it: the window, the drawing, and the Inspect
# and View windows. It carries GLFW, GL, glad and Dear ImGui, none of which the
# core can reach, so RVD_BUILD_VIEWER=OFF leaves a build that needs none of
# them.
#
# RVD_VIEWER_DIR picks the checkout: ext/mesh-viewer once that submodule
# exists, the sibling checkout next to this project until then.
if(RVD_BUILD_VIEWER)
    if(EXISTS "${PROJECT_SOURCE_DIR}/ext/mesh-viewer/CMakeLists.txt")
        set(RVD_DEFAULT_VIEWER_DIR "${PROJECT_SOURCE_DIR}/ext/mesh-viewer")
    else()
        get_filename_component(RVD_DEFAULT_VIEWER_DIR
                               "${PROJECT_SOURCE_DIR}/../../mesh-viewer" ABSOLUTE)
    endif()
    set(RVD_VIEWER_DIR "${RVD_DEFAULT_VIEWER_DIR}" CACHE PATH
        "The mesh-viewer checkout to build the host against")
    if(NOT EXISTS "${RVD_VIEWER_DIR}/CMakeLists.txt")
        message(FATAL_ERROR "no mesh-viewer in ${RVD_VIEWER_DIR}: point "
                            "RVD_VIEWER_DIR at a checkout, or configure with "
                            "-DRVD_BUILD_VIEWER=OFF for the core and the tests")
    endif()
    add_subdirectory("${RVD_VIEWER_DIR}" mesh-viewer)

    add_executable(rvd
        src/main.cpp

        # app/    the host: what the viewer is shown, and RVD's windows
        src/app/display.cpp
        src/app/rvd_source.cpp
        src/app/windows.cpp
    )
    target_link_libraries(rvd PRIVATE rvd_core rvd_flags mesh_viewer::mesh_viewer)
endif()

# ------------------------------------------------------------------ the tests
#
# Standalone binaries in mesh-dev's style: no window, no data files, no test
# framework, one source file per suite, and mesh-dev's own check() harness
# (tests/test_support.h there). A suite is added by dropping a file into tests/
# and naming it here.
if(RVD_BUILD_TESTS)
    enable_testing()

    # GMP, for the oracle: every exact computation in kernel/ is checked
    # against rational arithmetic that is exact by construction (tests/exact.h).
    # Only the tests link it; the kernel itself never does.
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(GMPXX REQUIRED IMPORTED_TARGET gmpxx)

    set(RVD_TEST_SUITES
        domain
        points
    )

    # The kernel suites join once the header they test exists, so the build
    # stays green while that header is being written. CONFIGURE_DEPENDS makes
    # the build notice a new header and re-run this without being asked.
    file(GLOB RVD_KERNEL_HEADERS CONFIGURE_DEPENDS src/kernel/*.h)
    # The kernel suites join once what they test exists, so the build stays
    # green while it is being written. A suite tests either a whole header,
    # joining once the file exists, or a function added to one, joining once
    # the function's name appears in it. The kernel headers are configure
    # dependencies, so the build re-runs these checks by itself whenever one
    # of them changes.
    #
    #   suite      header        function
    set(RVD_KERNEL_SUITES
        "eft       eft           -"
        "expansion expansion     -"
        "fast_sum  expansion     fast_sum"
        "linear_sum expansion    linear_sum"
        "product   expansion     product"
        "orient2d  orient2d      -"
        "bounded   bounded       -"
        "orient2d_filter orient2d orient2d_counts"
        "orient3d  orient3d      -"
        "dynamic_expansion dynamic_expansion -"
        "dynamic_linear_sum dynamic_expansion linear_sum"
        "side1     side          side1"
        "side2     side          side2"
        "side3     side          side3"
        "side4     side          side4"
    )
    foreach(entry IN LISTS RVD_KERNEL_SUITES)
        string(REGEX REPLACE " +" ";" entry "${entry}")
        list(GET entry 0 suite)
        list(GET entry 1 header)
        list(GET entry 2 function)
        set(path "${PROJECT_SOURCE_DIR}/src/kernel/${header}.h")
        if(NOT EXISTS "${path}")
            message(STATUS "test suite ${suite} skipped: "
                           "src/kernel/${header}.h not written yet")
            continue()
        endif()
        set_property(DIRECTORY APPEND PROPERTY CMAKE_CONFIGURE_DEPENDS "${path}")
        if(NOT function STREQUAL "-")
            file(READ "${path}" text)
            # A whole name, so that two_product( does not count as product(.
            if(NOT text MATCHES "(^|[^A-Za-z0-9_])${function}[ \t\r\n]*\\(")
                message(STATUS "test suite ${suite} skipped: "
                               "${function}() not in src/kernel/${header}.h yet")
                continue()
            endif()
        endif()
        list(APPEND RVD_TEST_SUITES ${suite})
    endforeach()

    foreach(suite IN LISTS RVD_TEST_SUITES)
        add_executable(rvd-tests-${suite} tests/test_${suite}.cpp)
        target_include_directories(rvd-tests-${suite} PRIVATE
                                   ${RVD_MESH_DEV_DIR}/tests)
        target_link_libraries(rvd-tests-${suite} PRIVATE rvd_core rvd_flags
                                                         PkgConfig::GMPXX)
        # Asserts stay on in the tests, Release or not. The kernel is
        # header-only, so its preconditions (fast_two_sum's |a| >= |b|) are
        # compiled into each test and checked on every input the test makes:
        # a precondition broken by a caller often still gives an exact-looking
        # answer, and only the assert says so.
        target_compile_options(rvd-tests-${suite} PRIVATE -UNDEBUG)
        add_test(NAME ${suite} COMMAND rvd-tests-${suite})
    endforeach()
endif()
