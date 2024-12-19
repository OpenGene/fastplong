# CMakeLists.txt
cmake_minimum_required(VERSION 3.21)
project(fastplong VERSION 0.2.1)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Find dependencies provided by vcpkg (via vcpkg.cmake)
find_package(hwy CONFIG REQUIRED)
find_package(benchmark CONFIG REQUIRED)
find_package(libdeflate CONFIG REQUIRED)

if (MSVC)
    find_package(unofficial-isal CONFIG REQUIRED)
    set(
        FASTPLONG_LIBS
        hwy::hwy
        unofficial::isal::isal
        $<IF:$<TARGET_EXISTS:libdeflate::libdeflate_shared>,libdeflate::libdeflate_shared,libdeflate::libdeflate_static>
    )
else()
    set(
        FASTPLONG_LIBS
        hwy::hwy
        $<IF:$<TARGET_EXISTS:libdeflate::libdeflate_shared>,libdeflate::libdeflate_shared,libdeflate::libdeflate_static>
        -lisal
        -lpthread
    )
endif()

if (CMAKE_CXX_COMPILER_ARCHITECTURE_ID MATCHES "x86_64")
    add_compile_options(-march=haswell -maes)
    # add_link_options(-pg)
    ## HWY_AVX3 which is 512
    # SET(GCC_COVERAGE_COMPILE_FLAGS "-march=sapphirerapids")
endif()
if (MSVC)
    add_compile_options(/arch:AVX2)
endif()

set(
    FASTPLONG_SOURCES
    src/adaptertrimmer.cpp
    src/editdistance.cpp
    src/evaluator.cpp
    src/fastareader.cpp
    src/fastqreader.cpp
    src/filter.cpp
    src/filterresult.cpp
    src/jsonreporter.cpp
    src/htmlreporter.cpp
    src/nucleotidetree.cpp
    src/options.cpp
    src/polyx.cpp
    src/processor.cpp
    src/read.cpp
    src/readpool.cpp
    src/seprocessor.cpp
    src/sequence.cpp
    src/stats.cpp
    src/threadconfig.cpp
    src/writer.cpp
    src/writerthread.cpp
)

set(
    FASTPLONG_TEST_SOURCES
    test/adaptertrimmer_test.cpp
    test/fastareader_test.cpp
    test/filter_test.cpp
    test/nucleotidetree_test.cpp
    test/read_test.cpp
    test/evaluator_test.cpp
    test/fastqreader_test.cpp
    test/globals.cpp
    test/polyx_test.cpp
    test/sequence_test.cpp
    test/stats_test.cpp
)
include_directories(
    # intel mac install location
    /usr/local/include
    # arm64 mac install location
    /opt/homebrew/include
)
link_directories(
    # intel mac install location
    /usr/local/lib
    # arm64 mac install location
    /opt/homebrew/lib
)

add_library(fastplong_core ${FASTPLONG_SOURCES})
target_link_libraries(
    fastplong_core PRIVATE
    ${FASTPLONG_LIBS}
)

if (NOT MSVC)
    # benchmark target
    add_executable(fastplong_benchmarks
        benchmarks/sequence_benchmark.cpp
        test/globals.cpp
    )
    target_link_libraries(
        fastplong_benchmarks PRIVATE
        ${FASTPLONG_LIBS}
        fastplong_core
        benchmark::benchmark
        benchmark::benchmark_main
    )
endif()
# main target
add_executable(fastplong src/main.cpp)
target_link_libraries(
    fastplong PUBLIC
    fastplong_core
    ${FASTPLONG_LIBS}
)
# tests target
enable_testing()
find_package(GTest CONFIG REQUIRED)
configure_file(test/config.h.in ${CMAKE_BINARY_DIR}/include/config.h)
add_executable(fastplong_tests ${FASTPLONG_TEST_SOURCES})
target_include_directories(
    fastplong_tests
    PRIVATE ${CMAKE_BINARY_DIR}/include
)
target_link_libraries(
    fastplong_tests PRIVATE
    ${FASTPLONG_LIBS}
    fastplong_core
    GTest::gtest
    GTest::gtest_main
)
add_test(fastplong_tests fastplong_tests)

install(TARGETS fastplong
    CONFIGURATIONS Release
    RUNTIME
    DESTINATION bin
)
