# CMakeLists.txt
cmake_minimum_required(VERSION 3.21)
project(fastplong VERSION 0.2.1)
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
)
endif()

if (CMAKE_CXX_COMPILER_ARCHITECTURE_ID MATCHES "x86_64")
    add_compile_options(-march=haswell -maes)
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
)
link_directories(
    /usr/local/lib
    /opt/homebrew/lib
)

if (NOT MSVC)
    # benchmark target
    add_executable(fastplong_benchmarks)
    target_sources(
        fastplong_benchmarks PRIVATE
        benchmarks/sequence_benchmark.cpp
        test/globals.cpp
        ${FASTPLONG_SOURCES}
    )
    target_link_libraries(
        fastplong_benchmarks PRIVATE
        ${FASTPLONG_LIBS}
        benchmark::benchmark
        benchmark::benchmark_main
    )
endif()
# main target
add_executable(fastplong src/main.cpp ${FASTPLONG_SOURCES})
target_link_libraries(
    fastplong PUBLIC
    ${FASTPLONG_LIBS}
)
# tests target
enable_testing()
find_package(GTest CONFIG REQUIRED)
add_executable(fastplong_tests ${FASTPLONG_SOURCES} ${FASTPLONG_TEST_SOURCES})
set_property(TARGET fastplong_tests PROPERTY CXX_STANDARD 23)
target_link_libraries(
    fastplong_tests PRIVATE
    ${FASTPLONG_LIBS}
    GTest::gtest
    GTest::gtest_main
    -lpthread
)
add_test(fastplong_tests fastplong_tests)