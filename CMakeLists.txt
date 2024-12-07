# CMakeLists.txt
cmake_minimum_required(VERSION 3.21)
project(fastplong VERSION 0.2.1)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Find dependencies provided by vcpkg (via vcpkg.cmake)
find_package(hwy CONFIG REQUIRED)
find_package(benchmark CONFIG REQUIRED)
find_package(libdeflate CONFIG REQUIRED)

if (CMAKE_CXX_COMPILER_ARCHITECTURE_ID MATCHES "x86_64")
    add_compile_options(-march=haswell -maes)
    ## HWY_AVX3 which is 512
    # SET(GCC_COVERAGE_COMPILE_FLAGS "-march=sapphirerapids")
endif()
if (MSVC)
    add_compile_options(/arch:AVX2 /Gv)
endif()

set(
    FASTPLONG_SOURCES
    src/adaptertrimmer.cpp
    src/adaptertrimmer.h
    src/cmdline.h
    src/common.h
    src/editdistance.cpp
    src/editdistance.h
    src/evaluator.cpp
    src/evaluator.h
    src/fastareader.cpp
    src/fastareader.h
    src/fastqreader.cpp
    src/fastqreader.h
    src/filter.cpp
    src/filter.h
    src/filterresult.cpp
    src/filterresult.h
    src/jsonreporter.cpp
    src/jsonreporter.h
    src/htmlreporter.cpp
    src/htmlreporter.h
    src/knownadapters.h
    src/nucleotidetree.cpp
    src/nucleotidetree.h
    src/options.h
    src/options.cpp
    src/polyx.h
    src/polyx.cpp
    src/processor.h
    src/processor.cpp
    src/read.cpp
    src/read.h
    src/readpool.cpp
    src/readpool.h
    src/seprocessor.cpp
    src/seprocessor.h
    src/sequence.cpp
    src/sequence.h
    src/simdutil.h
    src/singleproducersingleconsumerlist.h
    src/stats.cpp
    src/stats.h
    src/threadconfig.cpp
    src/threadconfig.h
    src/util.h
    src/writer.cpp
    src/writer.h
    src/writerthread.cpp
    src/writerthread.h
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

add_library(fastplong_lib STATIC ${FASTPLONG_SOURCES})
target_link_libraries(
    fastplong_lib PUBLIC
    hwy::hwy
    $<IF:$<TARGET_EXISTS:libdeflate::libdeflate_shared>,libdeflate::libdeflate_shared,libdeflate::libdeflate_static>
    -lpthread
    -lisal
)
set_property(TARGET fastplong_lib PROPERTY CXX_STANDARD 23)

# benchmark target
add_executable(fastplong_benchmarks)
target_sources(
    fastplong_benchmarks PRIVATE
    benchmarks/sequence_benchmark.cpp
    src/read.cpp
    test/globals.cpp
    src/sequence.cpp
)
target_link_libraries(
    fastplong_benchmarks PRIVATE
    hwy::hwy benchmark::benchmark benchmark::benchmark_main
    fastplong_lib
)
set_property(TARGET fastplong_benchmarks PROPERTY CXX_STANDARD 23)

# main target
add_executable(fastplong src/main.cpp)
target_link_libraries(
    fastplong PRIVATE
    fastplong_lib
)
set_property(TARGET fastplong PROPERTY CXX_STANDARD 23)
# tests target
enable_testing()
find_package(GTest CONFIG REQUIRED)
add_executable(fastplong_tests ${FASTPLONG_TEST_SOURCES})
set_property(TARGET fastplong_tests PROPERTY CXX_STANDARD 23)
target_link_libraries(
    fastplong_tests PRIVATE
    fastplong_lib
    GTest::gtest GTest::gtest_main
    -lpthread
    -lisal
)
add_test(fastplong_tests fastplong_tests)

