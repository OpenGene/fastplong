#include <benchmark/benchmark.h>
#include "../src/sequence.h"

static void BM_SequenceReverseSerial(benchmark::State& state) {
  // Perform setup here
  auto input = new std::string("AAAAAAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC");
  for (auto _ : state) {
    // This code gets timed
    const auto output = Sequence::reverseComplement(input);
  }
  delete input;
}

static void BM_SequenceReverseHwy(benchmark::State& state) {
  // Perform setup here
  auto input = new std::string("AAAAAAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC");
  for (auto _ : state) {
    // This code gets timed
    const auto output = Sequence::reverseComplementHwy(input);
  }
  delete input;
}

// Register the function as a benchmark
BENCHMARK(BM_SequenceReverseSerial);
BENCHMARK(BM_SequenceReverseHwy);
// Run the benchmark
// BENCHMARK_MAIN();