#include <benchmark/benchmark.h>
#include "../src/sequence.h"


// serial implementation of reverseComplement to compare
static string reverseComplementSerial(string* origin) {
    string str(origin->length(), 0);
    int len = origin->length();
    for (int c = 0; c < origin->length(); c++)
    {
      char base = (*origin)[c];
      switch (base)
      {
      case 'A':
      case 'a':
        str[len - c - 1] = 'T';
        break;
      case 'T':
      case 't':
        str[len - c - 1] = 'A';
        break;
      case 'C':
      case 'c':
        str[len - c - 1] = 'G';
        break;
      case 'G':
      case 'g':
        str[len - c - 1] = 'C';
        break;
      default:
        str[len - c - 1] = 'N';
      }
    }
    return str;
}

static void BM_SequenceReverseSerial(benchmark::State& state) {
  // Perform setup here
  auto input = new std::string("AAAAAAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC");
  for (auto _ : state) {
    const auto output = reverseComplementSerial(input);
    // This code gets timed
  }
  delete input;
}

static void BM_SequenceReverseSIMD(benchmark::State& state) {
  // Perform setup here
  auto input = new std::string("AAAAAAAAAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC");
  for (auto _ : state) {
    // This code gets timed
    const auto output = Sequence::reverseComplement(input);
  }
  delete input;
}

// Register the function as a benchmark
BENCHMARK(BM_SequenceReverseSerial);
BENCHMARK(BM_SequenceReverseSIMD);
// Run the benchmark
// BENCHMARK_MAIN();