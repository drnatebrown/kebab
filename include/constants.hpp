#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

// VERSION
static constexpr const char* VERSION = "0.2.0";

// I/O
static constexpr size_t DEFAULT_BUFFER_SIZE = 64ULL * 1024ULL * 1024ULL; // 64MB

// ESTIMATE
static constexpr uint64_t HLL_SIZE = 20; // 2^20 bytes

// BUILD
static constexpr uint16_t DEFAULT_KMER_SIZE = 20;
static constexpr double DEFAULT_FP_RATE = 0.1;
static constexpr uint16_t DEFAULT_HASH_FUNCS = 0;
static constexpr uint64_t DEFAULT_EXPECTED_KMERS = 0;

// SCAN
static constexpr uint64_t DEFAULT_MIN_MEM_LENGTH = 20;
static constexpr bool DEFAULT_SORT_FRAGMENTS = false;
static constexpr bool DEFAULT_REMOVE_OVERLAPS = false;

#endif