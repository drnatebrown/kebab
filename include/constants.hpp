#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP
// VERSION
static constexpr const char* VERSION = "1.0.1";

// I/O
static constexpr size_t DEFAULT_BUFFER_SIZE = 64ULL * 1024ULL * 1024ULL; // 64MB
static constexpr const char* KEBAB_FILE_SUFFIX = ".kbb";

// K-mer Mode
enum class KmerMode {
    BOTH_STRANDS,          // Include both forward and reverse complement
    CANONICAL_ONLY,        // Use canonical form of each k-mer
    FORWARD_ONLY           // Only use forward k-mers
};
inline bool use_build_rev_comp(KmerMode mode) { return mode == KmerMode::BOTH_STRANDS || mode == KmerMode::CANONICAL_ONLY; }
inline bool use_scan_rev_comp(KmerMode mode) { return mode == KmerMode::CANONICAL_ONLY; }
static constexpr bool DEFAULT_REVERSE_COMPLEMENT = true;

// Filter Size Mode
enum class FilterSizeMode {
    NEXT_POWER_OF_TWO,
    PREVIOUS_POWER_OF_TWO,
    EXACT
};
inline bool use_shift_filter(FilterSizeMode mode) { return mode == FilterSizeMode::NEXT_POWER_OF_TWO || mode == FilterSizeMode::PREVIOUS_POWER_OF_TWO; }
static constexpr double ROUND_THRESHOLD = 0.10; // 10% tolerance to round to the nearest power of two despite mode

// ESTIMATE
static constexpr uint64_t HLL_SIZE = 20; // 2^20 bytes

// LATENCY HIDING
static constexpr uint64_t PREFETCH_DISTANCE = 32; // prefetch this many read operations on the bloom filter

// BUILD
static constexpr uint16_t DEFAULT_KMER_SIZE = 20;
static constexpr KmerMode DEFAULT_KMER_MODE = KmerMode::CANONICAL_ONLY;
static constexpr double DEFAULT_FP_RATE = 0.1;
static constexpr uint16_t DEFAULT_HASH_FUNCS = 0; // 0 means optimal number of hash functions
static constexpr uint64_t DEFAULT_EXPECTED_KMERS = 0; // 0 means use hyperloglog to estimate number of k-mers
static constexpr uint16_t DEFAULT_BUILD_THREADS = 8; // overridden by call to omp_get_max_threads()
static constexpr FilterSizeMode DEFAULT_FILTER_SIZE_MODE = FilterSizeMode::PREVIOUS_POWER_OF_TWO;

// SCAN
static constexpr uint64_t DEFAULT_MIN_MEM_LENGTH = 25;
static constexpr uint16_t DEFAULT_TOP_T = 0; // 0 means no top-t filtering
static constexpr bool DEFAULT_SORT_FRAGMENTS = false;
static constexpr bool DEFAULT_REMOVE_OVERLAPS = false;
static constexpr bool DEFAULT_PREFETCH = true;
static constexpr uint16_t DEFAULT_SCAN_THREADS = 8; // overridden by call to omp_get_max_threads()

#endif
