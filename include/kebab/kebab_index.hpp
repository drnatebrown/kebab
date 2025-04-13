#ifndef KEBAB_INDEX_HPP
#define KEBAB_INDEX_HPP

#include "kebab/nt_hash.hpp"
#include "kebab/bloom_filter.hpp"

#include "external/kseq.h"

#include "constants.hpp"

#include <string>
#include <vector>
#include <fstream>
#include <cstdint>

namespace kebab {

struct Fragment {
    size_t start;
    size_t length;

    bool operator<(const Fragment& other) const {
        return length > other.length; // > for descending order
    }
};

template<typename Filter = ShiftFilter>
class KebabIndex {
public:
    KebabIndex(size_t k, size_t expected_kmers, double fp_rate, size_t num_hashes = DEFAULT_HASH_FUNCS, KmerMode kmer_mode = DEFAULT_KMER_MODE, FilterSizeMode filter_size_mode = DEFAULT_FILTER_SIZE_MODE);
    explicit KebabIndex(std::istream& in);

    size_t get_k() const { return k; }

    void add_sequence(const char* seq, size_t len);
    std::vector<Fragment> scan_read(const char* seq, size_t len, uint64_t min_mem_length, bool remove_overlaps = DEFAULT_REMOVE_OVERLAPS, bool prefetch = DEFAULT_PREFETCH);
    std::string get_stats() const;
    
    void save(std::ostream& out) const;
    void load(std::istream& in);

private:
    size_t k;
    KmerMode kmer_mode;
    bool build_rev_comp;
    bool scan_rev_comp;
    Filter bf;

    struct PendingKmer {
        PrefetchInfo prefetch_info;
        size_t pos;

        PendingKmer(size_t num_hashes) : prefetch_info(num_hashes), pos(0) {}
    };

    std::vector<Fragment> scan_read(const char* seq, size_t len, NtHash<>& scan_hasher, uint64_t min_mem_length, bool remove_overlaps = DEFAULT_REMOVE_OVERLAPS);
    std::vector<Fragment> scan_read_prefetch(const char* seq, size_t len, NtHash<>& scan_hasher, uint64_t min_mem_length, bool remove_overlaps = DEFAULT_REMOVE_OVERLAPS);

    uint64_t scan_hash(const NtHash<>& hasher) const {
        return (scan_rev_comp) ? hasher.hash_canonical() : hasher.hash();
    }
};

} // namespace kebab

#endif // KEBAB_INDEX_HPP
