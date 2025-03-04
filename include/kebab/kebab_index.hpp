#ifndef KEBAB_INDEX_HPP
#define KEBAB_INDEX_HPP

#include "kebab/nt_hash.hpp"
#include "kebab/bloom_filter.hpp"

#include "external/kseq.h"

#include "constants.hpp"

#include <string>
#include <vector>
#include <fstream>

namespace kebab {

struct Fragment {
    size_t start;
    size_t length;

    bool operator<(const Fragment& other) const {
        return length > other.length; // > for descending order
    }
};

class KebabIndex {
public:
    KebabIndex(size_t k, size_t expected_kmers, double fp_rate, size_t num_hashes = DEFAULT_HASH_FUNCS);
    explicit KebabIndex(std::istream& in);

    size_t get_k() const { return k; }

    void add_sequence(const char* seq, size_t len);
    std::vector<Fragment> scan_read(const char* seq, size_t len, uint64_t min_mem_length, bool remove_overlaps = DEFAULT_REMOVE_OVERLAPS);

    std::string get_stats() const;

    std::string get_file_extension() const { return FILE_EXTENSION; }
    void save(std::ostream& out) const;
    void load(std::istream& in);

private:
    static constexpr const char* FILE_EXTENSION = ".kbb";

    size_t k;
    NtHash<> hasher;
    BloomFilter<> bf;
};

} // namespace kebab

#endif // KEBAB_INDEX_HPP
