#include "kebab/kebab_index.hpp"

namespace kebab {

KebabIndex::KebabIndex(size_t k, size_t expected_kmers,double fp_rate, size_t num_hashes)
    : k(k)
    , hasher(k)
    , bf(expected_kmers, fp_rate, num_hashes)
{
}

KebabIndex::KebabIndex(std::istream& in) 
    : k(0)
    , hasher()
    , bf()
{
    load(in);
}

void KebabIndex::add_sequence(const char* seq, size_t len) {
    hasher.set_sequence(seq, len);
    for (size_t i = 0; i < len - k + 1; ++i) {
        bf.add(hasher.hash());
        hasher.unsafe_roll();
    }
    // do {
    //     bf.add(hasher.hash());
    // } while (hasher.roll());
}

std::vector<Fragment> KebabIndex::scan_read(const char* seq, size_t len, uint64_t min_mem_length) {
    if (min_mem_length < k) {
        throw std::invalid_argument("min_mem_length (" + std::to_string(min_mem_length) + ") must be greater than or equal to k (" + std::to_string(k) + ")");
    }

    hasher.set_sequence(seq, len);
    std::vector<Fragment> fragments;

    size_t leftEnd = 0;
    
    for (size_t i = k - 1; i < len; ++i) {
        if (!bf.contains(hasher.hash())) {
            if (i - leftEnd >= min_mem_length) {
                fragments.push_back({leftEnd, i - leftEnd});
            }
            leftEnd = i - k + 2; // +2 to move past the offending k-mer
        }
        hasher.unsafe_roll();
    }

    if ((len - leftEnd) >= min_mem_length) {
        fragments.push_back({leftEnd, len - leftEnd});
    }

    return fragments;
}

std::string KebabIndex::get_stats() const {
    return "KebabIndex Stats:\n"
           "\tk: " + std::to_string(k) + "\n" 
           + bf.get_stats();
}

void KebabIndex::save(std::ostream& out) const {
    out.write(reinterpret_cast<const char*>(&k), sizeof(k));
    bf.save(out);
}

void KebabIndex::load(std::istream& in) {
    in.read(reinterpret_cast<char*>(&k), sizeof(k));
    hasher = NtHash<>(k);
    bf.load(in);
}

} // namespace kebab
