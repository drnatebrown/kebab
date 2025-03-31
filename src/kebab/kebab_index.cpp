#include "kebab/kebab_index.hpp"

namespace kebab {

KebabIndex::KebabIndex(size_t k, size_t expected_kmers, double fp_rate, size_t num_hashes, bool canonical)
    : k(k)
    , canonical(canonical)
    , hasher(k, canonical)
    , bf(expected_kmers, fp_rate, num_hashes)
{
}

KebabIndex::KebabIndex(std::istream& in) 
    : k(0)
    , canonical(DEFAULT_CANONICAL)
    , hasher(0, DEFAULT_CANONICAL)
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

std::vector<Fragment> KebabIndex::scan_read(const char* seq, size_t len, uint64_t min_mem_length, bool remove_overlaps) {
    if (min_mem_length < k) {
        throw std::invalid_argument("min_mem_length (" + std::to_string(min_mem_length) + ") must be greater than or equal to k (" + std::to_string(k) + ")");
    }

    hasher.set_sequence(seq, len);
    std::vector<Fragment> fragments;
    // fragments.reserve(remove_overlaps ? std::ceil(static_cast<double>(len) / min_mem_length) : len - k + 1);

    size_t start = 0;
    size_t last_frag_end = 0;

    // end is exclusive
    auto update_fragments = [&](size_t frag_end) {
        if (frag_end - start >= min_mem_length) {
            // Check if overlaps the last fragment
            if (remove_overlaps && start < last_frag_end) {
                fragments.back().length += frag_end - last_frag_end;
            } else {
                fragments.push_back({start, frag_end - start});
            }
            last_frag_end = frag_end;
        }
    };

    // k-mer identified by position of last character
    for (size_t i = k - 1; i < len; ++i) {
        if (!bf.contains(hasher.hash())) {
            update_fragments(i);
            start = i - k + 2; // i - (k - 1) + 1 -> move to start of k-mer, plus one to move past the offending k-mer
        }
        hasher.unsafe_roll();
    }
    update_fragments(len);

    return fragments;
}

std::string KebabIndex::get_stats() const {
    return  "\tk: " + std::to_string(k) + "\n" 
            + bf.get_stats();
}

void KebabIndex::save(std::ostream& out) const {
    out.write(reinterpret_cast<const char*>(&k), sizeof(k));
    out.write(reinterpret_cast<const char*>(&canonical), sizeof(canonical));
    bf.save(out);
}

void KebabIndex::load(std::istream& in) {
    in.read(reinterpret_cast<char*>(&k), sizeof(k));
    in.read(reinterpret_cast<char*>(&canonical), sizeof(canonical));
    hasher = NtHash<>(k, canonical);
    bf.load(in);
}

} // namespace kebab
