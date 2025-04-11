#include "kebab/kebab_index.hpp"

namespace kebab {

template<typename Filter>
KebabIndex<Filter>::KebabIndex(size_t k, size_t expected_kmers, double fp_rate, size_t num_hashes, KmerMode kmer_mode, FilterSizeMode filter_size_mode)
    : k(k)
    , kmer_mode(kmer_mode)
    , build_rev_comp(use_build_rev_comp(kmer_mode))
    , scan_rev_comp(use_scan_rev_comp(kmer_mode))
    , build_hasher(k, build_rev_comp)
    , scan_hasher(k, scan_rev_comp)
    , bf(expected_kmers, fp_rate, num_hashes, filter_size_mode)
{
}

template<typename Filter>
KebabIndex<Filter>::KebabIndex(std::istream& in) 
    : k(0)
    , kmer_mode(DEFAULT_KMER_MODE)
    , build_rev_comp(use_build_rev_comp(kmer_mode))
    , scan_rev_comp(use_scan_rev_comp(kmer_mode))
    , build_hasher(0, build_rev_comp)
    , scan_hasher(0, scan_rev_comp)
    , bf()
{
    load(in);
}

template<typename Filter>
void KebabIndex<Filter>::add_sequence(const char* seq, size_t len) {
    build_hasher.set_sequence(seq, len);
    for (size_t i = 0; i < len - k + 1; ++i) {
        switch (kmer_mode) {
            case KmerMode::FORWARD_ONLY:
                bf.add(build_hasher.hash());
                break;
            case KmerMode::BOTH_STRANDS:
                bf.add(build_hasher.hash());
                bf.add(build_hasher.hash_rc());
                break;
            case KmerMode::CANONICAL_ONLY:
                bf.add(build_hasher.hash_canonical());
                break;
        }
        build_hasher.unsafe_roll();
    }
    // do {
    //     bf.add(hasher.hash());
    // } while (hasher.roll());
}

template<typename Filter>
std::vector<Fragment> KebabIndex<Filter>::scan_read(const char* seq, size_t len, uint64_t min_mem_length, bool remove_overlaps) {
    if (min_mem_length < k) {
        throw std::invalid_argument("min_mem_length (" + std::to_string(min_mem_length) + ") must be greater than or equal to k (" + std::to_string(k) + ")");
    }

    scan_hasher.set_sequence(seq, len);
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
        if (!bf.contains(scan_rev_comp ? scan_hasher.hash_canonical() : scan_hasher.hash())) {
            update_fragments(i);
            start = i - k + 2; // i - (k - 1) + 1 -> move to start of k-mer, plus one to move past the offending k-mer
        }
        scan_hasher.unsafe_roll();
    }
    update_fragments(len);

    return fragments;
}

template<typename Filter>
std::string KebabIndex<Filter>::get_stats() const {
    return  "\tk: " + std::to_string(k) + "\n" 
            + bf.get_stats();
}

template<typename Filter>
void KebabIndex<Filter>::save(std::ostream& out) const {
    out.write(reinterpret_cast<const char*>(&k), sizeof(k));
    out.write(reinterpret_cast<const char*>(&kmer_mode), sizeof(kmer_mode));
    bf.save(out);
}

template<typename Filter>
void KebabIndex<Filter>::load(std::istream& in) {
    in.read(reinterpret_cast<char*>(&k), sizeof(k));
    in.read(reinterpret_cast<char*>(&kmer_mode), sizeof(kmer_mode));
    build_rev_comp = use_build_rev_comp(kmer_mode);
    scan_rev_comp = use_scan_rev_comp(kmer_mode);
    build_hasher = NtHash<>(k, build_rev_comp);
    scan_hasher = NtHash<>(k, scan_rev_comp);
    bf.load(in);
}

// Explicit instantiation
template class KebabIndex<ShiftFilter>;
template class KebabIndex<ModFilter>;

} // namespace kebab
