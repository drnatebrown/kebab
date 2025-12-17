#include "kebab/kebab_index.hpp"

namespace kebab {

template<typename Filter>
KebabIndex<Filter>::KebabIndex(size_t k, size_t expected_kmers, double fp_rate, size_t num_hashes, KmerMode kmer_mode, FilterSizeMode filter_size_mode)
    : k(k)
    , kmer_mode(kmer_mode)
    , build_rev_comp(use_build_rev_comp(kmer_mode))
    , scan_rev_comp(use_scan_rev_comp(kmer_mode))
    , bf(expected_kmers, fp_rate, num_hashes, filter_size_mode)
{
}

template<typename Filter>
KebabIndex<Filter>::KebabIndex(std::istream& in) 
    : k(0)
    , kmer_mode(DEFAULT_KMER_MODE)
    , build_rev_comp(use_build_rev_comp(kmer_mode))
    , scan_rev_comp(use_scan_rev_comp(kmer_mode))
    , bf()
{
    load(in);
}

template<typename Filter>
void KebabIndex<Filter>::add_sequence(const char* seq, size_t len) {
    thread_local static NtHash<> build_hasher(k, build_rev_comp);

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
std::vector<Fragment> KebabIndex<Filter>::scan_read(const char* seq, size_t len, uint64_t min_mem_length, bool remove_overlaps, bool prefetch) {
    thread_local static NtHash<> scan_hasher(k, scan_rev_comp);

    if (prefetch) {
        return scan_read_prefetch(seq, len, scan_hasher, min_mem_length, remove_overlaps);
    }
    else {
        return scan_read(seq, len, scan_hasher, min_mem_length, remove_overlaps);
    }
}

template<typename Filter>
std::vector<Fragment> KebabIndex<Filter>::scan_read(const char* seq, size_t len, NtHash<>& scan_hasher, uint64_t min_mem_length, bool remove_overlaps) {
    if (min_mem_length <= k) {
        throw std::invalid_argument("min_mem_length (" + std::to_string(min_mem_length) + ") must be greater than k (" + std::to_string(k) + ")");
    }

    scan_hasher.set_sequence(seq, len);
    std::vector<Fragment> fragments;
    // fragments.reserve(remove_overlaps ? std::ceil(static_cast<double>(len) / min_mem_length) : len - k + 1);

    size_t start = 0;
    size_t last_frag_end = 0;
    size_t i = 0;
    size_t numb_of_kmers = len-k+1;
    size_t buffer_capacity = min_mem_length-k+1;

    std::vector<uint64_t> buffer;
    buffer.reserve(buffer_capacity);

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

    auto fill_buffer = [&]() {
      for (size_t bi = 0; bi < buffer_capacity && i < numb_of_kmers; bi++, i++ ){
        buffer[bi] = scan_hash(scan_hasher);
        scan_hasher.unsafe_roll();
      }
    };

    // back scanning of the buffer, decides on the start position
    auto buffer_scan = [&]() {
      for (size_t bi = buffer_capacity-1; bi<buffer_capacity; bi--) {
        if (!bf.contains(buffer[bi])) {
          // (buffer_capacity-bi-1) gives the number of kmers that passed the back check
          start = i - (buffer_capacity-bi-1);
          break;
        }
      }
    };
 
    // decides on the end position and reports fragment once end is reached
    auto forward_scan = [&]() {
      for (; i < numb_of_kmers; i++) {
        if (!bf.contains(scan_hash(scan_hasher))) break;
        scan_hasher.unsafe_roll();
      }
      scan_hasher.unsafe_roll();
      update_fragments(i + (k-1));
      i++; // move past the offending kmer
    };

    while (i < numb_of_kmers) {
      start = i;
      fill_buffer();
      buffer_scan();
      forward_scan();
    }

    return fragments;
}

template<typename Filter>
std::vector<Fragment> KebabIndex<Filter>::scan_read_prefetch(const char* seq, size_t len, NtHash<>& scan_hasher, uint64_t min_mem_length, bool remove_overlaps) {
    if (min_mem_length <= k) {
        throw std::invalid_argument("min_mem_length (" + std::to_string(min_mem_length) + ") must be greater than k (" + std::to_string(k) + ")");
    }

    // Based on number of hashes to adequately spread out work done when prefetching
    const size_t NUM_PREFETCH_KMERS = PREFETCH_DISTANCE/bf.get_num_hashes();
    std::vector<PendingKmer> pending_kmers(NUM_PREFETCH_KMERS, PendingKmer(bf.get_num_hashes()));
    size_t pending_head = 0;
    size_t pending_tail = 0;
    size_t pending_count = 0;

    scan_hasher.set_sequence(seq, len);
    std::vector<Fragment> fragments;

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

    auto remove_pending_kmer = [&]() {
        if (!bf.check_prefetch(pending_kmers[pending_head].prefetch_info)) {
            update_fragments(pending_kmers[pending_head].pos);
            start = pending_kmers[pending_head].pos - k + 2;
        }
        pending_head = (pending_head + 1) % NUM_PREFETCH_KMERS;
        --pending_count;
    };

    auto add_pending_kmer = [&](size_t pos) {
        bf.prefetch_words(scan_hash(scan_hasher), pending_kmers[pending_tail].prefetch_info);
        pending_kmers[pending_tail].pos = pos;
        pending_tail = (pending_tail + 1) % NUM_PREFETCH_KMERS;
        ++pending_count;
    };

    // Prefetch initial k-mers
    add_pending_kmer(k - 1); // first doesn't need to be rolled
    for (size_t i = k; i < k - 1 + NUM_PREFETCH_KMERS && i < len; ++i) {
        scan_hasher.unsafe_roll();
        add_pending_kmer(i);
    }

    // Check fetched k-mer, prefetch next k-mer
    for (size_t i = k - 1 + NUM_PREFETCH_KMERS; i < len; ++i) {
        scan_hasher.unsafe_roll();
        remove_pending_kmer();
        add_pending_kmer(i);
    }

    // Check remaining pending k-mers
    while (pending_count > 0) {
        remove_pending_kmer();
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
    bf.load(in);
}

// Explicit instantiation
template class KebabIndex<ShiftFilter>;
template class KebabIndex<ModFilter>;

} // namespace kebab
