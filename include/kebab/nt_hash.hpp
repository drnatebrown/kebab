#ifndef NTHASH_H
#define NTHASH_H

#include <cstdint>
#include <cstddef>
#include <climits>
#include <algorithm>
#include "constants.hpp"

namespace kebab {

template<typename T = uint64_t>
class NtHash {
public:
    explicit NtHash(size_t k, bool canonical = DEFAULT_CANONICAL) noexcept;
    NtHash() noexcept : k(0), canonical(DEFAULT_CANONICAL), seq(nullptr), len(0), pos(0), hash_val(0), hash_val_rc(0) {}

    void set_sequence(const char* seq, size_t len) noexcept;

    [[nodiscard]] size_t get_k() const noexcept { return k; }
    [[nodiscard]] size_t get_pos() const noexcept { return pos; }
    [[nodiscard]] size_t get_len() const noexcept { return len; }

    [[nodiscard]] static T get_max_hash() noexcept { return static_cast<T>(-1); } // maximum hash value for T

    [[nodiscard]] bool roll() noexcept;
    void unsafe_roll() noexcept;
    [[nodiscard]] T hash() const noexcept { return canonical ? std::min(hash_val, hash_val_rc) : hash_val; }

private:
    size_t k;
    bool canonical;

    const char* seq;
    size_t len;
    size_t pos;

    T hash_val;
    T hash_val_rc;

    T rol_k_map[256];
    T rol_k_minus_one_map_rc[256];
    // T rol_k_map_rc[256];
    void init_rol_k_map() noexcept;
    void init_rol_k_minus_one_map_rc() noexcept;
    // void init_rol_k_map_rc() noexcept;

    static inline T rol(T v, size_t n) noexcept;
    static inline T ror(T v, size_t n) noexcept;
};

} // namespace kebab

#endif // NTHASH_H