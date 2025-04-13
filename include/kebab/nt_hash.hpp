#ifndef NTHASH_H
#define NTHASH_H

#include <cstdint>
#include <cstddef>
#include <climits>
#include <algorithm>
#include <mutex>
#include <unordered_map>
#include "constants.hpp"

namespace kebab {

template<typename T = uint64_t>
class NtHash {
public:
    explicit NtHash(size_t k, bool rev_comp = DEFAULT_REVERSE_COMPLEMENT) noexcept;
    NtHash() noexcept : k(0), rev_comp(DEFAULT_REVERSE_COMPLEMENT), seq(nullptr), len(0), pos(0), hash_val(0), hash_val_rc(0) {}

    void set_sequence(const char* seq, size_t len) noexcept;

    [[nodiscard]] size_t get_k() const noexcept { return k; }
    [[nodiscard]] size_t get_pos() const noexcept { return pos; }
    [[nodiscard]] size_t get_len() const noexcept { return len; }

    [[nodiscard]] static T get_max_hash() noexcept { return static_cast<T>(-1); } // maximum hash value for T

    [[nodiscard]] bool roll() noexcept;
    void unsafe_roll() noexcept;
    [[nodiscard]] T hash() const noexcept { return hash_val; }
    [[nodiscard]] T hash_rc() const noexcept { return hash_val_rc; }
    [[nodiscard]] T hash_canonical() const noexcept { return rev_comp ? std::min(hash_val, hash_val_rc) : hash_val; }

private:
    size_t k;
    bool rev_comp;

    const char* seq;
    size_t len;
    size_t pos;

    T hash_val;
    T hash_val_rc;

    static std::mutex init_mutex;
    static std::unordered_map<size_t, std::array<T, 256>> cached_rol_k_map;
    static std::unordered_map<size_t, std::array<T, 256>> cached_rol_k_map_rc;

    std::array<T, 256>* rol_k_map;
    std::array<T, 256>* rol_k_map_rc;
    void init_rol_k_map() noexcept;
    void init_rol_k_map_rc() noexcept;

    static inline T rol(T v, size_t n) noexcept;
    static inline T ror(T v, size_t n) noexcept;
};

} // namespace kebab

#endif // NTHASH_H