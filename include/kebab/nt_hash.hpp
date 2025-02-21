#ifndef NTHASH_H
#define NTHASH_H

#include <cstdint>
#include <cstddef>
#include <climits>

namespace kebab {

template<typename T = uint64_t>
class NtHash {
public:
    explicit NtHash(size_t k) noexcept;
    NtHash() noexcept : k(0), seq(nullptr), len(0), pos(0), hash_val(0) {}

    void set_sequence(const char* seq, size_t len) noexcept;

    [[nodiscard]] size_t get_k() const noexcept { return k; }
    [[nodiscard]] size_t get_pos() const noexcept { return pos; }
    [[nodiscard]] size_t get_len() const noexcept { return len; }

    [[nodiscard]] bool roll() noexcept;
    void unsafe_roll() noexcept;
    [[nodiscard]] T hash() const noexcept { return hash_val; }

private:
    size_t k;

    const char* seq;
    size_t len;
    size_t pos;

    T hash_val;

    T rol_k_map[256];
    void init_rol_k_map() noexcept;

    static inline T rol(T v, size_t n) noexcept;
};

} // namespace kebab

#endif // NTHASH_H