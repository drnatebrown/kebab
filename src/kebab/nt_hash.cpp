/* Implemented from the paper/codebase of https://github.com/bcgsc/ntHash */
#include "kebab/nt_hash.hpp"

namespace {
template<typename T>
struct NtMap;

template<>
struct NtMap<uint64_t> {
    static constexpr uint64_t A = 0x668C9689C1A9287CULL;
    static constexpr uint64_t C = 0x3260979910886E71ULL;
    static constexpr uint64_t G = 0x5BCAA0C13EE6F2BDULL;
    static constexpr uint64_t T = 0x93619763BF5F2651ULL;

    static constexpr uint64_t map[256] = {
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,A,0,C,0,0,0,G,0,0,0,0,0,0,0,0, 
        0,0,0,0,T,0,0,0,0,0,0,0,0,0,0,0,  
        0,A,0,C,0,0,0,G,0,0,0,0,0,0,0,0,  
        0,0,0,0,T,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    };

    static constexpr uint64_t map_rc[256] = {
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,T,0,G,0,0,0,C,0,0,0,0,0,0,0,0, 
        0,0,0,0,A,0,0,0,0,0,0,0,0,0,0,0,  
        0,T,0,G,0,0,0,C,0,0,0,0,0,0,0,0,  
        0,0,0,0,A,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    };   
};

template<>
struct NtMap<uint32_t> {
    static constexpr uint32_t A = 0xC1A9287C;
    static constexpr uint32_t C = 0x10886E71;
    static constexpr uint32_t G = 0x3EE6F2BD;
    static constexpr uint32_t T = 0xBF5F2651;

    static constexpr uint32_t map[256] = {
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,A,0,C,0,0,0,G,0,0,0,0,0,0,0,0, 
        0,0,0,0,T,0,0,0,0,0,0,0,0,0,0,0,  
        0,A,0,C,0,0,0,G,0,0,0,0,0,0,0,0,  
        0,0,0,0,T,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    };

    static constexpr uint32_t map_rc[256] = {
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,  
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, 
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,T,0,G,0,0,0,C,0,0,0,0,0,0,0,0, 
        0,0,0,0,A,0,0,0,0,0,0,0,0,0,0,0,  
        0,T,0,G,0,0,0,C,0,0,0,0,0,0,0,0,  
        0,0,0,0,A,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0
    };
};

template<typename T>
constexpr size_t BITS_IN_TYPE = CHAR_BIT * sizeof(T);
}

namespace kebab {

template<typename T>
NtHash<T>::NtHash(size_t k, bool rev_comp) noexcept
    : k(k)
    , rev_comp(rev_comp)
    , seq(nullptr)
    , len(0)
    , pos(0)
    , hash_val(0)
    , hash_val_rc(0)
    , rol_k_map{}
    , rol_k_map_rc{}
{
    init_rol_k_map();
    init_rol_k_map_rc();
}

template<typename T>
void NtHash<T>::set_sequence(const char* seq, size_t len) noexcept {
    this->seq = seq;
    this->len = len;
    this->pos = 0;
    this->hash_val = 0;

    // Initial hash
    if (len >= k) {
        // Except last, see below
        for (size_t i = 0; i < k - 1; ++i) {
            hash_val ^= rol(NtMap<T>::map[static_cast<uint8_t>(seq[i])], k - i - 1);
        }
        hash_val ^= NtMap<T>::map[static_cast<uint8_t>(seq[k - 1])]; // rol 0 = no shift

        if (rev_comp) {
            this->hash_val_rc = 0;
            hash_val_rc ^= NtMap<T>::map_rc[static_cast<uint8_t>(seq[0])]; // rol 0 = no shift
            for (size_t i = 1; i < k; ++i) {
                hash_val_rc ^= rol(NtMap<T>::map_rc[static_cast<uint8_t>(seq[i])], i);
            }
        }
    }
    else {
        pos = len;
    }
}

template<typename T>
bool NtHash<T>::roll() noexcept {
    if (pos + k >= len) {
        return false;
    }
    unsafe_roll();
    return true;
}

template<typename T>
void NtHash<T>::unsafe_roll() noexcept {
    uint8_t outgoing = static_cast<uint8_t>(seq[pos]);
    uint8_t incoming = static_cast<uint8_t>(seq[pos + k]);
    
    hash_val = rol(hash_val, 1);
    hash_val ^= rol_k_map[outgoing];
    hash_val ^= NtMap<T>::map[incoming];

    if (rev_comp) {
        hash_val_rc ^= NtMap<T>::map_rc[outgoing];
        hash_val_rc ^= rol_k_map_rc[incoming];
        hash_val_rc = ror(hash_val_rc, 1);
    }
    
    ++pos;
}

template<typename T>
void NtHash<T>::init_rol_k_map() noexcept {
    rol_k_map['A'] = rol(NtMap<T>::map['A'], k);
    rol_k_map['C'] = rol(NtMap<T>::map['C'], k);
    rol_k_map['G'] = rol(NtMap<T>::map['G'], k);
    rol_k_map['T'] = rol(NtMap<T>::map['T'], k);
    rol_k_map['a'] = rol_k_map['A'];
    rol_k_map['c'] = rol_k_map['C'];
    rol_k_map['g'] = rol_k_map['G'];
    rol_k_map['t'] = rol_k_map['T'];
}

template<typename T>
void NtHash<T>::init_rol_k_map_rc() noexcept {
    rol_k_map_rc['A'] = rol(NtMap<T>::map_rc['A'], k);
    rol_k_map_rc['C'] = rol(NtMap<T>::map_rc['C'], k);
    rol_k_map_rc['G'] = rol(NtMap<T>::map_rc['G'], k);
    rol_k_map_rc['T'] = rol(NtMap<T>::map_rc['T'], k);
    rol_k_map_rc['a'] = rol_k_map_rc['A'];
    rol_k_map_rc['c'] = rol_k_map_rc['C'];
    rol_k_map_rc['g'] = rol_k_map_rc['G'];
    rol_k_map_rc['t'] = rol_k_map_rc['T'];
}

template<typename T>
T NtHash<T>::rol(T v, size_t n) noexcept {
    return (v << n) | (v >> (BITS_IN_TYPE<T> - n));
}

template<typename T>
T NtHash<T>::ror(T v, size_t n) noexcept {
    return (v >> n) | (v << (BITS_IN_TYPE<T> - n));
}

// Explicit instantiation
template class NtHash<uint64_t>;
template class NtHash<uint32_t>;

} // namespace kebab
