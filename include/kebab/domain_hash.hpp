#ifndef KEBAB_HASH_HPP
#define KEBAB_HASH_HPP

#include <cstdint>
#include <cmath>

namespace kebab {

// =============================================
// Hash Functions
// =============================================

class MultiplyHash {
public:
    uint64_t operator()(uint64_t x, uint64_t seed) const {
        return x * seed;
    }
};

class NtManyHash {
public:
    uint64_t operator()(uint64_t x, uint64_t seed) const {
        x *= seed;
        x ^= x >> shift;
        return x;
    }

    uint64_t operator()(uint64_t x) const {
        return operator()(x, seed);
    }
private:
    static constexpr uint8_t shift = 27;
    static constexpr uint64_t seed = 0x90b45d39fb6da1fa;
};

class MurmurHash2 {
public:
    uint64_t operator()(uint64_t x, uint64_t seed) const {
        uint64_t h = seed ^ (len * m);
        uint64_t k = x;
        
        k *= m;
        k ^= k >> r;
        k *= m;
        
        h ^= k;
        h *= m;
        
        h ^= h >> r;
        h *= m;
        h ^= h >> r;
        
        return h;
    }   
private:
    static constexpr uint64_t m = 0xc6a4a7935bd1e995ULL;
    static constexpr uint8_t r = 47;
    static constexpr uint8_t len = 8;
};

// AES-NI Hash? TODO

// =============================================
// Domain Reducers
// =============================================

class ShiftReducer {
public:
    explicit ShiftReducer(size_t domain_size) 
        : shift(64 - std::floor(std::log2(domain_size))) {}
    
    uint64_t operator()(uint64_t hash) const {
        return hash >> shift;
    }
private:
    uint8_t shift;
};

class ModuloReducer {
public:
    explicit ModuloReducer(size_t domain_size) : domain_size(domain_size) {}
    
    uint64_t operator()(uint64_t hash) const {
        return hash % domain_size;
    }
private:
    size_t domain_size;
};

// =============================================
// Hash Function + Domain Reducer Combination
// =============================================

template<typename Hash=MultiplyHash, typename Reducer=ShiftReducer>
class DomainHashFunction {
public:
    DomainHashFunction()
        : hash_(Hash())
        , reducer_(Reducer(0)) {}

    DomainHashFunction(size_t domain_size)
        : hash_(Hash())
        , reducer_(Reducer(domain_size)) {}

    uint64_t operator()(uint64_t x, uint64_t seed) const {
        return reducer_(hash_(x, seed));
    }

    uint64_t hash(uint64_t x, uint64_t seed) const {
        return hash_(x, seed);
    }

    uint64_t reduce(uint64_t hash) const {
        return reducer_(hash);
    }

private:
    Hash hash_;
    Reducer reducer_;
};

// =============================================
// Common Hash Combinations
// =============================================

using MultiplyShift = DomainHashFunction<MultiplyHash, ShiftReducer>;
using MultiplyMod = DomainHashFunction<MultiplyHash, ModuloReducer>;
using NtManyShift = DomainHashFunction<NtManyHash, ShiftReducer>;
using NtManyMod = DomainHashFunction<NtManyHash, ModuloReducer>;
using MurmurShift = DomainHashFunction<MurmurHash2, ShiftReducer>;
using MurmurMod = DomainHashFunction<MurmurHash2, ModuloReducer>;

} // namespace kebab

#endif // KEBAB_HASH_HPP