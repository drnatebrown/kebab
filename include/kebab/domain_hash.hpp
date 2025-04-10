#ifndef KEBAB_HASH_HPP
#define KEBAB_HASH_HPP

#include <cstdint>
#include <cmath>
#include <stdexcept>
#include <variant>

namespace kebab {

// =============================================
// Hash and Reducer Types
// =============================================

enum class HashType {
    MULTIPLY,
    MURMUR
};

enum class ReducerType {
    SHIFT,
    MODULO
};

// =============================================
// Base Interfaces
// =============================================

class HashFunction {
public:
    virtual uint64_t operator()(uint64_t x, uint64_t seed) const = 0;
};

class DomainReducer {
public:
    virtual uint64_t reduce(uint64_t hash) const = 0;
};

// =============================================
// Hash Function Implementations
// =============================================

class MultiplyHash : public HashFunction {
public:
    uint64_t operator()(uint64_t x, uint64_t seed) const override {
        return x * seed;
    }
};

class NtManyHash : public HashFunction {
public:
    uint64_t operator()(uint64_t x, uint64_t seed) const override {
        return x * seed;
    }
};

class MurmurHash2 : public HashFunction {
public:
    uint64_t operator()(uint64_t x, uint64_t seed) const override {
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
    static constexpr int r = 47;
    static constexpr int len = 8;
};

// =============================================
// Domain Reducer Implementations
// =============================================

class ShiftReducer : public DomainReducer {
public:
    explicit ShiftReducer(size_t domain_size) 
        : shift(64 - std::floor(std::log2(domain_size))) {}
    
    uint64_t reduce(uint64_t hash) const override {
        return hash >> shift;
    }
private:
    uint32_t shift;
};

class ModuloReducer : public DomainReducer {
public:
    explicit ModuloReducer(size_t domain_size) : domain_size(domain_size) {}
    
    uint64_t reduce(uint64_t hash) const override {
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
    DomainHashFunction(size_t domain_size)
        : hash_(Hash())
        , reducer_(Reducer(domain_size)) {}

    uint64_t operator()(uint64_t x, uint64_t seed) const {
        return reducer_.reduce(hash_(x, seed));
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
using MurmurShift = DomainHashFunction<MurmurHash2, ShiftReducer>;
using MurmurMod = DomainHashFunction<MurmurHash2, ModuloReducer>;

// AES-NI Hash? TODO
} // namespace kebab

#endif // KEBAB_HASH_HPP