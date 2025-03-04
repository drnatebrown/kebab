#ifndef KEBAB_HASH_HPP
#define KEBAB_HASH_HPP

#include <cstdint>
#include <cmath>

namespace kebab {

class DomainHashFunction {
public:
    DomainHashFunction() : domain_size(0) {}
    explicit DomainHashFunction(size_t domain_size) : domain_size(domain_size) {}
    virtual uint64_t operator()(uint64_t x, uint64_t seed) const = 0;
protected:
    size_t domain_size;
};

// Shifts into range [0,m), best for power of 2 domains
class MultiplyShift : public DomainHashFunction {  
public:
    MultiplyShift() : DomainHashFunction(0), shift(0) {}
    explicit MultiplyShift(size_t domain_size) 
        : DomainHashFunction(domain_size)
        , shift(64 - std::floor(std::log2(domain_size))) {}

    uint64_t operator()(uint64_t x, uint64_t seed) const override {
        return ((x * seed) >> shift);
    }
private:
    uint32_t shift;
};

// Multiplies and takes modulo
class MultiplyModulo : public DomainHashFunction {
public:
    MultiplyShiftModulo() : DomainHashFunction(0) {}
    explicit MultiplyShiftModulo(size_t domain_size) : DomainHashFunction(domain_size) {}

    uint64_t operator()(uint64_t x, uint64_t seed) const override {
        return (x * seed) % domain_size;
    }
};

// MurmurHash2 optimized for 64-bit systems
struct MurmurMix2 : public DomainHashFunction {
    MurmurMix2() : DomainHashFunction(0) {}
    explicit MurmurMix2(size_t domain_size) : DomainHashFunction(domain_size) {}

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
        
        return h % domain_size;
    }   
private:
    static constexpr uint64_t m = 0xc6a4a7935bd1e995ULL;
    static constexpr int r = 47;
    static constexpr int len = 8;
};

} // namespace kebab

#endif // KEBAB_HASH_HPP