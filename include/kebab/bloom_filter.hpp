#ifndef BLOOM_FILTER_H
#define BLOOM_FILTER_H


#include <vector>
#include <cstdint>
#include <climits>
#include <cmath>
#include <iostream>

#include "constants.hpp"

#include "kebab/domain_hash.hpp"

namespace {

// Supports at most 32 hash functions
constexpr uint64_t SEEDS[] = {          
    0x153C67147CEBD9C1, 0xE9E9221977E2486E,
    0xBD2A5DE364F86CEC, 0xF53E63242C7C96CA,
    0xEA71F713607B8025, 0xDA1DC2E81860AC93,
    0x700FC578B9B89EFC, 0x7ED09A9433D0F542,
    0xED43BDEDBCF69432, 0x1D322B028A861DAA,
    0x6E8CDB8F04EE5FFD, 0xEC53221EFD3A5C53,
    0x01EE14F09892D967, 0xD6382ACCCBCF0420,
    0xD448F78598D09FBE, 0x922AA2623D2BF77A,
    0x4AF98D70BD02F4D9, 0xBE9A532696D539D9,
    0x57CB1CF8FA6F105D, 0x4347990C105CF57C,
    0xD5E6B9B31C51D5D6, 0x2196C4CF3D467371,
    0x78BD99C62BA864CD, 0x0B747BD60B9F2FB4,
    0xE636A63B15DC2C60, 0xE3D4C1379D7C2FF0,
    0x2B5C7FAF45C1B370, 0xFE0247B305095328,
    0xE4F3205AADABEA31, 0xD631A450CF4BA7BA,
    0x7E0034EEC6C9E610, 0xCAF71C56BB5D4B4D
};

inline uint64_t previous_power_of_two(uint64_t x) {
    if (x == 0) return 1;
    return 1ull << static_cast<uint64_t>(std::floor(std::log2(static_cast<double>(x))));
}

inline uint64_t next_power_of_two(uint64_t x) {
    if (x == 0) return 1;
    return 1ull << static_cast<uint64_t>(std::ceil(std::log2(static_cast<double>(x))));
}

} // namespace  

namespace kebab {

using word_t = uint64_t;
static constexpr size_t BITS_PER_WORD = sizeof(word_t) * CHAR_BIT;

inline size_t calculate_num_words(size_t size) {
    return std::ceil(static_cast<double>(size) / BITS_PER_WORD);
}

template<typename Hash = MultiplyShift, bool REUSE_FIRST_HASH = DEFAULT_REUSE_FIRST_HASH>
class BloomFilter {
public:
    BloomFilter() : num_elements(0), error_rate(0), bits(0), set_bits(0), filter(), num_hashes(0), hash() {}

    BloomFilter(size_t elements, double error_rate = DEFAULT_FP_RATE, size_t num_hashes = DEFAULT_HASH_FUNCS, FilterSizeMode filter_size_mode = DEFAULT_FILTER_SIZE_MODE) {
        init(elements, error_rate, num_hashes, filter_size_mode);
    }

    void add(uint64_t val) {
        if constexpr (REUSE_FIRST_HASH) {
            set_bit(hash.reduce(val));
            for (size_t i = 1; i < num_hashes; ++i) {
                set_bit(hash(val, SEEDS[i]));
            }
        } else {
            for (size_t i = 0; i < num_hashes; ++i) {
                set_bit(hash(val, SEEDS[i]));
            }
        }
    }

    bool contains(uint64_t val) const {
        if constexpr (REUSE_FIRST_HASH) {
            if (!check_bit(hash.reduce(val))) {
                return false;
            }
            for (size_t i = 1; i < num_hashes; ++i) {
                if (!check_bit(hash(val, SEEDS[i]))) {
                    return false;
                }
            }
        } else {
            for (size_t i = 0; i < num_hashes; ++i) {
                if (!check_bit(hash(val, SEEDS[i]))) {
                    return false;
                }
            }
        }

        return true;
    }

    std::string get_stats() const {
        double load_factor = static_cast<double>(set_bits) / bits;
        return "\tDesired FP Rate: " + std::to_string(error_rate) + "\n"
               "\tObserved FP Rate: " + std::to_string(std::pow(load_factor, num_hashes)) + "\n"
               "\t# Hashes: " + std::to_string(num_hashes) + "\n"
               "\t# Set Bits: " + std::to_string(set_bits) + "\n"
               "\t# Bits: " + std::to_string(bits) + "\n"
               "\tLoad: " + std::to_string(load_factor);
    }
    
    void save(std::ostream& out) const {
        out.write(reinterpret_cast<const char*>(&bits), sizeof(bits));
        out.write(reinterpret_cast<const char*>(&set_bits), sizeof(set_bits));

        out.write(reinterpret_cast<const char*>(filter.data()), filter.size() * sizeof(word_t));

        out.write(reinterpret_cast<const char*>(&num_hashes), sizeof(num_hashes));
    }

    void load(std::istream& in) {
        in.read(reinterpret_cast<char*>(&bits), sizeof(bits));
        in.read(reinterpret_cast<char*>(&set_bits), sizeof(set_bits));

        filter = std::vector<word_t>(calculate_num_words(bits), 0ULL);
        in.read(reinterpret_cast<char*>(filter.data()), filter.size() * sizeof(word_t));

        in.read(reinterpret_cast<char*>(&num_hashes), sizeof(num_hashes));
        hash = Hash(bits);
    }

private:
    size_t num_elements;
    double error_rate;
    
    size_t bits;
    size_t set_bits;
    std::vector<word_t> filter;

    size_t num_hashes;
    Hash hash;

    void init(size_t elements, double error_rate, size_t num_hashes, FilterSizeMode filter_size_mode) {
        num_elements = elements;
        this->error_rate = error_rate;
        validate_params();

        bits = (num_hashes == 0) 
            ? optimal_bits(elements, error_rate) 
            : optimal_bits(elements, error_rate, num_hashes);

        switch (filter_size_mode) {
            case FilterSizeMode::NEXT_POWER_OF_TWO:
                bits = next_power_of_two(bits);
                break;
            case FilterSizeMode::PREVIOUS_POWER_OF_TWO:
                bits = previous_power_of_two(bits);
                break;
            default:
                break;
        }

        set_bits = 0;
        filter = std::vector<word_t>(calculate_num_words(bits), 0ULL);

        this->num_hashes = (num_hashes == 0) ? optimal_hashes(error_rate) : num_hashes;
        hash = Hash(bits);
        validate_num_hashes();
    }

    size_t optimal_hashes(double error_rate) const {
        // k = -ln(p) / ln(2)
        double k = -std::log(error_rate) / std::log(2);

        size_t k_ceil = static_cast<size_t>(std::ceil(k));
        size_t k_floor = static_cast<size_t>(std::floor(k));
        if (k_floor == 0) {
            return k_ceil;
        }

        // fp = (1 - e^(-k * n / m))^k
        auto fp = [this](size_t k) {
            return std::pow(1-std::exp(-k*num_elements/bits), k);
        };

        double p_ceil = fp(k_ceil);
        double p_floor = fp(k_floor);

        return (p_ceil < p_floor) ? k_ceil : k_floor;
    }

    static size_t optimal_bits(size_t elements, double error_rate, size_t num_hashes) {
        // m = (-k * n) / ln(1-p^(1/k))
        return static_cast<size_t>(((-1.0 * num_hashes * elements) / (std::log(1-std::pow(error_rate, 1.0/num_hashes)))));
    }

    static size_t optimal_bits(size_t elements, double error_rate) {
        // m = (-n ln(p)) / (ln(2))^2
        return static_cast<size_t>((-1.0 *elements * std::log(error_rate)) / (std::log(2) * std::log(2)));
    } 

    void validate_params() const {
        if (this->error_rate <= 0 || this->error_rate >= 1) {
            throw std::invalid_argument("Error rate must be between 0 and 1, not " + std::to_string(this->error_rate));
        }
        if (this->num_elements == 0) {
            throw std::invalid_argument("Estimated number of elements must be greater than 0, not " + std::to_string(this->num_elements));
        }
    }

    void validate_num_hashes() const {
        if (this->num_hashes >= std::size(SEEDS)) {
            throw std::invalid_argument("Number of hashes must be less than the number of seeds (" + std::to_string(std::size(SEEDS)) + "), not " + std::to_string(this->num_hashes));
        }
    }

    word_t* get_word(uint64_t hash_val) noexcept {
        return &filter[hash_val / BITS_PER_WORD];
    }

    word_t get_word(uint64_t hash_val) const noexcept {
        return filter[hash_val / BITS_PER_WORD];
    }

    static constexpr word_t get_bit_mask(uint64_t hash_val) noexcept {
        return word_t{1} << (hash_val % BITS_PER_WORD);
    }

    void set_bit(uint64_t hash_val) {
        word_t* word = get_word(hash_val);
        word_t bit_mask = get_bit_mask(hash_val);

        if (!(*word & bit_mask)) {
            *word |= bit_mask;
            set_bits++;
        }
    }

    bool check_bit(uint64_t hash_val) const {
        word_t word = get_word(hash_val);
        word_t bit_mask = get_bit_mask(hash_val);
        return word & bit_mask;
    }
};

using ModFilter = BloomFilter<MultiplyMod>;
using ShiftFilter = BloomFilter<MultiplyShift>;

} // namespace kebab

#endif // BLOOM_FILTER_H

