#include <string>
#include <stdio.h>
#include <fstream>
#include <unistd.h>
#include <chrono>

#include "external/kseq.h"
#include "external/CLI11.hpp"
#include "external/hll/hll.h"

#include "kebab/kebab_index.hpp"
#include "kebab/nt_hash.hpp"

#include "constants.hpp"

/* =============================== UTILITIES =============================== */

KSEQ_INIT(int, read)

kseq_t* open_fasta(const std::string& fasta_file, FILE** fp) {
    *fp = fopen(fasta_file.c_str(), "r");
    if (!*fp) {
        if (errno == ENOENT) {
            std::cerr << "File not found: " << fasta_file << std::endl;
        } else {
            std::cerr << "Error opening file " << fasta_file << ": " 
                      << strerror(errno) << std::endl;
        }
        return nullptr;
    }
    return kseq_init(fileno(*fp));
}

/* =============================== ESTIMATE =============================== */

uint64_t card_estimate(const std::string& fasta_file, uint16_t kmer_size, KmerMode kmer_mode) {
    const auto start_time = std::chrono::steady_clock::now();

    FILE* fp;
    kseq_t* seq = open_fasta(fasta_file, &fp);

    // For progress
    size_t file_size = std::filesystem::file_size(fasta_file);
    size_t bytes_processed = 0;

    kebab::NtHash hasher(kmer_size, use_build_rev_comp(kmer_mode));
    kebab::NtManyHash rehasher; // Used only for canonical mode to rehash the value

    hll::hll_t hll(HLL_SIZE);

    int64_t l = 0;
    while ((l = kseq_read(seq)) >= 0) {
        hasher.set_sequence(seq->seq.s, l);
        for (size_t i = 0; i < static_cast<size_t>(l) - kmer_size + 1; ++i) {
            switch (kmer_mode) {
                case KmerMode::FORWARD_ONLY:
                    hll.add(hasher.hash());
                    break;
                case KmerMode::BOTH_STRANDS:
                    hll.add(hasher.hash());
                    hll.add(hasher.hash_rc());
                    break;
                case KmerMode::CANONICAL_ONLY:
                    hll.add(rehasher(hasher.hash_canonical())); // hashes again, since canonical biases estimate lower
                    break;
            }
            hasher.unsafe_roll();
        }

        bytes_processed += l + seq->name.l + seq->comment.l + 2;  // header, seq, and newlines
        std::cerr << "\rEstimating Cardinality: " 
                  << std::fixed << std::setprecision(2) << std::setw(6) 
                  << (bytes_processed * 100.0 / file_size) << "%" << std::flush;
    }
    const auto end_time = std::chrono::steady_clock::now();
    const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);   
    std::cerr << "\rEstimating Cardinality: 100.00% [" << std::fixed << std::setprecision(2) 
              << (elapsed.count() / 1000.0) << "s]" << std::endl;

    // TODO: ADD STATS
    std::cerr << "\tEstimate: " << static_cast<uint64_t>(std::ceil(hll.report())) << std::endl;
    std::cerr << "\tError Bounds: " << hll.est_err() << std::endl;

    kseq_destroy(seq);
    fclose(fp);

    return static_cast<uint64_t>(std::ceil(hll.report()));
}

/* =============================== BUILD =============================== */

struct BuildParams {
    std::string fasta_file;
    std::string output_file;
    uint16_t kmer_size = DEFAULT_KMER_SIZE;
    double fp_rate = DEFAULT_FP_RATE;
    uint16_t hash_funcs = DEFAULT_HASH_FUNCS;
    uint64_t expected_kmers = DEFAULT_EXPECTED_KMERS;
    KmerMode kmer_mode = DEFAULT_KMER_MODE;
    FilterSizeMode filter_size_mode = DEFAULT_FILTER_SIZE_MODE;
};

struct SavedOptions {
    FilterSizeMode filter_size_mode = DEFAULT_FILTER_SIZE_MODE;
};

void save_options(std::ostream& out, const BuildParams& params) {
    out.write(reinterpret_cast<const char*>(&params.filter_size_mode), sizeof(params.filter_size_mode));
}

void load_options(std::istream& in, SavedOptions& options) {
    in.read(reinterpret_cast<char*>(&options.filter_size_mode), sizeof(options.filter_size_mode));
}

template<typename Index>
void populate_index(const BuildParams& params) {
    uint64_t num_expected_kmers = params.expected_kmers;
    if (num_expected_kmers == 0) {
        num_expected_kmers = card_estimate(params.fasta_file, params.kmer_size, params.kmer_mode);
    }

    const auto start_time = std::chrono::steady_clock::now();

    Index index(params.kmer_size, num_expected_kmers, params.fp_rate, params.hash_funcs, params.kmer_mode, params.filter_size_mode);

    FILE* fp;
    kseq_t* seq = open_fasta(params.fasta_file, &fp);

    // For progress
    size_t file_size = std::filesystem::file_size(params.fasta_file);
    size_t bytes_processed = 0;

    int64_t l = 0;
    while ((l = kseq_read(seq)) >= 0) {
        index.add_sequence(seq->seq.s, l);

        bytes_processed += l + seq->name.l + seq->comment.l + 2;  // header, seq, and newlines
        std::cerr << "\rIndexing: " 
                  << std::fixed << std::setprecision(2) << std::setw(6) 
                  << (bytes_processed * 100.0 / file_size) << "%" << std::flush;
    }
    const auto end_time = std::chrono::steady_clock::now();
    const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cerr << "\rIndexing: 100.00% [" << std::fixed << std::setprecision(2) 
              << (elapsed.count() / 1000.0) << "s]" << std::endl;

    kseq_destroy(seq);
    fclose(fp);

    std::cerr << index.get_stats() << std::endl;

    std::ofstream out(params.output_file + FILE_EXTENSION);
    save_options(out, params);
    index.save(out);
}

void build_index(const BuildParams& params) {
    if (use_shift_filter(params.filter_size_mode)) {
        populate_index<kebab::KebabIndex<kebab::ShiftFilter>>(params);
    } else {
        populate_index<kebab::KebabIndex<kebab::ModFilter>>(params);
    }
}

/* =============================== SCAN =============================== */

struct ScanParams {
    std::string fasta_file;
    std::string index_file;
    std::string output_file;
    uint64_t min_mem_length = DEFAULT_MIN_MEM_LENGTH;
    bool sort_fragments = DEFAULT_SORT_FRAGMENTS;
    bool remove_overlaps = DEFAULT_REMOVE_OVERLAPS;
};

template<typename Index>
void process_reads(const ScanParams& params, std::ifstream& index_stream) {
    Index index(index_stream);

    FILE* fp;
    kseq_t* seq = open_fasta(params.fasta_file, &fp);
    int64_t l = 0;

    // For maximum performance, use system buffer size
    FILE* out = fopen(params.output_file.c_str(), "w");
    int fd = fileno(out);
    long long buffer_size = fpathconf(fd, _PC_REC_XFER_ALIGN);
    if (buffer_size <= 0) {
        buffer_size = DEFAULT_BUFFER_SIZE;
    }
    setvbuf(out, nullptr, _IOFBF, buffer_size);

    while ((l = kseq_read(seq)) >= 0) {
        std::vector<kebab::Fragment> fragments = index.scan_read(seq->seq.s, l, params.min_mem_length, params.remove_overlaps);

        if (params.sort_fragments) {
            std::sort(fragments.begin(), fragments.end());
        }

        for (const auto& fragment : fragments) {
            // use 1-based inclusive
            fprintf(out, ">%s:%zu-%zu\n", seq->name.s, fragment.start + 1, fragment.start + fragment.length);
            fwrite(seq->seq.s + fragment.start, 1, fragment.length, out);
            fputc('\n', out);
        }
    }

    kseq_destroy(seq);
    fclose(fp);
    fclose(out);
}

void scan_reads(const ScanParams& params) {
    std::ifstream index_stream(params.index_file);
    
    SavedOptions options;
    load_options(index_stream, options);
    
    if (use_shift_filter(options.filter_size_mode)) {
        process_reads<kebab::KebabIndex<kebab::ShiftFilter>>(params, index_stream);
    } else {
        process_reads<kebab::KebabIndex<kebab::ModFilter>>(params, index_stream);
    }
}

/* =============================== MAIN =============================== */

int main(int argc, char** argv) {
    CLI::App app{"KeBaB: K-mer Based Breaking"};
    app.require_subcommand(1);
    app.set_version_flag("--version", "KeBaB " + std::string(VERSION));
    app.failure_message(CLI::FailureMessage::help);

    // BUILD COMMAND
    auto build = app.add_subcommand("build", "Build a KeBaB index");

    BuildParams build_params;
    bool no_filter_rounding = false;

    build->add_option("fasta", build_params.fasta_file, "Input FASTA file")->required();
    build->add_option("-o,--output", build_params.output_file, "Output prefix for .kbb index file")->required();
    build->add_option("-k,--kmer-size", build_params.kmer_size, "k-mer size")
        ->default_val(DEFAULT_KMER_SIZE)
        ->check(CLI::PositiveNumber);
    build->add_option("-m,--expected-kmers", build_params.expected_kmers, "Expected number of k-mers (if not provided, will be estimated)")
        ->check(CLI::NonNegativeNumber)
        ->default_val(DEFAULT_EXPECTED_KMERS);
    build->add_option("-e,--fp-rate", build_params.fp_rate, "False positive rate (between 0 and 1)")
        ->default_val(DEFAULT_FP_RATE)
        ->check(CLI::Range(0.0, 1.0))
        ->type_name("FLOAT");
    // build->add_option("-d,--kmer-freq", kmer_freq, "k-mer sampling rate")->default_val(1);
    build->add_option("-f,--hash-funcs", build_params.hash_funcs, "Number of hash functions")
        ->check(CLI::PositiveNumber);
    build->add_option("--kmer-mode", build_params.kmer_mode, "K-mer strand mode")
        ->default_val(DEFAULT_KMER_MODE)
        ->transform(CLI::CheckedTransformer(std::map<std::string, KmerMode>{
            {"forward", KmerMode::FORWARD_ONLY},
            {"both", KmerMode::BOTH_STRANDS},
            {"canonical", KmerMode::CANONICAL_ONLY}
        }));
    build->add_flag("--no-rounding", no_filter_rounding, "Don't round to power of 2 for filter size (slower)");

    // SCAN COMMAND
    auto scan = app.add_subcommand("scan", "Breaks sequences into fragments using KeBaB index");

    ScanParams scan_params;

    scan->add_option("fasta", scan_params.fasta_file, "Patterns FASTA file")->required();
    scan->add_option("-i,--index", scan_params.index_file, "KeBaB index file")->required();
    scan->add_option("-o,--output", scan_params.output_file, "Output FASTA file")->required();
    scan->add_option("-l,--mem-length", scan_params.min_mem_length, "Minimum MEM length")
        ->default_val(DEFAULT_MIN_MEM_LENGTH)
        ->check(CLI::PositiveNumber);
    scan->add_flag("-s,--sort", scan_params.sort_fragments, "Sort fragments");
    scan->add_flag("-r,--remove-overlaps", scan_params.remove_overlaps, "Merge overlapping fragments");

    try {
        app.parse(argc, argv);
        
        if (build->parsed()) { 
            if (no_filter_rounding) {
                build_params.filter_size_mode = FilterSizeMode::EXACT;
            }
            build_index(build_params);
        }
        if (scan->parsed()) {
            scan_reads(scan_params);
        }

    } catch (const CLI::ParseError &e) {
        return app.exit(e);
    }

    return 0;
}