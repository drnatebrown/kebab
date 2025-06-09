#include <string>
#include <stdio.h>
#include <fstream>
#include <unistd.h>
#include <chrono>
#include <omp.h>

#include "external/kseq.h"
#include "external/CLI11.hpp"
#include "external/hll/hll.h"

#include "kebab/kebab_index.hpp"
#include "kebab/nt_hash.hpp"

#include "constants.hpp"
#include "util.hpp"

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

// Stores sequence information for multi-threaded processing
struct SeqInfo {
    char* seq_content;
    char* seq_name;
    int64_t seq_len;
    int64_t seq_name_len;
    int64_t seq_comment_len;
};

template<typename ProcessFunc>
void process_sequences(kseq_t* seq, uint16_t threads, ProcessFunc process_func) {
    #pragma omp parallel
    {
        SeqInfo seq_info;
        while (true) {
            #pragma omp critical(read_seq)
            {
                seq_info.seq_len = kseq_read(seq);
                if (seq_info.seq_len >= 0) {
                    if (threads > 1) {
                        seq_info.seq_content = strdup(seq->seq.s);
                        seq_info.seq_name = strdup(seq->name.s);  // Copy name for multi-threaded
                    } 
                    // avoid extra copy for single thread
                    else {
                        seq_info.seq_content = seq->seq.s;
                        seq_info.seq_name = seq->name.s;
                    }
                    seq_info.seq_name_len = seq->name.l;
                    seq_info.seq_comment_len = seq->comment.l;
                }
            }
            if (seq_info.seq_len <= 0) {
                break;
            }

            process_func(seq_info);
        }

        if (threads > 1) {
            free(seq_info.seq_content);
            free(seq_info.seq_name);
        }
    }
}

size_t bytes_read(const SeqInfo& seq_info) {
    // seq, header, comment, and newlines
    return seq_info.seq_len + seq_info.seq_name_len + seq_info.seq_comment_len + 2;
}

/* =============================== ESTIMATE =============================== */

uint64_t card_estimate(const std::string& fasta_file, uint16_t kmer_size, KmerMode kmer_mode, uint16_t threads) {
    const auto start_time = std::chrono::steady_clock::now();

    FILE* fp;
    kseq_t* seq = open_fasta(fasta_file, &fp);

    // For progress
    size_t file_size = std::filesystem::file_size(fasta_file);
    size_t bytes_processed = 0;

    kebab::NtManyHash rehasher; // Used only for canonical mode to rehash the value

    hll::hll_t hll(HLL_SIZE);

    auto cardinality_step = [&](const SeqInfo& seq_info) {
        thread_local static kebab::NtHash hasher(kmer_size, use_build_rev_comp(kmer_mode));

        auto add_kmer = [&]() {
            switch (kmer_mode) {
                case KmerMode::FORWARD_ONLY:
                    hll.add(hasher.hash());
                    break;
                case KmerMode::BOTH_STRANDS:
                    hll.add(hasher.hash());
                    hll.add(hasher.hash_rc());
                    break;
                case KmerMode::CANONICAL_ONLY:
                    // hashes again, since canonical biases estimate lower
                    hll.add(rehasher(hasher.hash_canonical()));
                    break;
            }
        };

        hasher.set_sequence(seq_info.seq_content, seq_info.seq_len);
        add_kmer();
        for (size_t i = 1; i < static_cast<size_t>(seq_info.seq_len) - kmer_size + 1; ++i) {
            hasher.unsafe_roll();
            add_kmer();
        }

        #pragma omp critical(update_progress)
        {
            bytes_processed += bytes_read(seq_info);
            std::cerr << "\rEstimating Cardinality: " 
                    << std::fixed << std::setprecision(2) << std::setw(6) 
                    << (bytes_processed * 100.0 / file_size) << "%" << std::flush;
        }
    };

    process_sequences(seq, threads, cardinality_step);
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
    KmerMode kmer_mode = DEFAULT_KMER_MODE;
    double fp_rate = DEFAULT_FP_RATE;
    uint16_t hash_funcs = DEFAULT_HASH_FUNCS;
    uint64_t expected_kmers = DEFAULT_EXPECTED_KMERS;
    uint16_t threads = DEFAULT_BUILD_THREADS;
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
        num_expected_kmers = card_estimate(params.fasta_file, params.kmer_size, params.kmer_mode, params.threads);
    }

    const auto start_time = std::chrono::steady_clock::now();

    Index index(params.kmer_size, num_expected_kmers, params.fp_rate, params.hash_funcs, params.kmer_mode, params.filter_size_mode);

    FILE* fp;
    kseq_t* seq = open_fasta(params.fasta_file, &fp);

    // For progress
    size_t file_size = std::filesystem::file_size(params.fasta_file);
    size_t bytes_processed = 0;

    auto add_sequence_step = [&](const SeqInfo& seq_info) {
        index.add_sequence(seq_info.seq_content, seq_info.seq_len);

        #pragma omp critical(update_progress)
        {
            bytes_processed += bytes_read(seq_info);
            std::cerr << "\rIndexing: " 
                    << std::fixed << std::setprecision(2) << std::setw(6) 
                    << (bytes_processed * 100.0 / file_size) << "%" << std::flush;
        }
    };

    process_sequences(seq, params.threads, add_sequence_step);

    const auto end_time = std::chrono::steady_clock::now();
    const auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time);
    std::cerr << "\rIndexing: 100.00% [" << std::fixed << std::setprecision(2) 
              << (elapsed.count() / 1000.0) << "s]" << std::endl;

    kseq_destroy(seq);
    fclose(fp);

    std::cerr << index.get_stats() << std::endl;

    std::ofstream out(params.output_file + KEBAB_FILE_SUFFIX);
    save_options(out, params);
    index.save(out);
}

void build_index(const BuildParams& params) {
    if (params.fp_rate <= 0 || params.fp_rate >= 1) {
        error_exit("Desired false positive rate (" + std::to_string(params.fp_rate) + ") must be between 0 and 1");
    }
    if (params.hash_funcs > std::size(SEEDS)) {
        error_exit("Number of hashes (" + std::to_string(params.hash_funcs) + ") must be less than the number of seeds (" + std::to_string(std::size(SEEDS)) + ")");
    }

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
    bool prefetch = DEFAULT_PREFETCH;
    uint16_t threads = DEFAULT_SCAN_THREADS;
};

template<typename Index>
void filter_reads(const ScanParams& params, std::ifstream& index_stream) {
    Index index(index_stream);
    if (params.min_mem_length <= index.get_k()) {
        error_exit("min_mem_length (" + std::to_string(params.min_mem_length) + ") must be greater than k (" + std::to_string(index.get_k()) + ")");
    }

    FILE* fp;
    kseq_t* seq = open_fasta(params.fasta_file, &fp);

    // For maximum performance, use system buffer size
    FILE* out = fopen(params.output_file.c_str(), "w");
    int fd = fileno(out);
    long long buffer_size = fpathconf(fd, _PC_REC_XFER_ALIGN);
    if (buffer_size <= 0) {
        buffer_size = DEFAULT_BUFFER_SIZE;
    }
    setvbuf(out, nullptr, _IOFBF, buffer_size);

    auto filter_read_step = [&](const SeqInfo& seq_info) {
        thread_local static std::vector<kebab::Fragment> fragments;

        fragments.clear();
        fragments = index.scan_read(seq_info.seq_content, seq_info.seq_len, params.min_mem_length, params.remove_overlaps, params.prefetch);
        
        if (params.sort_fragments) {
            std::sort(fragments.begin(), fragments.end());
        }
        
        #pragma omp critical(write_fragments)
        {
            for (const auto& fragment : fragments) {
                // use 1-based inclusive
                fprintf(out, ">%s:%zu-%zu\n", seq_info.seq_name, fragment.start + 1, fragment.start + fragment.length);
                fwrite(seq_info.seq_content + fragment.start, 1, fragment.length, out);
                fputc('\n', out);
            }
        }
    };

    process_sequences(seq, params.threads, filter_read_step);

    kseq_destroy(seq);
    fclose(fp);
    fclose(out);
}

void scan_reads(const ScanParams& params) {
    std::ifstream index_stream(params.index_file);
    
    SavedOptions options;
    load_options(index_stream, options);
    
    if (use_shift_filter(options.filter_size_mode)) {
        filter_reads<kebab::KebabIndex<kebab::ShiftFilter>>(params, index_stream);
    } else {
        filter_reads<kebab::KebabIndex<kebab::ModFilter>>(params, index_stream);
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
    build_params.threads = omp_get_max_threads();

    build->add_option("fasta", build_params.fasta_file, "Input FASTA file")->required();
    build->add_option("-o,--output", build_params.output_file, "Output prefix for index file, [PREFIX]" + std::string(KEBAB_FILE_SUFFIX))->required();
    build->add_option("-k,--kmer-size", build_params.kmer_size, "K-mer size used to populate the index")
        ->default_val(DEFAULT_KMER_SIZE)
        ->check(CLI::PositiveNumber);
    build->add_option("--kmer-mode", build_params.kmer_mode, "K-mer strands to include in the index")
        ->default_val(DEFAULT_KMER_MODE)
        ->transform(CLI::CheckedTransformer(std::map<std::string, KmerMode>{
            {"forward", KmerMode::FORWARD_ONLY},
            {"both", KmerMode::BOTH_STRANDS},
            {"canonical", KmerMode::CANONICAL_ONLY}
        }));
    build->add_option("-m,--expected-kmers", build_params.expected_kmers, "Expected number of k-mers (otherwise estimated)")
        ->check(CLI::NonNegativeNumber);
    build->add_option("-e,--fp-rate", build_params.fp_rate, "Desired false positive rate (between 0 and 1)")
        ->default_val(DEFAULT_FP_RATE)
        ->check(CLI::Range(0.0, 1.0))
        ->type_name("FLOAT");
    // build->add_option("-d,--kmer-freq", kmer_freq, "k-mer sampling rate")->default_val(1);
    build->add_option("-f,--hash-funcs", build_params.hash_funcs, "Number of hash functions (otherwise set to minimize index size)")
        ->check(CLI::PositiveNumber);
    build->add_option("-t,--threads", build_params.threads, "Number of threads to use")
        ->default_val(build_params.threads)
        ->check(CLI::PositiveNumber);
    build->add_flag("--no-rounding", no_filter_rounding, "Don't round to power of 2 for filter size (slower)");

    // SCAN COMMAND
    auto scan = app.add_subcommand("scan", "Breaks sequences into fragments using KeBaB index");

    ScanParams scan_params;
    bool no_prefetch = false;
    scan_params.threads = (DEFAULT_PREFETCH) ? omp_get_num_procs() : omp_get_max_threads();

    scan->add_option("fasta", scan_params.fasta_file, "Patterns FASTA file")->required();
    scan->add_option("-i,--index", scan_params.index_file, "KeBaB index file")->required();
    scan->add_option("-o,--output", scan_params.output_file, "Output FASTA file")->required();
    scan->add_option("-l,--mem-length", scan_params.min_mem_length, "Minimum MEM length (must be greater than k-mer size of index)")
        ->default_val(DEFAULT_MIN_MEM_LENGTH)
        ->check(CLI::PositiveNumber);
    scan->add_flag("-s,--sort", scan_params.sort_fragments, "Sort fragments by length");
    scan->add_flag("-r,--remove-overlaps", scan_params.remove_overlaps, "Merge overlapping fragments");
    scan->add_option("-t,--threads", scan_params.threads, "Number of threads to use")
        ->default_val(scan_params.threads)
        ->check(CLI::PositiveNumber);
    scan->add_flag("--no-prefetch", no_prefetch, "Don't prefetch k-mers to avoid latency");

    try {
        app.parse(argc, argv);
        
        if (build->parsed()) { 
            if (no_filter_rounding) {
                build_params.filter_size_mode = FilterSizeMode::EXACT;
            }
            omp_set_num_threads(build_params.threads);
            build_index(build_params);
        }
        if (scan->parsed()) {
            if (no_prefetch) {
                scan_params.prefetch = false;
                scan_params.threads = (scan->count("--threads") > 0) ? scan_params.threads : omp_get_max_threads();
            }
            omp_set_num_threads(scan_params.threads);
            scan_reads(scan_params);
        }

    } catch (const CLI::ParseError &e) {
        return app.exit(e);
    }

    return 0;
}
