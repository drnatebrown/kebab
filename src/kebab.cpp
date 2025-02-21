#include <string>
#include <stdio.h>
#include <fstream>
#include <unistd.h>

#include "external/kseq.h"
#include "external/CLI11.hpp"
#include "kebab/kebab_index.hpp"

static constexpr size_t DEFAULT_BUFFER_SIZE = 64ULL * 1024ULL * 1024ULL; // 64MB
KSEQ_INIT(int, read)

constexpr uint16_t DEFAULT_KMER_SIZE = 20;
constexpr double DEFAULT_FP_RATE = 0.01;
constexpr uint16_t DEFAULT_HASH_FUNCS = 0;

struct BuildParams {
    std::string fasta_file;
    std::string output_file;
    uint16_t kmer_size = DEFAULT_KMER_SIZE;
    double fp_rate = DEFAULT_FP_RATE;
    uint16_t hash_funcs = DEFAULT_HASH_FUNCS;
    uint64_t expected_kmers;
};

void build_index(const BuildParams& params) {
    // temp for testing
    // uint64_t expected_kmers = 445759753ULL;
    kebab::KebabIndex index(params.kmer_size, params.expected_kmers, params.fp_rate, params.hash_funcs);

    FILE* fp = fopen(params.fasta_file.c_str(), "r");
    if (!fp) {
        if (errno == ENOENT) {
            std::cerr << "File not found: " << params.fasta_file << std::endl;
        } else {
            std::cerr << "Error opening file " << params.fasta_file << ": " 
                      << strerror(errno) << std::endl;
        }
        return;
    }
    kseq_t* seq = kseq_init(fileno(fp));

    // For progress
    size_t file_size = std::filesystem::file_size(params.fasta_file);
    size_t bytes_processed = 0;

    int64_t l = 0;
    while ((l = kseq_read(seq)) >= 0) {
        index.add_sequence(seq->seq.s, l);

        bytes_processed += l + seq->name.l + seq->comment.l + 2;  // header, seq, and newlines
        std::cerr << "\rFASTA Processed: " 
                  << std::fixed << std::setprecision(2) << std::setw(6) 
                  << (bytes_processed * 100.0 / file_size) << "%" << std::flush;
    }
    std::cerr << "\rFASTA Processed: 100.00%" << std::endl;

    kseq_destroy(seq);
    fclose(fp);

    std::cerr << index.get_stats() << std::endl;

    std::ofstream out(params.output_file);
    index.save(out);
}

constexpr uint64_t DEFAULT_MIN_MEM_LENGTH = 20;
constexpr bool DEFAULT_SORT_FRAGMENTS = false;

struct ScanParams {
    std::string fasta_file;
    std::string index_file;
    std::string output_file;
    bool sort_fragments = DEFAULT_SORT_FRAGMENTS;
    uint64_t min_mem_length = DEFAULT_MIN_MEM_LENGTH;
};

void scan_reads(const ScanParams& params) {
    std::ifstream index_stream(params.index_file);
    kebab::KebabIndex index(index_stream);

    FILE* fp = fopen(params.fasta_file.c_str(), "r");
    if (!fp) {
        if (errno == ENOENT) {
            std::cerr << "File not found: " << params.fasta_file << std::endl;
        } else {
            std::cerr << "Error opening file " << params.fasta_file << ": " 
                      << strerror(errno) << std::endl;
        }
        return;
    }
    kseq_t* seq = kseq_init(fileno(fp));
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
        std::vector<kebab::Fragment> fragments = index.scan_read(seq->seq.s, l, params.min_mem_length);

        if (params.sort_fragments) {
            std::sort(fragments.begin(), fragments.end());
        }

        for (const auto& fragment : fragments) {
            fprintf(out, ">%s:[%zu-%zu]\n", seq->name.s, fragment.start, fragment.start + fragment.length - 1);
            fwrite(seq->seq.s + fragment.start, 1, fragment.length, out);
            fputc('\n', out);
        }
    }

    kseq_destroy(seq);
    fclose(fp);
    fclose(out);
}

int main(int argc, char** argv) {
    CLI::App app{"KeBaB: K-mer Based Breaking"};
    app.require_subcommand(1);
    app.set_version_flag("--version", "KeBaB 0.1");
    app.failure_message(CLI::FailureMessage::help);

    // BUILD COMMAND
    auto build = app.add_subcommand("build", "Build a KeBaB index");

    // std::string fasta_file;
    // std::string build_output;
    // uint16_t kmer_size = 20;
    // double fp_rate = 0.01;
    // // uint16_t kmer_freq = 1;
    // uint16_t hash_funcs = 0;  // 0 indicates not specified

    BuildParams build_params;

    build->add_option("fasta", build_params.fasta_file, "Input FASTA file")->required();
    build->add_option("-o,--output", build_params.output_file, "Output prefix for .kbb index file")->required();

    // TO REMOVE
    build->add_option("-m,--expected-kmers", build_params.expected_kmers, "Expected number of k-mers")
        ->check(CLI::PositiveNumber)
        ->required();

    build->add_option("-k,--kmer-size", build_params.kmer_size, "k-mer size")
        ->default_val(DEFAULT_KMER_SIZE)
        ->check(CLI::PositiveNumber);
    build->add_option("-e,--fp-rate", build_params.fp_rate, "False positive rate (between 0 and 1)")
        ->default_val(DEFAULT_FP_RATE)
        ->check(CLI::Range(0.0, 1.0))
        ->type_name("FLOAT");
    // build->add_option("-d,--kmer-freq", kmer_freq, "k-mer sampling rate")->default_val(1);
    build->add_option("-f,--hash-funcs", build_params.hash_funcs, "Number of hash functions")
        ->check(CLI::PositiveNumber);

    // SCAN COMMAND
    auto scan = app.add_subcommand("scan", "Breaks sequences into fragments using KeBaB index");

    // std::string read_file;
    // std::string index_file;
    // std::string scan_output;
    // bool sort_fragments = false;
    // uint64_t min_mem_length = 20;

    ScanParams scan_params;

    scan->add_option("fasta", scan_params.fasta_file, "Patterns FASTA file")->required();
    scan->add_option("-i,--index", scan_params.index_file, "KeBaB index file")->required();
    scan->add_option("-o,--output", scan_params.output_file, "Output FASTA file")->required();
    scan->add_flag("-s,--sort", scan_params.sort_fragments, "Sort fragments")->default_val(DEFAULT_SORT_FRAGMENTS);
    scan->add_option("-l,--mem-length", scan_params.min_mem_length, "Minimum MEM length")
        ->default_val(DEFAULT_MIN_MEM_LENGTH)
        ->check(CLI::PositiveNumber);

    try {
        app.parse(argc, argv);
        
        if (build->parsed()) {
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