![Image](https://github.com/user-attachments/assets/ceccfb5b-b557-435b-b009-1aa7f3424bb2)  
_**K**-m**e**r **Ba**sed **B**reaking for finding long maximal exact matches._  
_(Release Version 1.0.1)_
  
**KeBaB** breaks nucleotide DNA patterns into pseudo-MEM fragments, filtering out subsequences which cannot overlap any maximal exact match (MEM) of some minimum match length, $L$. Using fragments results in faster MEM-finding queries. The full paper can be found on [arXiv](https://arxiv.org/abs/2502.20338).
  
![Image](https://github.com/user-attachments/assets/2443f05e-e1f1-4ffe-a58c-367071a924a7)
  
## How-to
### Compile
Creates `./kebab` executable:
```
git clone https://github.com/drnatebrown/kebab.git
cd kebab
make
```
### Build
Build a KeBaB index (bloom filter).
```
Usage: ./kebab build [OPTIONS] fasta

Positionals:
  fasta TEXT REQUIRED         Input FASTA file

Options:
  -h,--help                   Print this help message and exit
  -o,--output TEXT REQUIRED   Output prefix for index file, [PREFIX].kbb
  -k,--kmer-size UINT:POSITIVE [20] 
                              K-mer size used to populate the index
  --kmer-mode ENUM:value in {both->0,canonical->1,forward->2} OR {0,1,2} [1] 
                              K-mer strands to include in the index
  -m,--expected-kmers UINT:POSITIVE
                              Expected number of k-mers (otherwise estimated)
  -e,--fp-rate FLOAT:FLOAT in [0 - 1] [0.1] 
                              Desired false positive rate (between 0 and 1)
  -f,--hash-funcs UINT:POSITIVE
                              Number of hash functions (otherwise set to minimize index size)
  -t,--threads UINT:POSITIVE [8] 
                              Number of threads to use
  --no-rounding               Don't round to power of 2 for filter size (slower)
```
Note that a chosen ``-k`` affects which minimum MEM lengths are valid (see below).
### Scan
Breaks sequences into fragments using KeBaB index. Fragments use ``[SEQ]:[START]-[END]`` notation where the range is 1-based and inclusive.
```
Usage: ./kebab scan [OPTIONS] fasta

Positionals:
  fasta TEXT REQUIRED         Patterns FASTA file

Options:
  -h,--help                   Print this help message and exit
  -i,--index TEXT REQUIRED    KeBaB index file
  -o,--output TEXT REQUIRED   Output FASTA file
  -l,--mem-length UINT:POSITIVE [25] 
                              Minimum MEM length (must be greater than k-mer size of index)
  --top-t UINT:POSITIVE       Keep only top-t longest fragments
  -s,--sort                   Sort fragments by length
  -r,--remove-overlaps        Merge overlapping fragments
  -t,--threads UINT:POSITIVE [8] 
                              Number of threads to use
  --no-prefetch               Don't prefetch k-mers to avoid latency
```
To ensure fragments support early stopping (e.g., top t-MEMs), use -s and **do not use** -r.
## Example Usage
### Using KeBaB
```
./kebab build -k 20 -e 0.1 -f 1 -o ~/data/ref_index ~/data/ref.fa
./kebab scan -o ~/data/reads.frag.fa -i ~/data/ref_index.kbb -l 40 ~/data/reads.fa
```
In this example, all MEMs of length at least 40 are contained in the fragments found by breaking reads using KeBaB.
### MEM Finding
The fragments output by a KeBaB scan can be used with MEM-finding tools such as [ropebwt3](https://github.com/lh3/ropebwt3):
```
./ropebwt3 mem -l 40 ~/data/ropebwt3_index.fmd ~/data/reads.frag.fa > ~/data/reads.frag.mems
```
To verify correctness, ``ropefix`` fixes output to match that of running ropebwt3 alone, removing any fragment based notation:
```
./ropefix ~/data/reads.frag.mems > ~/data/reads.mems
```

## Thirdparty

KeBaB utilizes the following third-party libraries:

* [kseq.h](https://lh3lh3.users.sourceforge.net/kseq.shtml) - FASTA parser
* [CLI11](https://github.com/CLIUtils/CLI11) - Command line parser
* [hll](https://github.com/mindis/hll) - HyperLogLog implementation

The k-mer hash implementation is based on the original codebase/paper of [ntHash](https://github.com/bcgsc/ntHash).

## Academic Use
If using this tool or its ideas in an academic setting, please cite:
>Brown, N. K., Depuydt, L., Zakeri, M., Alhadi, A., Allam, N., Begleiter, D., Karpagavalli, N., Khajjayam, S., Wahed, H., Gagie, T., Langmead, B. (2025). KeBaB: $k$-mer based breaking for finding long MEMs. arXiv. arXiv:2502.20338

