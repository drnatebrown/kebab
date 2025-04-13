# KeBaB
_**K**-m**e**r **Ba**sed **B**reaking for finding maximal exact matches._

Breaks patterns into fragments, removing sequence content which cannot overlap any maximal exact match (MEM) of some minimum length.

## How-to
### Compile
Creates `./kebab` executable:
```
git clone https://github.com/drnatebrown/KeBaB.git
cd KeBaB
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
  -m,--expected-kmers UINT:NONNEGATIVE [0] 
                              Expected number of k-mers (otherwise estimated)
  -e,--fp-rate FLOAT:FLOAT in [0 - 1] [0.1] 
                              Desired false positive rate (between 0 and 1)
  -f,--hash-funcs UINT:POSITIVE
                              Number of hash functions (otherwise set to minimize index size)
  --kmer-mode ENUM:value in {both->0,canonical->1,forward->2} OR {0,1,2} [1] 
                              K-mer strands to include in the index
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
  -l,--mem-length UINT:POSITIVE [20] 
                              Minimum MEM length (must be >= k-mer size of index)
  -s,--sort                   Sort fragments by length
  -r,--remove-overlaps        Merge overlapping fragments
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
To verify correctness, ``ropefix.sh`` fixes output to match that of running ropebwt3 alone, removing any fragment based notation:
```
./ropefix.sh ~/data/reads.frag.mems > ~/data/reads.mems
