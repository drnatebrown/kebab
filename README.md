# KeBaB
**K**-m**e**r **Ba**sed **B**reaking for finding super-maximal exact matches.

## Compile
Creates `kebab` executable:
```
git clone https://github.com/drnatebrown/KeBaB.git
cd KeBaB
make
```
## Build
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
```
## Scan
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
```
## Example
```
./kebab build -k 20 -e 0.1 -f 1 -o ~/data/ref_index ~/data/ref.fa
./kebab scan -o ~/data/pos_reads.frag.fa -i ~/data/ref_index.kbb -l 40 ~/data/pos_reads.fa
```
Use `-s` to ensure MEMs are sorted for early stopping approaches (top t-MEMs).
## Usage
The fragments output by a KeBaB scan can be used with MEM-finding tools such as [ropebwt3](https://github.com/lh3/ropebwt3):
```
./ropebwt3 mem -l40 ~/data/ropebwt3_index.fmd ~/data/pos_reads.frag.fa > ~/data/pos_reads.frag.mems
```
The ``ropefix.sh`` script fixes the output to match that of running ropebwt3 alone by rectifying differences due to the fragment notation:
```
./ropefix.sh ~/data/pos_reads.frag.mems > ~/data/pos_reads.mems
