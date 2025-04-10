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
  -o,--output TEXT REQUIRED   Output prefix for .kbb index file
  -m,--expected-kmers UINT:NONNEGATIVE [0] 
                              Expected number of k-mers (if not provided, will be estimated)
  -k,--kmer-size UINT:POSITIVE [20] 
                              k-mer size
  -e,--fp-rate FLOAT:FLOAT in [0 - 1] [0.1] 
                              False positive rate (between 0 and 1)
  -f,--hash-funcs UINT:POSITIVE
                              Number of hash functions
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
  -s,--sort [0]               Sort fragments
  -r,--remove-overlaps [0]    Merge overlapping fragments
  -l,--mem-length UINT:POSITIVE [20] 
                              Minimum MEM length
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
The ``ropefix.sh`` script removes the fragment notation ``[SEQ]:[START]-[END]`` from ropebwt3's output to use global coordinates.
```
./ropefix.sh ~/data/pos_reads.frag.mems > ~/data/pos_reads.mems
