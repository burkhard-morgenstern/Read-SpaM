# Read-SpaM

### Compilation
cd into the root directory (containing the 'Makefile') and type:

```	make ```

### Run

```	./readspam [options] <filelist> ```

### Filelist

The program takes a plain text file containing the relative paths to each input
dataset. To create your 'filelist' simply type:

``` ls -1 path/to/input/* > filelist ```

This will list each file in specified directory, one file per line.

### Options
```
        -h: print this help and exit
        -k \<integer\>: pattern weight (default 12)
        -l \<integer\>: don't care positions (default 100)
        -t \<integer\>: numer of threads (default: 10)
        -s \<integer\>: the minimum score of a spaced-word match to be considered homologous (default: 0)
```

### Sequence format:

Sequence must be in FASTA format. All protein sequences of one proteome must be contained in one FASTA file.

Example:
```
    >Read1
    ATAGTAGATGAT..
    >Read2
    ATAGTAGATGAT..
    >Read3
    ATGATGATGATGATG..
    ..
```
### Citation:
```
Scientific publications using Read-SpaM should cite:

A.K. Lau, C.-A. Leimeister, B. Morgenstern (2019)
Read-SpaM: assembly-free and alignment-free comparison of bacterial genomes with low sequencing coverage
bioRxiv, doi:10.1101/550632

https://www.biorxiv.org/content/10.1101/550632v1

```

### Paper Abstract:
```
In many fields of biomedical research, it is important to estimate phylogenetic distances between taxa based on low-coverage sequencing reads. Major applications are, for example, phylogeny reconstruction, species identification from small sequencing samples, or bacterial strain typing in medical diagnostics. Herein, we adapt our previously developed software program Filtered Spaced-Word Matches (FSWM) for alignment-free phylogeny reconstruction to work on unassembled reads; we call this implementation Read-SpaM. Test runs on simulated reads from bacterial genomes show that our approach can estimate phylogenetic distances with high accuracy, even for large evolutionary distances and for very low sequencing coverage.
```

### Contact:
bmorgen@gwdg.de