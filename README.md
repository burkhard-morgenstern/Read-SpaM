# Read-SpaM


Usage: ./fswm [options] \<filelist\>

format:

\<filelist\>: A plain text file, specifying the relative path to each input dataset.

To create the filelist simply run:

ls ./path/to/input/* > filelist

(assuming all your proteome files (*.fasta, *.faa, etc.) are stored in the input-folder.)

Sequence must be in FASTA format. Each genome must be in its own FASTA file.
There can be multiple header in each FASTA file

>Read1
ATAGTAGATGAT..
>Read2
ATAGTAGATGAT..
>Read3
ATGATGATGATGATG..
..
         
Options:
-h: print this help and exit
-k \<integer\>: pattern weight (default 12)
-l \<integer\>: don't care positions (default 100)
-t \<integer\>: numer of threads (default: 10)
-s \<integer\>: the minimum score of a spaced-word match to be considered homologous (default: 0)


Publication: https://www.biorxiv.org/content/10.1101/550632v1
