# Find RNAs with nontemplated repetitive tails

This is a workflow designed to identify RNAs that have nontemplated repetitive
sequences appended to their 3' termini.
```
USAGE
find_nontemplated_tails.sh -s <sample-name> -g <reference-genome> -r <repeat> -n
<min_number_repeats> -l <min-nonrep-length> [options]

EXAMPLE
find_nontemplated_tails.sh -s SRR11059947 -g
/n/groups/kennedy/pagano/pUG_RNAs/references/SARS-CoV-2/genome_for_tophat/NC_045512.2.fa -r GT -n 10
-l 20 [options]

REQUIRED ARGUMENTS
-s, --sample-name <string>          The name or prefix given to the fastq input file(s) (not including
                                    the extensions). Note: Input files must be gz compressed and
                                    have the extension fastq.gz if read-type is single-end or
                                    _1.fastq.gz and _2.fastq.gz if read-type is paired-end. If data
                                    is stored locally, then set --sample-name to full path and file
                                    prefix (e.g. PATH/TO/LOCAL/DATA/YY1000). If data is to be
                                    downloaded from NCBI SRA, then set --sample-name to SRR#, ERR#,
                                    or DRR# (e.g. SRR11059947).

-g, --reference-genome <string>     Reference genome with .fa extension. For larger genomes (e.g.
                                    human), a bowtie2 genome index must be in the same directory as
                                    the fasta file.

-b, --build-genome-index <string>   This should be set to \"true\" for small genomes (e.g. viral
                                    genomes) and \"false\" for larger genomes (e.g. human genomes).
                                    Default: true. Possible values: {true, false}

-r, --repeat <string>               Repeat to find.

-n, --min-number-repeats <integer>  Minimum number of repeats a sequence must have to be considered
                                    a candidate.

-l, --min-nonrep-length <integer>   Minimum length of nonrepetitive portion of read for read to be
                                    considered a candidate.

-c, --max-cores <integer>           Maximum number of cores to allocate. This workflow requires, at a
                                    minimum, 4 cores. Default: 4.

-m, --max-memory <integer>          Maximum amount of memory (in Gb) to allocate. Default: 4.


ADVANCED ARGUMENTS
-f, --featureCounts <string>        Count the number of candidates mapping to each gene. Requires a
                                    gtf/gff file be specified with the --annotation argument.
                                    Default: false. Possible values: {true, false}

-a, --annotation <string>           Annotation file.

--no-novel-juncs <string>           By default, TopHat will look for novel splice junctions in
                                    reads. If FIND NONTEMPLATED TAILS fails because of a TopHat
                                    error \"Error: Splice sequence indexing failed with err =1\",
                                    then set --no-novel-juncs to true. Otherwise --no-novel-juncs
                                    should be set to false. Default: false. Possible values: {true,
                                    false}

-h, --help                          Gives this help message.
```
