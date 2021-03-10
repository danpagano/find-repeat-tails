#!/bin/bash
####################################################
# Name		: find_nontemplated_tails.sh
# Version   : 0.0.1
# Authors	: Dan Pagano, Yuhan Fei, Scott Kennedy 
# Copyright : Dan Pagano, Yuhan Fei, Scott Kennedy 
# License   : GNU General Public License
####################################################

####################################################
# CONFIGURATION

# specify paths to modules, programs, and files to be used in this script
SCRIPT_DIR=$(readlink -f $0)
SCRIPT_DIR=${SCRIPT_DIR%/*}
export SCRIPT_DIR

# NOTE: this is O2-specific!
# NOTE: loading of particluar modules can also load dependencies. For instance, loading trimgalore also loads python, fastqc, and cutadapt  
module load sratoolkit/2.9.0
module load trimgalore/0.4.5
module load bowtie2/2.3.4.3
module load tophat/2.1.1
module load samtools/1.9
module load star/2.7.3a
PATH_TO_SEQTK=/n/groups/kennedy/pagano/modules/seqtk
PATH_TO_SUBREAD=/n/groups/kennedy/pagano/modules/subread-2.0.0/bin


# check that modules and paths were loaded/set properly
if ! [ $(which fastq-dump) = /n/app/sratoolkit/2.9.0/bin/fastq-dump ]; then echo "ERROR: sratoolkit was not properly loaded"; exit 1; fi
if ! [ $(which trim_galore) = /n/app/TrimGalore/0.4.5/trim_galore ]; then echo "ERROR: trimgalore was not properly loaded"; exit 1; fi
if ! [ $(which bowtie2) = /n/app/bowtie2/2.3.4.3/bin/bowtie2 ]; then echo "ERROR: bowtie2 was not properly loaded"; exit 1; fi
if ! [ $(which tophat) = /n/app/tophat/2.1.1/tophat ]; then echo "ERROR: tophat was not properly loaded"; exit 1; fi
if ! [ $(which samtools) = /n/app/samtools/1.9/bin/samtools ]; then echo "ERROR: samtools was not properly loaded"; exit 1; fi
if ! [ $(which STAR) = /n/app/star/2.7.3a/bin/Linux_x86_64/STAR ]; then echo "ERROR: STAR was not properly loaded"; exit 1; fi
if [ -z "$PATH_TO_SEQTK" ]; then echo "ERROR: path to seqtk was not properly set"; exit 1; fi
if [ -z "$PATH_TO_SUBREAD" ]; then echo "ERROR: path to subread was not properly set"; exit 1; fi

# export paths
export PATH_TO_SEQTK
export PATH_TO_SUBREAD

# default command line configurations
SAMPLE=
REFERENCE_GENOME=
BUILD_GENOME_INDEX=true
REPEAT=
MIN_NUM_REPEATS=10
MIN_NONREPETITIVE_LENGTH=20
MAX_CORES=4
MAX_MEM=4
FEATURE_COUNTS=false
ANNOTATION=
NO_NOVEL_JUNCS=false
COMMAND=

# END CONFIGURATION
####################################################

# check getopt mode
getopt -T
# the above command should give an exit status equal to 4 if the proper version of getopt is being used.
if [ $? -ne 4 ]; then echo "ERROR: Requires enhanced getopt, obtain new version."; exit 1; fi
# the above command outputs the exit status of the last command. -ne means not equal to, so, it reads if 4 doesn't equal 4, then echo and exit.

shopt -s -o nounset

# usage
USAGE="###########################
  FIND NONTEMPLATED TAILS
###########################

DESCRIPTION
FIND NONTEMPLATED TAILS is a workflow designed to identify RNAs that have nontemplated repetitive
sequences appended to their 3' termini.

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
" 

	SCRIPT="find_nontemplated_tails.sh"
	OPTSTRING="s:g:b:r:n:l:c:m:f:a:Vh"
	LOPTSTRING="sample-name:,reference-genome:,build-genome-index:,repeat:,min-number-repeats:,min-nonrep-length:,max-cores:,max-memory:,featureCounts:,annotation:,no-novel-juncs:,verbose,help"

	RESULT=$(getopt -n "$SCRIPT" -o "$OPTSTRING" -l "$LOPTSTRING" -- "$@")
	if [ $? -ne 0 ]
	then
		# parsing error, show usage
		echo "$USAGE" 
		exit 1
	fi

	eval set -- "$RESULT"
	while [ true ] ; do
		case "$1" in
			-s|--sample-name) 
				shift 
				SAMPLE="$1"
			;;
			-g|--reference-genome)
				shift
				REFERENCE_GENOME="$1"
			;;
			-b|--build-genome-index)
				shift
				BUILD_GENOME_INDEX="$1"
			;;
			-r|--repeat)
				shift
				REPEAT="$1"
			;;
			-n|--min-number-repeats)
				shift
				MIN_NUM_REPEATS="$1"
			;;	
			-l|--min-nonrep-length)
				shift
				MIN_NONREPETITIVE_LENGTH="$1"
			;;
			-c|--max-cores)
				shift
				MAX_CORES="$1"
			;;
			-m|--max-memory)
				shift
				MAX_MEM="$1"
			;;				
			-f|--featureCounts)
				shift
				FEATURE_COUNTS="$1"
			;;
			-a|--annotation)
				shift
				ANNOTATION="$1"
			;;
			--no-novel-juncs)
				shift
				NO_NOVEL_JUNCS="$1"
			;;
			-h|--help)
				echo "$USAGE"
				exit 0
			;;
			--)
				shift
				break
			;;
		esac
		shift
	done

	# check that --sample-name argument has been set
	if [ -z "$SAMPLE" ]
	then
		echo -e "$USAGE \n"
		echo -e "ERROR: Need to specify sample name" 
		exit 1
	fi

	# check that --sample-name argument doesn't end in fastq or gz
	if [[ $(echo "$SAMPLE" | sed 's/.*\(..\)/\1/') == "gz" ]] ||
	   [[ $(echo "$SAMPLE" | sed 's/.*\(.....\)/\1/') == "fastq" ]]
	then
		echo -e "$USAGE \n"
		echo -e "ERROR: Sample name should include prefix only and cannot contain fastq or gz extension" 
		exit 1
	fi		

	# set SAMPLE_PREFIX variable
	SAMPLE_PREFIX=$(echo "$SAMPLE" | awk -F '/' '{print $NF}')

	# check that --reference-genome argument has been set
	if [ -z "$REFERENCE_GENOME" ]
	then
		echo -e "$USAGE \n"
		echo -e "ERROR: Need to specify a reference genome" 
		exit 1
	fi

	# set REF_SEQ_PREFIX variable
	REF_SEQ_PREFIX=$(echo "$REFERENCE_GENOME" | awk -F '/' '{print $NF}' | sed 's/.[^.]*$//')

	# check that --build-genome-index is set to true or false
	if [ $BUILD_GENOME_INDEX != "true" ] && 
	   [ $BUILD_GENOME_INDEX != "false" ]
	then
		echo -e "$USAGE \n"
		echo -e "ERROR: --build-genome-index was given the value "$BUILD_GENOME_INDEX". Possible values: {true, false}"
		exit 1
	fi

	# check that --repeat argument has been set
	if [ -z "$REPEAT" ]
	then
		echo -e "$USAGE \n"
		echo -e "ERROR: Need to specify a reference genome" 
		exit 1
	fi

	# covert --repeat argument to uppercase
	REPEAT=$(printf '%s\n' "$REPEAT" | awk '{ print toupper($0) }')

	# check that --repeat argument contains only standard IUPAC nucleotides
	REPEAT_CHARACTERS=$(echo $REPEAT | sed -e "s/./\0\n/g" | sort -u)

	for i in $(echo $REPEAT_CHARACTERS)
	do
		if [[ "$i" != A ]] &&
		   [[ "$i" != T ]] &&	
		   [[ "$i" != C ]] &&
		   [[ "$i" != G ]] &&
		   [[ "$i" != U ]]
		then
	    	echo -e "$USAGE \n"
	    	echo -e "ERROR: The repeat can only be comprised of combinations of A, T, C, G, and U"
	    	exit 1
		fi
	done

	# check that --min-number-repeats argument has been set to an integer
	if ! [[ "$MIN_NUM_REPEATS" =~ ^[0-9]+$ ]]
	then 
    	echo -e "$USAGE \n"
    	echo -e "ERROR: --min-number-repeats must be an integer"  
		exit 1
	fi

	# check that --min-nonrep-length argument has been set to an integer
	if ! [[ "$MIN_NONREPETITIVE_LENGTH" =~ ^[0-9]+$ ]]
	then 
    	echo -e "$USAGE \n"
    	echo -e "ERROR: --min-nonrep-length must be an integer"  
		exit 1
	fi

	# check that --max-cores argument has been set to an integer
	if ! [[ "$MAX_CORES" =~ ^[0-9]+$ ]]
	then 
    	echo -e "$USAGE \n"
    	echo -e "ERROR: --max-cores must be an integer"  
		exit 1
	fi

	# check that --max-cores argument is greater than 3
	if [[ "$MAX_CORES" < 4 ]]
	then 
    	echo -e "$USAGE \n"
    	echo -e "ERROR: --max-cores must be set to 4 or greater"  
		exit 1
	fi

	# check that --max-memory argument has been set to an integer
	if ! [[ "$MAX_MEM" =~ ^[0-9]+$ ]]
	then 
    	echo -e "$USAGE \n"
    	echo -e "ERROR: --max-memory must be an integer"  
		exit 1
	fi

	# check that --featureCounts is set to true or false
	if [ $FEATURE_COUNTS != "true" ] && 
	   [ $FEATURE_COUNTS != "false" ]
	then
		echo -e "$USAGE \n"
		echo -e "ERROR: --featureCounts was given the value "$FEATURE_COUNTS". Possible values: {true, false}"
		exit 1
	fi

	# check that --annotation argument has been set if --featureCounts is set to true
	if [ $FEATURE_COUNTS = "true" ] && 
	   [ -z "$ANNOTATION" ]
	then
		echo -e "$USAGE \n"
		echo -e "ERROR: Need to specify an annotation file when --featureCounts is set to \"true\""
		exit 1
	fi

	# check that --annotation argument, if set, is ends with gtf or gff extension
	if [ ! -z "$ANNOTATION" ] &&
	   [[ $(echo "$ANNOTATION" | sed 's/.*\(...\)/\1/') != "gtf" ]] &&
	   [[ $(echo "$ANNOTATION" | sed 's/.*\(...\)/\1/') != "gff" ]]
	then
		echo -e "$USAGE \n"
		echo -e "ERROR: Annotation file must have a .gtf or .gff extension" 
		exit 1
	fi		

	# check that --no-novel-juncs is set to true or false
	if [ $NO_NOVEL_JUNCS != "true" ] && 
	   [ $NO_NOVEL_JUNCS != "false" ]
	then
		echo -e "$USAGE \n"
		echo -e "ERROR: --no-novel-juncs was given the value "$NO_NOVEL_JUNCS". Possible values: {true, false}"
		exit 1
	fi

	# set command to appropriate workflow given length of 1 repeat unit
	if [ ${#REPEAT} == 1 ]; then COMMAND=""$SCRIPT_DIR"/find_pNs.sh"; fi
	if [ ${#REPEAT} == 2 ]; then COMMAND=""$SCRIPT_DIR"/find_pNNs.sh"; fi
	if [ ${#REPEAT} == 3 ]; then COMMAND=""$SCRIPT_DIR"/find_pNNNs.sh"; fi
	if [ ${#REPEAT} == 4 ]; then COMMAND=""$SCRIPT_DIR"/find_pNNNNs.sh"; fi
	if [ ${#REPEAT} == 5 ]; then COMMAND=""$SCRIPT_DIR"/find_pNNNNNs.sh"; fi

	if [ -z "$COMMAND" ]
	then
		echo -e "$USAGE \n"
		echo -e "ERROR: Repeat unit length exceeds 5nt. At the moment, this script can only accommodate repeats of 5nt or less." 
		exit 1
	fi 	

	export SAMPLE
	export REFERENCE_GENOME
	export BUILD_GENOME_INDEX
	export REPEAT
	export MIN_NUM_REPEATS
	export MIN_NONREPETITIVE_LENGTH
	export MAX_CORES
	export MAX_MEM
	export FEATURE_COUNTS
	export ANNOTATION
	export NO_NOVEL_JUNCS
	export SAMPLE_PREFIX
	export REF_SEQ_PREFIX

	#run command
	sbatch "$COMMAND"

	STATUS=$?

	exit $STATUS