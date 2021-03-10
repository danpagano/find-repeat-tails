#!/bin/bash
#SBATCH -p short
#SBATCH -t 0-04:00
#SBATCH -c 4
#SBATCH --mem=4G
#SBATCH -o find_pNNs_%j.out
#SBATCH -e find_pNNs_%j.err 
#SBATCH--mail-type=FAIL

# Print command line inputs and configuration to log file
	printf -- "%s\n" "FIND NONTEMPLATED "$REPEAT" TAILS LOG FILE" > "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
	printf -- "%s\n" "[`date`] " "" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt

	printf -- "%s\n" "-----COMMAND-LINE-INPUTS-----" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
	printf -- "%s\n" "Sample: "$SAMPLE_PREFIX"" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
	printf -- "%s\n" "Reference genome: "$REFERENCE_GENOME"" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
	printf -- "%s\n" "Build reference genome index: "$BUILD_GENOME_INDEX"" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
	printf -- "%s\n" "Repeat to find: "$REPEAT"" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
	printf -- "%s\n" "Minimum number of repeats: "$MIN_NUM_REPEATS"" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
	printf -- "%s\n" "Minimum nonrepetitive length of read: "$MIN_NONREPETITIVE_LENGTH"" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
	printf -- "%s\n" "Run featureCounts: "$FEATURE_COUNTS"" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
	if [ ! -z "$ANNOTATION" ]; then printf -- "%s\n" "Annotation file for featureCounts: "$ANNOTATION"" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt; fi
	printf -- "%s\n" "Run TopHat in --no-novel-juncs mode: "$NO_NOVEL_JUNCS"" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
	printf -- "%s\n" "Max cores: "$MAX_CORES"" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
	printf -- "%s\n" "Max memory (Gb): "$MAX_MEM"" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt

	printf -- "%s\n" " " "-------CONFIGURATIONS--------" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
	printf -- "%s\n" "Python: "$(which python)"" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
	printf -- "%s\n" "SRA Toolkit: "$(which fastq-dump)"" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
	printf -- "%s\n" "Trim Galore!: "$(which trim_galore)"" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
	printf -- "%s\n" "Cutadapt: "$(which cutadapt)"" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
	printf -- "%s\n" "Bowtie 2: "$(which bowtie2)"" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
	printf -- "%s\n" "TopHat: "$(which tophat)"" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
	printf -- "%s\n" "Samtools: "$(which samtools)"" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
	printf -- "%s\n" "STAR: "$(which STAR)"" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
	printf -- "%s\n" "Path to seqtk: "$PATH_TO_SEQTK"" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt 
	printf -- "%s\n" "Path to subread: "$PATH_TO_SUBREAD"" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt 
	printf -- "%s\n" " " >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt 

# Make analysis directory if it doesn't already exist 
	if [ -d ./"$SAMPLE_PREFIX"_find_nontemplated_tails ]
	then 
		echo "[`date`] Found "$SAMPLE_PREFIX"_find_nontemplated_tails directory" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
	else
		echo "[`date`] Making "$SAMPLE_PREFIX"_find_nontemplated_tails directory" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
		mkdir ./"$SAMPLE_PREFIX"_find_nontemplated_tails
	fi

# Change to analysis directory
	mv "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt ./"$SAMPLE_PREFIX"_find_nontemplated_tails
	cd ./"$SAMPLE_PREFIX"_find_nontemplated_tails

# Find fastq input files, or download if neccessary, and determine read-type
	if [ -f "$SAMPLE".fastq.gz ]
	then
		echo "[`date`] Found "$SAMPLE_PREFIX" compressed fastq input file" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
		echo "[`date`] Read-type was determined to be single-end. Running analysis in single-end mode" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
		READ_TYPE=single-end
	else
		if [ -f "$SAMPLE"_1.fastq.gz ] &&
		   [ -f "$SAMPLE"_2.fastq.gz ]
		then
			echo "[`date`] Found "$SAMPLE_PREFIX" compressed fastq input files" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
			echo "[`date`] Read-type was determined to be paired-end. Running analysis in paired-end mode" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
			READ_TYPE=paired-end
		else
			if [[ $(echo "$SAMPLE" | cut -c 1-3) == "SRR" ]] ||
	           [[ $(echo "$SAMPLE" | cut -c 1-3) == "ERR" ]] ||
	           [[ $(echo "$SAMPLE" | cut -c 1-3) == "DRR" ]]
	        then
	        	echo "[`date`] Downloading "$SAMPLE_PREFIX" from NCBI SRA" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
	        	fastq-dump --outdir ./ --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 --clip "$SAMPLE" 
	        	echo "[`date`] Download complete" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt

				if [ -f "$SAMPLE"_pass_1.fastq.gz ] &&
				   [ -f "$SAMPLE"_pass_2.fastq.gz ]	
				then
					READ_TYPE=paired-end
					mv "$SAMPLE"_pass_1.fastq.gz "$SAMPLE"_1.fastq.gz
					mv "$SAMPLE"_pass_2.fastq.gz "$SAMPLE"_2.fastq.gz
					echo "[`date`] Read-type was determined to be paired-end. Running analysis in paired-end mode" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
				else
					READ_TYPE=single-end
					mv "$SAMPLE"_pass.fastq.gz "$SAMPLE".fastq.gz
					echo "[`date`] Read-type was determined to be single-end. Running analysis in single-end mode" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
				fi
			else
				echo "[`date`] Unable to find "$SAMPLE_PREFIX" compressed fastq input files locally or on NCBI SRA. Check that the input files are properly named and that the --sample-name argument has been properly set" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt	
				mv "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt "$SAMPLE_PREFIX"_vda_log_$(date +"%Y-%m-%d_%H:%M:%S").txt
				echo "[`date`] Unable to find "$SAMPLE_PREFIX" compressed fastq input files locally or on NCBI SRA. Check that the input files are properly named and that the --sample-name argument has been properly set"
				exit 1
			fi
		fi
	fi

	if [ $READ_TYPE != "single-end" ] && [ $READ_TYPE != "paired-end" ]
	then
		echo "[`date`] ERROR: read-type is not set correctly. In theory, it shouldn't be possible to get this error message. If you do, then something is inherently wrong with this script." >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
		mv "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt "$SAMPLE_PREFIX"_vda_log_$(date +"%Y-%m-%d_%H:%M:%S").txt
		echo "[`date`] ERROR: read-type is not set correctly. In theory, it shouldn't be possible to get this error message. If you do, then something is inherently wrong with this script."
		exit 1
	fi

# Make repeat tails and heads with command line input variables
	# convert U to T
	REPEAT=$(echo "$REPEAT" | sed 's/U/T/g')

	# make REPEAT_REVCOMP variable
	paste -d '\n' <(echo \>REPEAT) <(echo "$REPEAT") > repeat.fa
	REPEAT_REVCOMP=$("$PATH_TO_SEQTK"/seqtk seq -r repeat.fa | grep -v ">")
	rm repeat.fa

	# make variables for different frames of the tail
	REPEAT_FRAME_1="$REPEAT"
	REPEAT_FRAME_2=$(echo "$REPEAT""$REPEAT" | cut -c 2-3)

	# make tail variable for each repeat frame
	for i in $(seq 1 "$MIN_NUM_REPEATS")
	do
		echo "$REPEAT_FRAME_1" | tr -d '\n' >> "$MIN_NUM_REPEATS"_"$REPEAT_FRAME_1"
	done
	TAIL_1=$(cat "$MIN_NUM_REPEATS"_"$REPEAT_FRAME_1")

	for i in $(seq 1 "$MIN_NUM_REPEATS")
	do
		echo "$REPEAT_FRAME_2" | tr -d '\n' >> "$MIN_NUM_REPEATS"_"$REPEAT_FRAME_2"
	done
	TAIL_2=$(cat "$MIN_NUM_REPEATS"_"$REPEAT_FRAME_2")

	# make head variable for each repeat frame
	echo \>"$REPEAT_FRAME_1" | cat - "$MIN_NUM_REPEATS"_"$REPEAT_FRAME_1" > "$MIN_NUM_REPEATS"_"$REPEAT_FRAME_1".fa
	HEAD_1=$("$PATH_TO_SEQTK"/seqtk seq -r "$MIN_NUM_REPEATS"_"$REPEAT_FRAME_1".fa | grep -v ">")

	echo \>"$REPEAT_FRAME_2" | cat - "$MIN_NUM_REPEATS"_"$REPEAT_FRAME_2" > "$MIN_NUM_REPEATS"_"$REPEAT_FRAME_2".fa
	HEAD_2=$("$PATH_TO_SEQTK"/seqtk seq -r "$MIN_NUM_REPEATS"_"$REPEAT_FRAME_2".fa | grep -v ">")

	printf -- "%s\n" " " >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt 
	printf -- "%s\n" "Repeat: "$REPEAT"" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
	printf -- "%s\n" "Repeat reverse complement: "$REPEAT_REVCOMP"" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
	printf -- "%s\n" "Repeat tail frame 1: "$TAIL_1"" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
	printf -- "%s\n" "Repeat tail frame 2: "$TAIL_2"" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
	printf -- "%s\n" "Repeat head frame 1: "$HEAD_1"" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
	printf -- "%s\n" "Repeat head frame 2: "$HEAD_2"" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt 
	printf -- "%s\n" " " >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt 

	rm "$MIN_NUM_REPEATS"_"$REPEAT_FRAME_1"
	rm "$MIN_NUM_REPEATS"_"$REPEAT_FRAME_2"
	rm "$MIN_NUM_REPEATS"_"$REPEAT_FRAME_1".fa
	rm "$MIN_NUM_REPEATS"_"$REPEAT_FRAME_2".fa

	# make variables to be used with cutadapt for repeat removal
    let MIN_TAIL_LENGTH="$MIN_NUM_REPEATS"*${#REPEAT}
	i=0
	n=150
	let m=${n}-1
	let MAX_NUM_REPEATS=${n}/${#REPEAT}
	let REMOVE_UP_TO=${n}-"$MIN_TAIL_LENGTH"
	# set a variable to largest tail size for each repeat frame
	for i in $(seq 1 "$MAX_NUM_REPEATS")
	do
		echo "$REPEAT_FRAME_1" | tr -d '\n' >> "$MAX_NUM_REPEATS"_"$REPEAT_FRAME_1"
	done
	MAX_TAIL_1=$(cat "$MAX_NUM_REPEATS"_"$REPEAT_FRAME_1")

	for i in $(seq 1 "$MAX_NUM_REPEATS")
	do
		echo "$REPEAT_FRAME_2" | tr -d '\n' >> "$MAX_NUM_REPEATS"_"$REPEAT_FRAME_2"
	done
	MAX_TAIL_2=$(cat "$MAX_NUM_REPEATS"_"$REPEAT_FRAME_2")

	# set a variable to largest head size for each repeat frame
	echo \>"$REPEAT_FRAME_1" | cat - "$MAX_NUM_REPEATS"_"$REPEAT_FRAME_1" > "$MAX_NUM_REPEATS"_"$REPEAT_FRAME_1".fa
	MAX_HEAD_1=$("$PATH_TO_SEQTK"/seqtk seq -r "$MAX_NUM_REPEATS"_"$REPEAT_FRAME_1".fa | grep -v ">")

	echo \>"$REPEAT_FRAME_2" | cat - "$MAX_NUM_REPEATS"_"$REPEAT_FRAME_2" > "$MAX_NUM_REPEATS"_"$REPEAT_FRAME_2".fa
	MAX_HEAD_2=$("$PATH_TO_SEQTK"/seqtk seq -r "$MAX_NUM_REPEATS"_"$REPEAT_FRAME_2".fa | grep -v ">")

	printf -- "%s\n" "Max tail frame 1: "$MAX_TAIL_1"" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
	printf -- "%s\n" "Max tail frame 2: "$MAX_TAIL_2"" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
	printf -- "%s\n" "Max head frame 1: "$MAX_HEAD_1"" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
	printf -- "%s\n" "Max head frame 2: "$MAX_HEAD_2"" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt 
	printf -- "%s\n" " " >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt 

	rm "$MAX_NUM_REPEATS"_"$REPEAT_FRAME_1"
	rm "$MAX_NUM_REPEATS"_"$REPEAT_FRAME_2"
	rm "$MAX_NUM_REPEATS"_"$REPEAT_FRAME_1".fa
	rm "$MAX_NUM_REPEATS"_"$REPEAT_FRAME_2".fa

# Build genome index
	if [ $BUILD_GENOME_INDEX = "true" ]
	then
		if [ -f references/"$REF_SEQ_PREFIX"_"$REPEAT"_tail/"$REF_SEQ_PREFIX"_dummied.1.bt2 ] &&
		   [ -f references/"$REF_SEQ_PREFIX"_"$REPEAT"_tail/"$REF_SEQ_PREFIX"_dummied.2.bt2 ] &&
		   [ -f references/"$REF_SEQ_PREFIX"_"$REPEAT"_tail/"$REF_SEQ_PREFIX"_dummied.3.bt2 ] &&
		   [ -f references/"$REF_SEQ_PREFIX"_"$REPEAT"_tail/"$REF_SEQ_PREFIX"_dummied.4.bt2 ] &&
		   [ -f references/"$REF_SEQ_PREFIX"_"$REPEAT"_tail/"$REF_SEQ_PREFIX"_dummied.rev.1.bt2 ] &&
		   [ -f references/"$REF_SEQ_PREFIX"_"$REPEAT"_tail/"$REF_SEQ_PREFIX"_dummied.rev.2.bt2 ]
		then 
			echo "[`date`] Genome index for "$REF_SEQ_PREFIX"_"$REPEAT"_tail already exists" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
			REF_FOR_TOPHAT=$(pwd)/references/"$REF_SEQ_PREFIX"_"$REPEAT"_tail/"$REF_SEQ_PREFIX"_dummied
		else 
			if [ ! -d references/"$REF_SEQ_PREFIX"_"$REPEAT"_tail ]
			then 
				echo "[`date`] Making "$REF_SEQ_PREFIX"_"$REPEAT"_tail reference directory" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt 
				mkdir -p references/"$REF_SEQ_PREFIX"_"$REPEAT"_tail
			fi
			cd references/"$REF_SEQ_PREFIX"_"$REPEAT"_tail	
			echo "[`date`] Adding dummy "$REPEAT"-repeat-containing chromosome to "$REF_SEQ_PREFIX" genome" >> ./../../"$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt 
			printf -- "%s\n" ""\>dummy_chr_"$REPEAT" TACATCGAGACGGAGCAGGCCCGTTCCAACATTCGTGAGGACGTTGCAAAGGAGGAAGTTCGTCATATTA"$TAIL_1""" | cat "$REFERENCE_GENOME" - > "$REF_SEQ_PREFIX"_dummied.fa
			echo "[`date`] Building genome index" >> ./../../"$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt 
			bowtie2-build "$REF_SEQ_PREFIX"_dummied.fa "$REF_SEQ_PREFIX"_dummied
			REF_FOR_TOPHAT=$(pwd)/"$REF_SEQ_PREFIX"_dummied
			cd ./../../
		fi
	else
		if [ $BUILD_GENOME_INDEX = "false" ]
		then
			REF_FOR_TOPHAT=$(echo "$REFERENCE_GENOME" | sed 's/.[^.]*$//')
		fi
	fi

###################################################################################################
###################################################################################################

# paired-end read workflow
	if [ $READ_TYPE = "paired-end" ]
	then
		#0.make analysis directories
		if [ ! -d "$SAMPLE_PREFIX"_1 ]; then mkdir -p "$SAMPLE_PREFIX"_1; fi
		if [ ! -d "$SAMPLE_PREFIX"_2 ]; then mkdir -p "$SAMPLE_PREFIX"_2; fi	

		#1.remove adapters and low quality bases
		if [ -f "$SAMPLE_PREFIX"_1/01-TrimReads/"$SAMPLE_PREFIX"_1_trimmed.fq ] &&
		   [ -f "$SAMPLE_PREFIX"_2/01-TrimReads/"$SAMPLE_PREFIX"_2_trimmed.fq ]
		then 
			echo "[`date`] Read trimming already complete" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
		else 
			if [ -d "$SAMPLE_PREFIX"_1/01-TrimReads ] ||
			   [ -d "$SAMPLE_PREFIX"_2/01-TrimReads ]
			then
				echo "[`date`] Read trimming is incomplete... Removing "$SAMPLE_PREFIX"_1/01-TrimReads and "$SAMPLE_PREFIX"_2/01-TrimReads directories and starting over." >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
				rm -rf "$SAMPLE_PREFIX"_1/01-TrimReads
				rm -rf "$SAMPLE_PREFIX"_2/01-TrimReads
			fi
			echo "[`date`] Trimming reads" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
			mkdir -p "$SAMPLE_PREFIX"_1/01-TrimReads
			mkdir -p "$SAMPLE_PREFIX"_2/01-TrimReads
			
			# trimming reads with mostly default parameters
			trim_galore -q 20 --dont_gzip --trim-n -o ./"$SAMPLE_PREFIX"_1/01-TrimReads "$SAMPLE"_1.fastq.gz &
			trim_galore -q 20 --dont_gzip --trim-n -o ./"$SAMPLE_PREFIX"_2/01-TrimReads "$SAMPLE"_2.fastq.gz &
			wait

			if [ -f "$SAMPLE_PREFIX"_1/01-TrimReads/"$SAMPLE_PREFIX"_1_trimmed.fq ] &&
			   [ -f "$SAMPLE_PREFIX"_2/01-TrimReads/"$SAMPLE_PREFIX"_2_trimmed.fq ]
			then 
				echo "[`date`] Read trimming complete" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
			else 
				echo "[`date`] ERROR: Read trimming failed. Check log for details" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
				mv "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt "$SAMPLE_PREFIX"_log_$(date +"%Y-%m-%d_%H:%M:%S").txt
				echo "[`date`] ERROR: Read trimming failed. Check log for details"
				exit 1
			fi
		fi

		# align trimmed reads to genome
		if [ -f "$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/star/"$SAMPLE_PREFIX"_Log.final.out ] &&
		   [ -f "$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/star/"$SAMPLE_PREFIX"_Log.final.out ]
		then 
			echo "[`date`] Trimmed read alignemnt to "$REF_SEQ_PREFIX" genome already complete" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
		else 
			echo "[`date`] Aligning trimmed reads to "$REF_SEQ_PREFIX" genome" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
			if [ ! -d references/"$REF_SEQ_PREFIX"_STAR_genomeDir ]; then mkdir -p references/"$REF_SEQ_PREFIX"_STAR_genomeDir; fi
			if [ ! -d "$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/star ]; then mkdir -p "$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/star; fi
			if [ ! -d "$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/star ]; then mkdir -p "$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/star; fi

			if [ ! -f references/"$REF_SEQ_PREFIX"_STAR_genomeDir/SAindex ] &&
			   [ ! -f references/"$REF_SEQ_PREFIX"_STAR_genomeDir/SA ] &&
			   [ ! -f references/"$REF_SEQ_PREFIX"_STAR_genomeDir/genomeParameters.txt ] &&
			   [ ! -f references/"$REF_SEQ_PREFIX"_STAR_genomeDir/Genome ] &&
			   [ ! -f references/"$REF_SEQ_PREFIX"_STAR_genomeDir/chrStart.txt ] &&
			   [ ! -f references/"$REF_SEQ_PREFIX"_STAR_genomeDir/chrNameLength.txt ] &&
			   [ ! -f references/"$REF_SEQ_PREFIX"_STAR_genomeDir/chrName.txt ] &&
			   [ ! -f references/"$REF_SEQ_PREFIX"_STAR_genomeDir/chrLength.txt ]
			then
				REF_SEQ_LENGTH=$("$PATH_TO_SEQTK"/seqtk comp "$REFERENCE_GENOME" | awk -F "\t" '{print $2}' | paste -sd+ - | bc)
				Nbases=$(echo "(l("$REF_SEQ_LENGTH")/l(2))/2 -1" | bc -l)
				printf -v NbasesRounded %.0f "$Nbases"

				if [[ "$NbasesRounded" -lt 14 ]]
				then 
					genomeSAindexNbases="$NbasesRounded"
				else
					genomeSAindexNbases=14
				fi 

				STAR \
				--runMode genomeGenerate \
				--genomeSAindexNbases "$genomeSAindexNbases" \
				--genomeDir ./references/"$REF_SEQ_PREFIX"_STAR_genomeDir/ \
				--genomeFastaFiles "$REFERENCE_GENOME"

				mv Log.out ./references/"$REF_SEQ_PREFIX"_STAR_genomeDir/
			fi

			STAR \
			--runMode alignReads \
			--runThreadN "$MAX_CORES" \
			--genomeDir ./references/"$REF_SEQ_PREFIX"_STAR_genomeDir/ \
			--outSAMtype BAM Unsorted \
			--outFileNamePrefix ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/star/"$SAMPLE_PREFIX"_ \
			--readFilesIn ./"$SAMPLE_PREFIX"_1/01-TrimReads/"$SAMPLE_PREFIX"_1_trimmed.fq

			STAR \
			--runMode alignReads \
			--runThreadN "$MAX_CORES" \
			--genomeDir ./references/"$REF_SEQ_PREFIX"_STAR_genomeDir/ \
			--outSAMtype BAM Unsorted \
			--outFileNamePrefix ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/star/"$SAMPLE_PREFIX"_ \
			--readFilesIn ./"$SAMPLE_PREFIX"_2/01-TrimReads/"$SAMPLE_PREFIX"_2_trimmed.fq

			if [ -f "$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/star/"$SAMPLE_PREFIX"_Log.final.out ] &&
			   [ -f "$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/star/"$SAMPLE_PREFIX"_Log.final.out ]
			then 
				echo "[`date`] Trimmed read alignemnt to "$REF_SEQ_PREFIX" genome complete" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
			else 
				echo "[`date`] ERROR: Trimmed read alignemnt to "$REF_SEQ_PREFIX" genome failed. Check log for details" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
				mv "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt "$SAMPLE_PREFIX"_log_$(date +"%Y-%m-%d_%H:%M:%S").txt
				echo "[`date`] ERROR: Trimmed read alignemnt to "$REF_SEQ_PREFIX" genome failed. Check log for details"
				exit 1
			fi
		fi

		#2.Select reads with desired repeat and length
		if [ -f "$SAMPLE_PREFIX"_1/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT"_tail.fastq ] &&
		   [ -f "$SAMPLE_PREFIX"_1/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT_REVCOMP"_head.fastq ] &&
	       [ -f "$SAMPLE_PREFIX"_2/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT"_tail.fastq ] &&
	       [ -f "$SAMPLE_PREFIX"_2/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT_REVCOMP"_head.fastq ]
		then 
			echo "[`date`] Read selection already complete" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
		else 
			echo "[`date`] Selecting reads with at least "$MIN_NUM_REPEATS" repeats of "$REPEAT_FRAME_1"/"$REPEAT_FRAME_2" at their 3' end" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
			if [ ! -d "$SAMPLE_PREFIX"_1/02-GrabReads ]; then mkdir -p "$SAMPLE_PREFIX"_1/02-GrabReads; fi
			if [ ! -d "$SAMPLE_PREFIX"_2/02-GrabReads ]; then mkdir -p "$SAMPLE_PREFIX"_2/02-GrabReads; fi

			# adding dummy reads to the trimmed reads from the previous step
			DUMMY_READ_SEQ_TAIL=TACATCGAGACGGAGCAGGCCCGTTCCAACATTCGTGAGGACGTTGCAAAGGAGGAAGTTCGTCATATTA"$TAIL_1"
			DUMMY_READ_SEQ_HEAD="$HEAD_1"TAATATGACGAACTTCCTCCTTTGCAACGTCCTCACGAATGTTGGAACGGGCCTGCTCCGTCTCGATGTA
			DUMMY_READ_LEN=${#DUMMY_READ_SEQ_TAIL}

			for i in $(seq 1 "$DUMMY_READ_LEN")
			do
				echo I | tr -d '\n' >> DUMMY_READ_QUALITY
			done
			DUMMY_READ_QUAL=$(cat DUMMY_READ_QUALITY)
			rm DUMMY_READ_QUALITY

			printf -- "%s\n" ""@dummy"$REPEAT".1 "$DUMMY_READ_SEQ_TAIL" + "$DUMMY_READ_QUAL" @dummy"$REPEAT_REVCOMP".1 "$DUMMY_READ_SEQ_HEAD" + "$DUMMY_READ_QUAL""" | cat ./"$SAMPLE_PREFIX"_1/01-TrimReads/"$SAMPLE_PREFIX"_1_trimmed.fq - > ./"$SAMPLE_PREFIX"_1/02-GrabReads/"$SAMPLE_PREFIX"_1_trimmed_"$REPEAT"_dummied.fq &
			printf -- "%s\n" ""@dummy"$REPEAT".2 "$DUMMY_READ_SEQ_TAIL" + "$DUMMY_READ_QUAL" @dummy"$REPEAT_REVCOMP".2 "$DUMMY_READ_SEQ_HEAD" + "$DUMMY_READ_QUAL""" | cat ./"$SAMPLE_PREFIX"_2/01-TrimReads/"$SAMPLE_PREFIX"_2_trimmed.fq - > ./"$SAMPLE_PREFIX"_2/02-GrabReads/"$SAMPLE_PREFIX"_2_trimmed_"$REPEAT"_dummied.fq &
			wait

			# grabbing reads with repeats, frame 1
			grep -i --no-group-separator ${TAIL_1}$ -B 1 -A 2 ./"$SAMPLE_PREFIX"_1/02-GrabReads/"$SAMPLE_PREFIX"_1_trimmed_"$REPEAT"_dummied.fq > ./"$SAMPLE_PREFIX"_1/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT"_tail.fastq &
			grep -i --no-group-separator ^${HEAD_1} -B 1 -A 2 ./"$SAMPLE_PREFIX"_1/02-GrabReads/"$SAMPLE_PREFIX"_1_trimmed_"$REPEAT"_dummied.fq > ./"$SAMPLE_PREFIX"_1/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT_REVCOMP"_head.fastq &
			grep -i --no-group-separator ${TAIL_1}$ -B 1 -A 2 ./"$SAMPLE_PREFIX"_2/02-GrabReads/"$SAMPLE_PREFIX"_2_trimmed_"$REPEAT"_dummied.fq > ./"$SAMPLE_PREFIX"_2/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT"_tail.fastq &
			grep -i --no-group-separator ^${HEAD_1} -B 1 -A 2 ./"$SAMPLE_PREFIX"_2/02-GrabReads/"$SAMPLE_PREFIX"_2_trimmed_"$REPEAT"_dummied.fq > ./"$SAMPLE_PREFIX"_2/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT_REVCOMP"_head.fastq &
			wait

			# grabbing reads with repeats, frame 2
			grep -i --no-group-separator ${TAIL_2}$ -B 1 -A 2 ./"$SAMPLE_PREFIX"_1/02-GrabReads/"$SAMPLE_PREFIX"_1_trimmed_"$REPEAT"_dummied.fq >> ./"$SAMPLE_PREFIX"_1/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT"_tail.fastq &
			grep -i --no-group-separator ^${HEAD_2} -B 1 -A 2 ./"$SAMPLE_PREFIX"_1/02-GrabReads/"$SAMPLE_PREFIX"_1_trimmed_"$REPEAT"_dummied.fq >> ./"$SAMPLE_PREFIX"_1/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT_REVCOMP"_head.fastq &
			grep -i --no-group-separator ${TAIL_2}$ -B 1 -A 2 ./"$SAMPLE_PREFIX"_2/02-GrabReads/"$SAMPLE_PREFIX"_2_trimmed_"$REPEAT"_dummied.fq >> ./"$SAMPLE_PREFIX"_2/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT"_tail.fastq &
			grep -i --no-group-separator ^${HEAD_2} -B 1 -A 2 ./"$SAMPLE_PREFIX"_2/02-GrabReads/"$SAMPLE_PREFIX"_2_trimmed_"$REPEAT"_dummied.fq >> ./"$SAMPLE_PREFIX"_2/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT_REVCOMP"_head.fastq &
			wait

			if [ -f "$SAMPLE_PREFIX"_1/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT"_tail.fastq ] &&
			   [ -f "$SAMPLE_PREFIX"_1/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT_REVCOMP"_head.fastq ] &&
	           [ -f "$SAMPLE_PREFIX"_2/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT"_tail.fastq ] &&
	           [ -f "$SAMPLE_PREFIX"_2/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT_REVCOMP"_head.fastq ]
			then 
				echo "[`date`] Read selection complete" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
			else 
				echo "[`date`] ERROR: Read selection failed. Check log for details" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
				mv "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt "$SAMPLE_PREFIX"_log_$(date +"%Y-%m-%d_%H:%M:%S").txt
				echo "[`date`] ERROR: Read selection failed. Check log for details"
				exit 1
			fi
		fi

		#3.Remove repeats from reads
		if [ -f "$SAMPLE_PREFIX"_1/03-RemoveRepeats/"$MIN_TAIL_LENGTH"_"$REPEAT"_tail.fastq ] &&
		   [ -f "$SAMPLE_PREFIX"_1/03-RemoveRepeats/"$MIN_TAIL_LENGTH"_"$REPEAT_REVCOMP"_head.fastq ] &&
	       [ -f "$SAMPLE_PREFIX"_2/03-RemoveRepeats/"$MIN_TAIL_LENGTH"_"$REPEAT"_tail.fastq ] &&
	       [ -f "$SAMPLE_PREFIX"_2/03-RemoveRepeats/"$MIN_TAIL_LENGTH"_"$REPEAT_REVCOMP"_head.fastq ]
		then 
			echo "[`date`] Repeat removal already complete" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
		else 
			echo "[`date`] Removing repeats from reads" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
			if [ ! -d "$SAMPLE_PREFIX"_1/03-RemoveRepeats ]; then mkdir -p "$SAMPLE_PREFIX"_1/03-RemoveRepeats; fi
			if [ ! -d "$SAMPLE_PREFIX"_2/03-RemoveRepeats ]; then mkdir -p "$SAMPLE_PREFIX"_2/03-RemoveRepeats; fi

			# remove repeats
			cutadapt -e 0 --minimum-length "$MIN_NONREPETITIVE_LENGTH" -a ${MAX_TAIL_1}$ -a ${MAX_TAIL_2}$ ./"$SAMPLE_PREFIX"_1/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT"_tail.fastq -o ./"$SAMPLE_PREFIX"_1/03-RemoveRepeats/${n}_"$REPEAT"_tail.fastq &
			cutadapt -e 0 --minimum-length "$MIN_NONREPETITIVE_LENGTH" -g ^${MAX_HEAD_1} -g ^${MAX_HEAD_2} ./"$SAMPLE_PREFIX"_1/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT_REVCOMP"_head.fastq -o ./"$SAMPLE_PREFIX"_1/03-RemoveRepeats/${n}_"$REPEAT_REVCOMP"_head.fastq &
			cutadapt -e 0 --minimum-length "$MIN_NONREPETITIVE_LENGTH" -a ${MAX_TAIL_1}$ -a ${MAX_TAIL_2}$ ./"$SAMPLE_PREFIX"_2/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT"_tail.fastq -o ./"$SAMPLE_PREFIX"_2/03-RemoveRepeats/${n}_"$REPEAT"_tail.fastq &
			cutadapt -e 0 --minimum-length "$MIN_NONREPETITIVE_LENGTH" -g ^${MAX_HEAD_1} -g ^${MAX_HEAD_2} ./"$SAMPLE_PREFIX"_2/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT_REVCOMP"_head.fastq -o ./"$SAMPLE_PREFIX"_2/03-RemoveRepeats/${n}_"$REPEAT_REVCOMP"_head.fastq &
			wait

			for i in $(seq 1 "$REMOVE_UP_TO")
			do
				cutadapt -e 0 --minimum-length "$MIN_NONREPETITIVE_LENGTH" -a ${MAX_TAIL_1:i}$ -a ${MAX_TAIL_2:i}$ ./"$SAMPLE_PREFIX"_1/03-RemoveRepeats/${n}_"$REPEAT"_tail.fastq -o ./"$SAMPLE_PREFIX"_1/03-RemoveRepeats/${m}_"$REPEAT"_tail.fastq &
				cutadapt -e 0 --minimum-length "$MIN_NONREPETITIVE_LENGTH" -g ^${MAX_HEAD_1:i} -g ^${MAX_HEAD_2:i} ./"$SAMPLE_PREFIX"_1/03-RemoveRepeats/${n}_"$REPEAT_REVCOMP"_head.fastq -o ./"$SAMPLE_PREFIX"_1/03-RemoveRepeats/${m}_"$REPEAT_REVCOMP"_head.fastq &
				cutadapt -e 0 --minimum-length "$MIN_NONREPETITIVE_LENGTH" -a ${MAX_TAIL_1:i}$ -a ${MAX_TAIL_2:i}$ ./"$SAMPLE_PREFIX"_2/03-RemoveRepeats/${n}_"$REPEAT"_tail.fastq -o ./"$SAMPLE_PREFIX"_2/03-RemoveRepeats/${m}_"$REPEAT"_tail.fastq &
				cutadapt -e 0 --minimum-length "$MIN_NONREPETITIVE_LENGTH" -g ^${MAX_HEAD_1:i} -g ^${MAX_HEAD_2:i} ./"$SAMPLE_PREFIX"_2/03-RemoveRepeats/${n}_"$REPEAT_REVCOMP"_head.fastq -o ./"$SAMPLE_PREFIX"_2/03-RemoveRepeats/${m}_"$REPEAT_REVCOMP"_head.fastq &
				wait
				m=$[$m-1] 
				n=$[$n-1] 
			done

			# restore n and m variables
			let i=0	
			let n="$MAX_NUM_REPEATS"*${#REPEAT}
			let m=${n}-1

			# remove temp files
			let start="$MIN_TAIL_LENGTH"+1
			let end=${MAX_NUM_REPEATS}*${#REPEAT}
			for i in $(seq "$start" "$end")
			do
				rm ./"$SAMPLE_PREFIX"_1/03-RemoveRepeats/${i}_"$REPEAT"_tail.fastq &
				rm ./"$SAMPLE_PREFIX"_1/03-RemoveRepeats/${i}_"$REPEAT_REVCOMP"_head.fastq &
				rm ./"$SAMPLE_PREFIX"_2/03-RemoveRepeats/${i}_"$REPEAT"_tail.fastq &
				rm ./"$SAMPLE_PREFIX"_2/03-RemoveRepeats/${i}_"$REPEAT_REVCOMP"_head.fastq &
				wait
			done

			if [ -f "$SAMPLE_PREFIX"_1/03-RemoveRepeats/"$MIN_TAIL_LENGTH"_"$REPEAT"_tail.fastq ] &&
			   [ -f "$SAMPLE_PREFIX"_1/03-RemoveRepeats/"$MIN_TAIL_LENGTH"_"$REPEAT_REVCOMP"_head.fastq ] &&
		       [ -f "$SAMPLE_PREFIX"_2/03-RemoveRepeats/"$MIN_TAIL_LENGTH"_"$REPEAT"_tail.fastq ] &&
		       [ -f "$SAMPLE_PREFIX"_2/03-RemoveRepeats/"$MIN_TAIL_LENGTH"_"$REPEAT_REVCOMP"_head.fastq ]
			then 
				echo "[`date`] Repeat removal complete" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
			else 
				echo "[`date`] ERROR: Repeat removal failed. Check log for details" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
				mv "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt "$SAMPLE_PREFIX"_log_$(date +"%Y-%m-%d_%H:%M:%S").txt
				echo "[`date`] ERROR: Repeat removal failed. Check log for details"
				exit 1
			fi
		fi

		#4.Align reads with repeats removed to genome
		if [ -f "$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT"-tail/accepted_hits.bam ] &&
		   [ -f "$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT_REVCOMP"-head/accepted_hits.bam ] &&
	       [ -f "$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT"-tail/accepted_hits.bam ] &&
	       [ -f "$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT_REVCOMP"-head/accepted_hits.bam ]
		then 
			echo "[`date`] Read alignment to "$REF_SEQ_PREFIX" already complete" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
		else 
			echo "[`date`] Aligning reads with repeats removed to "$REF_SEQ_PREFIX" genome" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
			if [ ! -d "$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT"-tail ]; then mkdir -p "$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT"-tail; fi
			if [ ! -d "$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT_REVCOMP"-head ]; then mkdir -p "$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT_REVCOMP"-head; fi
			if [ ! -d "$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT"-tail ]; then mkdir -p "$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT"-tail; fi
			if [ ! -d "$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT_REVCOMP"-head ]; then mkdir -p "$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT_REVCOMP"-head; fi

			# read alignment
			if [ $NO_NOVEL_JUNCS = "false" ]
			then
				tophat -p "$MAX_CORES" --read-mismatches 0 --read-gap-length 0 --max-deletion-length 0 --max-insertion-length 0 --read-edit-dist 0 --splice-mismatches 0 --segment-mismatches 0 -o ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT"-tail "$REF_FOR_TOPHAT" ./"$SAMPLE_PREFIX"_1/03-RemoveRepeats/"$MIN_TAIL_LENGTH"_"$REPEAT"_tail.fastq
				tophat -p "$MAX_CORES" --read-mismatches 0 --read-gap-length 0 --max-deletion-length 0 --max-insertion-length 0 --read-edit-dist 0 --splice-mismatches 0 --segment-mismatches 0 -o ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT_REVCOMP"-head "$REF_FOR_TOPHAT" ./"$SAMPLE_PREFIX"_1/03-RemoveRepeats/"$MIN_TAIL_LENGTH"_"$REPEAT_REVCOMP"_head.fastq
				
				tophat -p "$MAX_CORES" --read-mismatches 0 --read-gap-length 0 --max-deletion-length 0 --max-insertion-length 0 --read-edit-dist 0 --splice-mismatches 0 --segment-mismatches 0 -o ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT"-tail "$REF_FOR_TOPHAT" ./"$SAMPLE_PREFIX"_2/03-RemoveRepeats/"$MIN_TAIL_LENGTH"_"$REPEAT"_tail.fastq
				tophat -p "$MAX_CORES" --read-mismatches 0 --read-gap-length 0 --max-deletion-length 0 --max-insertion-length 0 --read-edit-dist 0 --splice-mismatches 0 --segment-mismatches 0 -o ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT_REVCOMP"-head "$REF_FOR_TOPHAT" ./"$SAMPLE_PREFIX"_2/03-RemoveRepeats/"$MIN_TAIL_LENGTH"_"$REPEAT_REVCOMP"_head.fastq
			else
				if [ $NO_NOVEL_JUNCS = "true" ]
				then
					tophat -p "$MAX_CORES" --no-novel-juncs --read-mismatches 0 --read-gap-length 0 --max-deletion-length 0 --max-insertion-length 0 --read-edit-dist 0 --splice-mismatches 0 --segment-mismatches 0 -o ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT"-tail "$REF_FOR_TOPHAT" ./"$SAMPLE_PREFIX"_1/03-RemoveRepeats/"$MIN_TAIL_LENGTH"_"$REPEAT"_tail.fastq
					tophat -p "$MAX_CORES" --no-novel-juncs --read-mismatches 0 --read-gap-length 0 --max-deletion-length 0 --max-insertion-length 0 --read-edit-dist 0 --splice-mismatches 0 --segment-mismatches 0 -o ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT_REVCOMP"-head "$REF_FOR_TOPHAT" ./"$SAMPLE_PREFIX"_1/03-RemoveRepeats/"$MIN_TAIL_LENGTH"_"$REPEAT_REVCOMP"_head.fastq
					
					tophat -p "$MAX_CORES" --no-novel-juncs --read-mismatches 0 --read-gap-length 0 --max-deletion-length 0 --max-insertion-length 0 --read-edit-dist 0 --splice-mismatches 0 --segment-mismatches 0 -o ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT"-tail "$REF_FOR_TOPHAT" ./"$SAMPLE_PREFIX"_2/03-RemoveRepeats/"$MIN_TAIL_LENGTH"_"$REPEAT"_tail.fastq
					tophat -p "$MAX_CORES" --no-novel-juncs --read-mismatches 0 --read-gap-length 0 --max-deletion-length 0 --max-insertion-length 0 --read-edit-dist 0 --splice-mismatches 0 --segment-mismatches 0 -o ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT_REVCOMP"-head "$REF_FOR_TOPHAT" ./"$SAMPLE_PREFIX"_2/03-RemoveRepeats/"$MIN_TAIL_LENGTH"_"$REPEAT_REVCOMP"_head.fastq
				fi
			fi

			if [ -f "$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT"-tail/accepted_hits.bam ] &&
			   [ -f "$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT_REVCOMP"-head/accepted_hits.bam ] &&
		       [ -f "$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT"-tail/accepted_hits.bam ] &&
		       [ -f "$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT_REVCOMP"-head/accepted_hits.bam ]
			then 
				echo "[`date`] Read alignment to "$REF_SEQ_PREFIX" complete" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
			else 
				echo "[`date`] ERROR: Read alignment to "$REF_SEQ_PREFIX" failed. Check log for details" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
				mv "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt "$SAMPLE_PREFIX"_log_$(date +"%Y-%m-%d_%H:%M:%S").txt
				echo "[`date`] ERROR: Read alignment to "$REF_SEQ_PREFIX" failed. Check log for details"
				exit 1
			fi
		fi

		#5.Add back repeats to aligned reads
		if [ -f "$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped_reads.fastq ] &&
		   [ -f "$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped_reads.fastq ] &&
	       [ -f "$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped_reads.fastq ] &&
	       [ -f "$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped_reads.fastq ]
		then 
			echo "[`date`] Adding back repeats already complete" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
		else 
			echo "[`date`] Adding back repeats to aligned reads" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
			if [ ! -d "$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/05-AddBackRepeats ]; then mkdir -p "$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/05-AddBackRepeats; fi
			if [ ! -d "$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/05-AddBackRepeats ]; then mkdir -p "$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/05-AddBackRepeats; fi

			# covert bam to sam
			samtools view ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT"-tail/accepted_hits.bam -o ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped.sam &
			samtools view ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT_REVCOMP"-head/accepted_hits.bam -o ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped.sam &
			samtools view ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT"-tail/accepted_hits.bam -o ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped.sam &
			samtools view ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT_REVCOMP"-head/accepted_hits.bam -o ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped.sam &
			wait

			# remove alignments to dummy chr that aren't dummy read
			awk -v rep="$REPEAT" 'BEGIN {FS="\t"; OFS="\t"} ($1 == "dummy"rep".1" && $3 == "dummy_chr_"rep"") || ($1 != "dummy"rep".1" && $3 != "dummy_chr_"rep"") {print $0}' ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped.sam > ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped_cleaned.sam &
			awk -v rep_rev="$REPEAT_REVCOMP" 'BEGIN {FS="\t"; OFS="\t"} ($1 == "dummy"rep_rev".1" && $3 == "dummy_chr_"rep"") || ($1 != "dummy"rep_rev".1" && $3 != "dummy_chr_"rep"") {print $0}' ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped.sam > ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped_cleaned.sam &
			awk -v rep="$REPEAT" 'BEGIN {FS="\t"; OFS="\t"} ($1 == "dummy"rep".2" && $3 == "dummy_chr_"rep"") || ($1 != "dummy"rep".2" && $3 != "dummy_chr_"rep"") {print $0}' ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped.sam > ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped_cleaned.sam &
			awk -v rep_rev="$REPEAT_REVCOMP" 'BEGIN {FS="\t"; OFS="\t"} ($1 == "dummy"rep_rev".2" && $3 == "dummy_chr_"rep"") || ($1 != "dummy"rep_rev".2" && $3 != "dummy_chr_"rep"") {print $0}' ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped.sam > ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped_cleaned.sam &
			wait

			# grab read names and remove duplicates
			awk '{print $1}' ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped_cleaned.sam | awk '!a[$0]++' > ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped_read_names &
			awk '{print $1}' ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped_cleaned.sam | awk '!a[$0]++' > ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped_read_names &
			awk '{print $1}' ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped_cleaned.sam | awk '!a[$0]++' > ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped_read_names &
			awk '{print $1}' ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped_cleaned.sam | awk '!a[$0]++' > ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped_read_names &
			wait

			# print fastq for given reads
			"$PATH_TO_SEQTK"/seqtk subseq ./"$SAMPLE_PREFIX"_1/02-GrabReads/"$SAMPLE_PREFIX"_1_trimmed_"$REPEAT"_dummied.fq ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped_read_names > ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped_reads.fastq &
			"$PATH_TO_SEQTK"/seqtk subseq ./"$SAMPLE_PREFIX"_1/02-GrabReads/"$SAMPLE_PREFIX"_1_trimmed_"$REPEAT"_dummied.fq ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped_read_names > ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped_reads.fastq &
			"$PATH_TO_SEQTK"/seqtk subseq ./"$SAMPLE_PREFIX"_2/02-GrabReads/"$SAMPLE_PREFIX"_2_trimmed_"$REPEAT"_dummied.fq ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped_read_names > ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped_reads.fastq &
			"$PATH_TO_SEQTK"/seqtk subseq ./"$SAMPLE_PREFIX"_2/02-GrabReads/"$SAMPLE_PREFIX"_2_trimmed_"$REPEAT"_dummied.fq ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped_read_names > ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped_reads.fastq &
			wait

			if [ -f "$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped_reads.fastq ] &&
			   [ -f "$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped_reads.fastq ] &&
		       [ -f "$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped_reads.fastq ] &&
		       [ -f "$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped_reads.fastq ]
			then 
				echo "[`date`] Adding back repeats complete" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
			else 
				echo "[`date`] ERROR: Adding back repeats failed. Check log for details" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
				mv "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt "$SAMPLE_PREFIX"_log_$(date +"%Y-%m-%d_%H:%M:%S").txt
				echo "[`date`] ERROR: Adding back repeats failed. Check log for details"
				exit 1
			fi
		fi

		#6.Align reads with repeats added back to genome
		if [ -f "$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT"-tail/accepted_hits.bam ] &&
		   [ -f "$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT_REVCOMP"-head/accepted_hits.bam ] &&
	       [ -f "$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT"-tail/accepted_hits.bam ] &&
	       [ -f "$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT_REVCOMP"-head/accepted_hits.bam ]
		then 
			echo "[`date`] Read alignment to "$REF_SEQ_PREFIX" already complete" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
		else 
			echo "[`date`] Aligning reads with repeats added back to "$REF_SEQ_PREFIX" genome" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
			if [ ! -d "$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT"-tail ]; then mkdir -p "$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT"-tail; fi
			if [ ! -d "$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT_REVCOMP"-head ]; then mkdir -p "$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT_REVCOMP"-head; fi
			if [ ! -d "$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT"-tail ]; then mkdir -p "$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT"-tail; fi
			if [ ! -d "$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT_REVCOMP"-head ]; then mkdir -p "$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT_REVCOMP"-head; fi

			# read alignment
			if [ $NO_NOVEL_JUNCS = "false" ]
			then
				tophat -p "$MAX_CORES" --read-mismatches 1 --read-gap-length 1 --max-deletion-length 1 --max-insertion-length 1 --read-edit-dist 1 --splice-mismatches 1 --segment-mismatches 1 -o ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT"-tail $REF_FOR_TOPHAT ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped_reads.fastq
				tophat -p "$MAX_CORES" --read-mismatches 1 --read-gap-length 1 --max-deletion-length 1 --max-insertion-length 1 --read-edit-dist 1 --splice-mismatches 1 --segment-mismatches 1 -o ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT_REVCOMP"-head $REF_FOR_TOPHAT ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped_reads.fastq
				
				tophat -p "$MAX_CORES" --read-mismatches 1 --read-gap-length 1 --max-deletion-length 1 --max-insertion-length 1 --read-edit-dist 1 --splice-mismatches 1 --segment-mismatches 1 -o ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT"-tail $REF_FOR_TOPHAT ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped_reads.fastq
				tophat -p "$MAX_CORES" --read-mismatches 1 --read-gap-length 1 --max-deletion-length 1 --max-insertion-length 1 --read-edit-dist 1 --splice-mismatches 1 --segment-mismatches 1 -o ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT_REVCOMP"-head $REF_FOR_TOPHAT ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped_reads.fastq
			else
				if [ $NO_NOVEL_JUNCS = "true" ]
				then
					tophat -p "$MAX_CORES" --no-novel-juncs --read-mismatches 1 --read-gap-length 1 --max-deletion-length 1 --max-insertion-length 1 --read-edit-dist 1 --splice-mismatches 1 --segment-mismatches 1 -o ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT"-tail $REF_FOR_TOPHAT ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped_reads.fastq
					tophat -p "$MAX_CORES" --no-novel-juncs --read-mismatches 1 --read-gap-length 1 --max-deletion-length 1 --max-insertion-length 1 --read-edit-dist 1 --splice-mismatches 1 --segment-mismatches 1 -o ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT_REVCOMP"-head $REF_FOR_TOPHAT ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped_reads.fastq
					
					tophat -p "$MAX_CORES" --no-novel-juncs --read-mismatches 1 --read-gap-length 1 --max-deletion-length 1 --max-insertion-length 1 --read-edit-dist 1 --splice-mismatches 1 --segment-mismatches 1 -o ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT"-tail $REF_FOR_TOPHAT ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped_reads.fastq
					tophat -p "$MAX_CORES" --no-novel-juncs --read-mismatches 1 --read-gap-length 1 --max-deletion-length 1 --max-insertion-length 1 --read-edit-dist 1 --splice-mismatches 1 --segment-mismatches 1 -o ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT_REVCOMP"-head $REF_FOR_TOPHAT ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped_reads.fastq
				fi
			fi

			if [ -f "$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT"-tail/accepted_hits.bam ] &&
			   [ -f "$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT_REVCOMP"-head/accepted_hits.bam ] &&
		       [ -f "$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT"-tail/accepted_hits.bam ] &&
		       [ -f "$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT_REVCOMP"-head/accepted_hits.bam ]
			then 
				echo "[`date`] Read alignment to "$REF_SEQ_PREFIX" complete" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
			else 
				echo "[`date`] ERROR: Read alignment to "$REF_SEQ_PREFIX" failed. Check log for details" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
				mv "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt "$SAMPLE_PREFIX"_log_$(date +"%Y-%m-%d_%H:%M:%S").txt
				echo "[`date`] ERROR: Read alignment to "$REF_SEQ_PREFIX" failed. Check log for details"
				exit 1
			fi
		fi

		#7.Print reads that failed to map and get mates
		if [ -f ./results/fastq/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates_read_names ] &&
		   [ -f ./results/fastq/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates.fastq ] &&
		   [ -f ./results/fastq/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates_read_names_with_mate_ids_stripped ] &&
		   [ -f ./results/fastq/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates_with_mates.fastq ]
		then 
			echo "[`date`] Printing nontemplated p"$REPEAT" tail reads already complete" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
		else 
			echo "[`date`] Printing nontemplated p"$REPEAT" tail reads and their mates" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
			if [ ! -d ./results/fastq ]; then mkdir -p ./results/fastq; fi
			if [ ! -d "$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/07-PrintReads ]; then mkdir -p "$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/07-PrintReads; fi
			if [ ! -d "$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/07-PrintReads ]; then mkdir -p "$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/07-PrintReads; fi

			# covert bam to sam
			samtools view ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT"-tail/unmapped.bam -o ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/07-PrintReads/"$REPEAT"_tail_unmapped.sam &
			samtools view ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT_REVCOMP"-head/unmapped.bam -o ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/07-PrintReads/"$REPEAT_REVCOMP"_head_unmapped.sam &
			samtools view ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT"-tail/unmapped.bam -o ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/07-PrintReads/"$REPEAT"_tail_unmapped.sam &
			samtools view ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT_REVCOMP"-head/unmapped.bam -o ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/07-PrintReads/"$REPEAT_REVCOMP"_head_unmapped.sam &
			wait

			# grab read names and remove duplicates
			awk '{print $1}' ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/07-PrintReads/"$REPEAT"_tail_unmapped.sam | awk '!a[$0]++' > ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/07-PrintReads/"$REPEAT"_tail_unmapped_read_names &
			awk '{print $1}' ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/07-PrintReads/"$REPEAT_REVCOMP"_head_unmapped.sam | awk '!a[$0]++' > ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/07-PrintReads/"$REPEAT_REVCOMP"_head_unmapped_read_names &
			awk '{print $1}' ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/07-PrintReads/"$REPEAT"_tail_unmapped.sam | awk '!a[$0]++' > ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/07-PrintReads/"$REPEAT"_tail_unmapped_read_names &
			awk '{print $1}' ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/07-PrintReads/"$REPEAT_REVCOMP"_head_unmapped.sam | awk '!a[$0]++' > ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/07-PrintReads/"$REPEAT_REVCOMP"_head_unmapped_read_names &
			wait

			sed s'/[1]$//' ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/07-PrintReads/"$REPEAT"_tail_unmapped_read_names > ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/07-PrintReads/"$REPEAT"_tail_unmapped_read_names_stripped &
			sed s'/[1]$//' ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/07-PrintReads/"$REPEAT_REVCOMP"_head_unmapped_read_names > ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/07-PrintReads/"$REPEAT_REVCOMP"_head_unmapped_read_names_stripped &
			sed s'/[2]$//' ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/07-PrintReads/"$REPEAT"_tail_unmapped_read_names > ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/07-PrintReads/"$REPEAT"_tail_unmapped_read_names_stripped &
			sed s'/[2]$//' ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/07-PrintReads/"$REPEAT_REVCOMP"_head_unmapped_read_names > ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/07-PrintReads/"$REPEAT_REVCOMP"_head_unmapped_read_names_stripped &
			wait

			# combine read names and remove duplicates
			cat ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/07-PrintReads/"$REPEAT"_tail_unmapped_read_names ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/07-PrintReads/"$REPEAT_REVCOMP"_head_unmapped_read_names ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/07-PrintReads/"$REPEAT"_tail_unmapped_read_names ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/07-PrintReads/"$REPEAT_REVCOMP"_head_unmapped_read_names | sort | uniq > ./results/fastq/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates_read_names &
			cat ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/07-PrintReads/"$REPEAT"_tail_unmapped_read_names_stripped ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/07-PrintReads/"$REPEAT_REVCOMP"_head_unmapped_read_names_stripped ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/07-PrintReads/"$REPEAT"_tail_unmapped_read_names_stripped ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/07-PrintReads/"$REPEAT_REVCOMP"_head_unmapped_read_names_stripped | sort | uniq > ./results/fastq/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates_read_names_with_mate_ids_stripped &
			wait

			# if read name files are empty (i.e. no candidates were identified), then exit script.
			if [[ $(wc ./results/fastq/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates_read_names | awk -F " " '{print $1}') = 0 ]]
			then
				# STATISTICS
				NUM_READS_BEFORE_TRIMMING_R1=$(grep "Total reads processed" ./"$SAMPLE_PREFIX"_1/01-TrimReads/"$SAMPLE_PREFIX"_1.fastq.gz_trimming_report.txt | awk '{print $4}' | sed 's/,//g')
				NUM_READS_BEFORE_TRIMMING_R2=$(grep "Total reads processed" ./"$SAMPLE_PREFIX"_2/01-TrimReads/"$SAMPLE_PREFIX"_2.fastq.gz_trimming_report.txt | awk '{print $4}' | sed 's/,//g')
				let NUM_READS_BEFORE_TRIMMING="$NUM_READS_BEFORE_TRIMMING_R1"+"$NUM_READS_BEFORE_TRIMMING_R2"

				NUM_READS_AFTER_TRIMMING_R1=$(grep "Number of input reads" ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/star/"$SAMPLE_PREFIX"_Log.final.out | awk -F "\t" '{print $NF}')
				NUM_READS_AFTER_TRIMMING_R2=$(grep "Number of input reads" ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/star/"$SAMPLE_PREFIX"_Log.final.out | awk -F "\t" '{print $NF}')
				let NUM_READS_AFTER_TRIMMING="$NUM_READS_AFTER_TRIMMING_R1"+"$NUM_READS_AFTER_TRIMMING_R2"

				PERCENT_REMAINING_AFTER_TRIM=$(awk -v TOTAL="$NUM_READS_BEFORE_TRIMMING" -v TRIM="$NUM_READS_AFTER_TRIMMING" 'BEGIN{printf("%.2f\n",(TRIM/TOTAL)*100)}')

				NUM_READS_THAT_UNIQUELY_MAP_TO_GENOME_R1=$(grep "Uniquely mapped reads number" ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/star/"$SAMPLE_PREFIX"_Log.final.out | awk -F "\t" '{print $NF}')
				NUM_READS_THAT_UNIQUELY_MAP_TO_GENOME_R2=$(grep "Uniquely mapped reads number" ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/star/"$SAMPLE_PREFIX"_Log.final.out | awk -F "\t" '{print $NF}')
				let NUM_READS_THAT_UNIQUELY_MAP_TO_GENOME="$NUM_READS_THAT_UNIQUELY_MAP_TO_GENOME_R1"+"$NUM_READS_THAT_UNIQUELY_MAP_TO_GENOME_R2"

				PERCENT_OF_TRIMMED_READS_THAT_MAP=$(awk -v TRIM="$NUM_READS_AFTER_TRIMMING" -v MAP="$NUM_READS_THAT_UNIQUELY_MAP_TO_GENOME" 'BEGIN{printf("%.2f\n",(MAP/TRIM)*100)}')

				NUM_READS_WITH_REPEAT_TAILS=$(cat \
				    ./"$SAMPLE_PREFIX"_1/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT"_tail.fastq \
				    ./"$SAMPLE_PREFIX"_1/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT_REVCOMP"_head.fastq \
				    ./"$SAMPLE_PREFIX"_2/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT"_tail.fastq \
				    ./"$SAMPLE_PREFIX"_2/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT_REVCOMP"_head.fastq \
				    | awk 'NR % 4 == 1' | grep -v dummy | sort | uniq | wc -l)

				NUM_READS_WITH_REPEAT_TAILS_THAT_MAP_TO_GENOME=$(cat \
					./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped_read_names \
					./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped_read_names \
					./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped_read_names \
					./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped_read_names \
					| grep -v dummy | wc -l)
					#read names from SAM files do not have read 1 or read 2 identifiers, so don't uniqify

				PERCENT_OF_REPEAT_READS_THAT_MAP=$(awk -v REPEATS="$NUM_READS_WITH_REPEAT_TAILS" -v MAP="$NUM_READS_WITH_REPEAT_TAILS_THAT_MAP_TO_GENOME" 'BEGIN{printf("%.2f\n",(MAP/REPEATS)*100)}')

				NUM_MAPPED_READS_NONTEMPLATED=$(wc -l ./results/fastq/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates_read_names | awk '{print $1}')

				PERCENT_OF_MAPPED_READS_THAT_ARE_NONTEMPLATED=$(awk -v MAP="$NUM_READS_WITH_REPEAT_TAILS_THAT_MAP_TO_GENOME" -v NONTEMP="$NUM_MAPPED_READS_NONTEMPLATED" 'BEGIN{printf("%.2f\n",(NONTEMP/MAP)*100)}')

				STATISTICS="=== STATISTICS ===

Total reads processed: 								$(printf "%'.f\n" $NUM_READS_BEFORE_TRIMMING)
Reads remaining after adapter and quality trim:		$(printf "%'.f\n" $NUM_READS_AFTER_TRIMMING) ("$PERCENT_REMAINING_AFTER_TRIM"% of total reads)
Reads that map to "$REF_SEQ_PREFIX" genome:				$(printf "%'.f\n" $NUM_READS_THAT_UNIQUELY_MAP_TO_GENOME) ("$PERCENT_OF_TRIMMED_READS_THAT_MAP"% of trimmed reads)

p"$REPEAT" reads with at least "$MIN_NUM_REPEATS" repeats: 				$(printf "%'.f\n" $NUM_READS_WITH_REPEAT_TAILS)
p"$REPEAT" reads that map to "$REF_SEQ_PREFIX" genome:			$(printf "%'.f\n" $NUM_READS_WITH_REPEAT_TAILS_THAT_MAP_TO_GENOME) ("$PERCENT_OF_REPEAT_READS_THAT_MAP"% of p"$REPEAT" reads)
p"$REPEAT" reads that are nontemplated: 					$(printf "%'.f\n" $NUM_MAPPED_READS_NONTEMPLATED) ("$PERCENT_OF_MAPPED_READS_THAT_ARE_NONTEMPLATED"% of p"$REPEAT" mapping reads)
"

				echo "[`date`] Failed to find any p"$REPEAT" reads that are nontemplated" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
				echo "[`date`] FIND NONTEMPLATED "$REPEAT" TAILS COMPLETE" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
				printf -- "%s\n" " " >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt 
				echo "$STATISTICS" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt

				mv "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log_final.txt
				echo "[`date`] FIND NONTEMPLATED "$REPEAT" TAILS COMPLETE. Check "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log_final.txt for summary"
				exit 0
			fi

			# grab reads
			echo | tr -d '\n' > ./results/fastq/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates.fastq
			echo | tr -d '\n' > ./results/fastq/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates_with_mates.fastq

			for i in $(cat ./results/fastq/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates_read_names)
			do
				grep -A 3 @"$i" ./"$SAMPLE_PREFIX"_1/02-GrabReads/"$SAMPLE_PREFIX"_1_trimmed_"$REPEAT"_dummied.fq >> ./results/fastq/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates.fastq
				grep -A 3 @"$i" ./"$SAMPLE_PREFIX"_2/02-GrabReads/"$SAMPLE_PREFIX"_2_trimmed_"$REPEAT"_dummied.fq >> ./results/fastq/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates.fastq
			done &
			for i in $(cat ./results/fastq/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates_read_names_with_mate_ids_stripped)
			do
				grep -A 3 @"$i"1 ./"$SAMPLE_PREFIX"_1/02-GrabReads/"$SAMPLE_PREFIX"_1_trimmed_"$REPEAT"_dummied.fq >> ./results/fastq/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates_with_mates.fastq
				grep -A 3 @"$i"2 ./"$SAMPLE_PREFIX"_2/02-GrabReads/"$SAMPLE_PREFIX"_2_trimmed_"$REPEAT"_dummied.fq >> ./results/fastq/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates_with_mates.fastq
			done &
			wait

			# copy 
			if [ ! -d ./../candidates ]; then mkdir -p ./../candidates; fi
			cp ./results/fastq/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates.fastq ./../candidates &
			cp ./results/fastq/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates_with_mates.fastq ./../candidates &
			wait

			if [ -f ./results/fastq/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates_read_names ] &&
			   [ -f ./results/fastq/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates.fastq ] &&
			   [ -f ./results/fastq/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates_read_names_with_mate_ids_stripped ] &&
			   [ -f ./results/fastq/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates_with_mates.fastq ]
			then 
				echo "[`date`] Printing nontemplated p"$REPEAT" tail reads complete" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
			else 
				echo "[`date`] ERROR: Printing nontemplated p"$REPEAT" tail reads failed. Check log for details" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
				mv "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt "$SAMPLE_PREFIX"_log_$(date +"%Y-%m-%d_%H:%M:%S").txt
				echo "[`date`] ERROR: Printing nontemplated p"$REPEAT" tail reads failed. Check log for details"
				exit 1
			fi
		fi
	
		#8.Align candidates to genome
		if [ -f ./results/alignments/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates.bam ] &&
		   [ -f ./results/alignments/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates.bam.bai ]
		then 
			echo "[`date`] p"$REPEAT" tail candidate read alignment to "$REF_SEQ_PREFIX" genome already complete" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
		else 
			echo "[`date`] Aligning nontemplated p"$REPEAT" tail candidates to "$REF_SEQ_PREFIX" genome" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
			if [ ! -d ./results/alignments ]; then mkdir -p ./results/alignments; fi
			if [ ! -d "$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/08-AlignCandidates ]; then mkdir -p "$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/08-AlignCandidates; fi
			if [ ! -d "$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/08-AlignCandidates ]; then mkdir -p "$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/08-AlignCandidates; fi

			# covert bam to sam
			samtools view ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT"-tail/accepted_hits.bam -o ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/08-AlignCandidates/"$REPEAT"_tail_mapped.sam &
			samtools view ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT_REVCOMP"-head/accepted_hits.bam -o ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/08-AlignCandidates/"$REPEAT_REVCOMP"_head_mapped.sam &
			samtools view ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT"-tail/accepted_hits.bam -o ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/08-AlignCandidates/"$REPEAT"_tail_mapped.sam &
			samtools view ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT_REVCOMP"-head/accepted_hits.bam -o ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/08-AlignCandidates/"$REPEAT_REVCOMP"_head_mapped.sam &
			wait

			# grab read names and remove duplicates
			awk '{print $1}' ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/08-AlignCandidates/"$REPEAT"_tail_mapped.sam | awk '!a[$0]++' > ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/08-AlignCandidates/"$REPEAT"_tail_mapped_read_names &
			awk '{print $1}' ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/08-AlignCandidates/"$REPEAT_REVCOMP"_head_mapped.sam | awk '!a[$0]++' > ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/08-AlignCandidates/"$REPEAT_REVCOMP"_head_mapped_read_names &
			awk '{print $1}' ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/08-AlignCandidates/"$REPEAT"_tail_mapped.sam | awk '!a[$0]++' > ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/08-AlignCandidates/"$REPEAT"_tail_mapped_read_names &
			awk '{print $1}' ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/08-AlignCandidates/"$REPEAT_REVCOMP"_head_mapped.sam | awk '!a[$0]++' > ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/08-AlignCandidates/"$REPEAT_REVCOMP"_head_mapped_read_names &
			wait

			# remove reads that aligned with repeats added back from the first alignment with reads with repeats removed
			for i in $(cat ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/08-AlignCandidates/"$REPEAT"_tail_mapped_read_names); do grep -v "$i" ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped_cleaned.sam > ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/08-AlignCandidates/"$REPEAT"_tail_mapped_cleaned_geno_encode_removed.sam; done &
			for i in $(cat ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/08-AlignCandidates/"$REPEAT_REVCOMP"_head_mapped_read_names); do grep -v "$i" ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped_cleaned.sam > ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/08-AlignCandidates/"$REPEAT_REVCOMP"_head_mapped_cleaned_geno_encode_removed.sam; done &
			for i in $(cat ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/08-AlignCandidates/"$REPEAT"_tail_mapped_read_names); do grep -v "$i" ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped_cleaned.sam > ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/08-AlignCandidates/"$REPEAT"_tail_mapped_cleaned_geno_encode_removed.sam; done &
			for i in $(cat ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/08-AlignCandidates/"$REPEAT_REVCOMP"_head_mapped_read_names); do grep -v "$i" ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped_cleaned.sam > ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/08-AlignCandidates/"$REPEAT_REVCOMP"_head_mapped_cleaned_geno_encode_removed.sam; done &
			wait

			# combine alignments and add header
			samtools view -H ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT"-tail/accepted_hits.bam | cat - ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/08-AlignCandidates/"$REPEAT"_tail_mapped_cleaned_geno_encode_removed.sam ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/08-AlignCandidates/"$REPEAT_REVCOMP"_head_mapped_cleaned_geno_encode_removed.sam ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/08-AlignCandidates/"$REPEAT"_tail_mapped_cleaned_geno_encode_removed.sam ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/08-AlignCandidates/"$REPEAT_REVCOMP"_head_mapped_cleaned_geno_encode_removed.sam > ./results/alignments/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates.sam

			# sort sam and output as bam
			samtools sort ./results/alignments/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates.sam -O BAM -o ./results/alignments/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates.bam
			
#			if [ $FEATURE_COUNTS = "true" ]
#			then
#				"$PATH_TO_SUBREAD"/featureCounts \
#				-a "$ANNOTATION" \
#				-o ./results/alignments/"$SAMPLE_PREFIX"_p"$REPEAT"_counts.txt \
#				-O \
#				./results/alignments/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates.bam
#			
#				awk '$7 != "0" {print $0}' ./results/alignments/"$SAMPLE_PREFIX"_p"$REPEAT"_counts.txt > ./results/alignments/"$SAMPLE_PREFIX"_p"$REPEAT"_non_zero_counts.txt
#			fi

			# index bam
			samtools index ./results/alignments/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates.bam

			# copy
			if [ ! -d ./../candidates ]; then mkdir -p ./../candidates; fi
			cp ./results/alignments/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates.bam ./../candidates &
			cp ./results/alignments/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates.bam.bai ./../candidates &
			wait

			if [ -f ./results/alignments/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates.bam ] &&
			   [ -f ./results/alignments/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates.bam.bai ]
			then 
				echo "[`date`] Aligning nontemplated p"$REPEAT" tail candidates to "$REF_SEQ_PREFIX" genome complete" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
			else 
				echo "[`date`] ERROR: Aligning nontemplated p"$REPEAT" tail candidates to "$REF_SEQ_PREFIX" genome failed. Check log for details" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
				mv "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt "$SAMPLE_PREFIX"_log_$(date +"%Y-%m-%d_%H:%M:%S").txt
				echo "[`date`] ERROR: Aligning nontemplated p"$REPEAT" tail candidates to "$REF_SEQ_PREFIX" genome failed. Check log for details"
				exit 1
			fi
		fi

		#9.Align candidates to genome with their mates
		if [ -f ./results/alignments/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates_with_mates.bam ] &&
		   [ -f ./results/alignments/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates_with_mates.bam.bai ]
		then 
			echo "[`date`] p"$REPEAT" tail candidate read alignment with mates to "$REF_SEQ_PREFIX" genome already complete" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
		else 
			echo "[`date`] Aligning nontemplated p"$REPEAT" tail candidates with mates to "$REF_SEQ_PREFIX" genome" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
			if [ ! -d ./results/alignments ]; then mkdir -p ./results/alignments; fi
			if [ ! -d "$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/09-AlignCandidatesWithMates-"$REPEAT"-tail ]; then mkdir -p "$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/09-AlignCandidatesWithMates-"$REPEAT"-tail; fi

			# remove repeats
			cutadapt -e 0 --minimum-length "$MIN_NONREPETITIVE_LENGTH" -a ${MAX_TAIL_1}$ -a ${MAX_TAIL_2}$ -g ^${MAX_HEAD_1} -g ^${MAX_HEAD_2} ./results/fastq/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates_with_mates.fastq -o ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/09-AlignCandidatesWithMates-"$REPEAT"-tail/${n}_"$REPEAT"_candidates_interleaved.fastq

			for i in $(seq 1 "$REMOVE_UP_TO")
			do
				cutadapt -e 0 --minimum-length "$MIN_NONREPETITIVE_LENGTH" -a ${MAX_TAIL_1:i}$ -a ${MAX_TAIL_2:i}$ -g ^${MAX_HEAD_1:i} -g ^${MAX_HEAD_2:i} ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/09-AlignCandidatesWithMates-"$REPEAT"-tail/${n}_"$REPEAT"_candidates_interleaved.fastq -o ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/09-AlignCandidatesWithMates-"$REPEAT"-tail/${m}_"$REPEAT"_candidates_interleaved.fastq
				m=$[$m-1] 
				n=$[$n-1] 
			done

			# remove temp files
			let start="$MIN_TAIL_LENGTH"+1
			let end=${MAX_NUM_REPEATS}*${#REPEAT}
			for i in $(seq "$start" "$end")
			do
				rm ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/09-AlignCandidatesWithMates-"$REPEAT"-tail/${i}_"$REPEAT"_candidates_interleaved.fastq
			done

			# split interleaved reads
			"$PATH_TO_SEQTK"/seqtk seq -1 ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/09-AlignCandidatesWithMates-"$REPEAT"-tail/"$MIN_TAIL_LENGTH"_"$REPEAT"_candidates_interleaved.fastq > ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/09-AlignCandidatesWithMates-"$REPEAT"-tail/"$MIN_TAIL_LENGTH"_"$REPEAT"_candidates_1.fastq &
			"$PATH_TO_SEQTK"/seqtk seq -2 ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/09-AlignCandidatesWithMates-"$REPEAT"-tail/"$MIN_TAIL_LENGTH"_"$REPEAT"_candidates_interleaved.fastq > ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/09-AlignCandidatesWithMates-"$REPEAT"-tail/"$MIN_TAIL_LENGTH"_"$REPEAT"_candidates_2.fastq &
			wait

			# read alignment
			if [ $NO_NOVEL_JUNCS = "false" ]
			then
				tophat -p "$MAX_CORES" --read-mismatches 1 --read-gap-length 1 --max-deletion-length 1 --max-insertion-length 1 --read-edit-dist 1 --splice-mismatches 1 --segment-mismatches 1 -o ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/09-AlignCandidatesWithMates-"$REPEAT"-tail "$REF_FOR_TOPHAT" ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/09-AlignCandidatesWithMates-"$REPEAT"-tail/"$MIN_TAIL_LENGTH"_"$REPEAT"_candidates_1.fastq ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/09-AlignCandidatesWithMates-"$REPEAT"-tail/"$MIN_TAIL_LENGTH"_"$REPEAT"_candidates_2.fastq
			else
				if [ $NO_NOVEL_JUNCS = "true" ]
				then
					tophat -p "$MAX_CORES" --no-novel-juncs --read-mismatches 1 --read-gap-length 1 --max-deletion-length 1 --max-insertion-length 1 --read-edit-dist 1 --splice-mismatches 1 --segment-mismatches 1 -o ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/09-AlignCandidatesWithMates-"$REPEAT"-tail "$REF_FOR_TOPHAT" ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/09-AlignCandidatesWithMates-"$REPEAT"-tail/"$MIN_TAIL_LENGTH"_"$REPEAT"_candidates_1.fastq ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/09-AlignCandidatesWithMates-"$REPEAT"-tail/"$MIN_TAIL_LENGTH"_"$REPEAT"_candidates_2.fastq
				fi
			fi

			# copy
			cp ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/09-AlignCandidatesWithMates-"$REPEAT"-tail/accepted_hits.bam ./results/alignments/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates_with_mates.bam

			# index
			samtools index ./results/alignments/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates_with_mates.bam

			# copy
			if [ ! -d ./../candidates ]; then mkdir -p ./../candidates; fi
			cp ./results/alignments/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates_with_mates.bam ./../candidates &
			cp ./results/alignments/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates_with_mates.bam.bai ./../candidates &
			wait

			if [ -f ./results/alignments/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates_with_mates.bam ] &&
			   [ -f ./results/alignments/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates_with_mates.bam.bai ]
			then 
				echo "[`date`] Aligning nontemplated p"$REPEAT" tail candidates with mates to "$REF_SEQ_PREFIX" genome complete" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
			else 
				echo "[`date`] ERROR: Aligning nontemplated p"$REPEAT" tail candidates with mates to "$REF_SEQ_PREFIX" genome failed. Check log for details" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
				mv "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log_$(date +"%Y-%m-%d_%H:%M:%S").txt
				echo "[`date`] ERROR: Aligning nontemplated p"$REPEAT" tail candidates with mates to "$REF_SEQ_PREFIX" genome failed. Check log for details"
				exit 1
			fi
		fi
	
		# STATISTICS
		NUM_READS_BEFORE_TRIMMING_R1=$(grep "Total reads processed" ./"$SAMPLE_PREFIX"_1/01-TrimReads/"$SAMPLE_PREFIX"_1.fastq.gz_trimming_report.txt | awk '{print $4}' | sed 's/,//g')
		NUM_READS_BEFORE_TRIMMING_R2=$(grep "Total reads processed" ./"$SAMPLE_PREFIX"_2/01-TrimReads/"$SAMPLE_PREFIX"_2.fastq.gz_trimming_report.txt | awk '{print $4}' | sed 's/,//g')
		let NUM_READS_BEFORE_TRIMMING="$NUM_READS_BEFORE_TRIMMING_R1"+"$NUM_READS_BEFORE_TRIMMING_R2"

		NUM_READS_AFTER_TRIMMING_R1=$(grep "Number of input reads" ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/star/"$SAMPLE_PREFIX"_Log.final.out | awk -F "\t" '{print $NF}')
		NUM_READS_AFTER_TRIMMING_R2=$(grep "Number of input reads" ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/star/"$SAMPLE_PREFIX"_Log.final.out | awk -F "\t" '{print $NF}')
		let NUM_READS_AFTER_TRIMMING="$NUM_READS_AFTER_TRIMMING_R1"+"$NUM_READS_AFTER_TRIMMING_R2"

		PERCENT_REMAINING_AFTER_TRIM=$(awk -v TOTAL="$NUM_READS_BEFORE_TRIMMING" -v TRIM="$NUM_READS_AFTER_TRIMMING" 'BEGIN{printf("%.2f\n",(TRIM/TOTAL)*100)}')

		NUM_READS_THAT_UNIQUELY_MAP_TO_GENOME_R1=$(grep "Uniquely mapped reads number" ./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/star/"$SAMPLE_PREFIX"_Log.final.out | awk -F "\t" '{print $NF}')
		NUM_READS_THAT_UNIQUELY_MAP_TO_GENOME_R2=$(grep "Uniquely mapped reads number" ./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/star/"$SAMPLE_PREFIX"_Log.final.out | awk -F "\t" '{print $NF}')
		let NUM_READS_THAT_UNIQUELY_MAP_TO_GENOME="$NUM_READS_THAT_UNIQUELY_MAP_TO_GENOME_R1"+"$NUM_READS_THAT_UNIQUELY_MAP_TO_GENOME_R2"

		PERCENT_OF_TRIMMED_READS_THAT_MAP=$(awk -v TRIM="$NUM_READS_AFTER_TRIMMING" -v MAP="$NUM_READS_THAT_UNIQUELY_MAP_TO_GENOME" 'BEGIN{printf("%.2f\n",(MAP/TRIM)*100)}')

		NUM_READS_WITH_REPEAT_TAILS=$(cat \
		    ./"$SAMPLE_PREFIX"_1/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT"_tail.fastq \
		    ./"$SAMPLE_PREFIX"_1/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT_REVCOMP"_head.fastq \
		    ./"$SAMPLE_PREFIX"_2/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT"_tail.fastq \
		    ./"$SAMPLE_PREFIX"_2/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT_REVCOMP"_head.fastq \
		    | awk 'NR % 4 == 1' | grep -v dummy | sort | uniq | wc -l)

		NUM_READS_WITH_REPEAT_TAILS_THAT_MAP_TO_GENOME=$(cat \
			./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped_read_names \
			./"$SAMPLE_PREFIX"_1/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped_read_names \
			./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped_read_names \
			./"$SAMPLE_PREFIX"_2/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped_read_names \
			| grep -v dummy | wc -l)
			#read names from SAM files do not have read 1 or read 2 identifiers, so don't uniqify

		PERCENT_OF_REPEAT_READS_THAT_MAP=$(awk -v REPEATS="$NUM_READS_WITH_REPEAT_TAILS" -v MAP="$NUM_READS_WITH_REPEAT_TAILS_THAT_MAP_TO_GENOME" 'BEGIN{printf("%.2f\n",(MAP/REPEATS)*100)}')

		NUM_MAPPED_READS_NONTEMPLATED=$(wc -l ./results/fastq/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates_read_names | awk '{print $1}')

		PERCENT_OF_MAPPED_READS_THAT_ARE_NONTEMPLATED=$(awk -v MAP="$NUM_READS_WITH_REPEAT_TAILS_THAT_MAP_TO_GENOME" -v NONTEMP="$NUM_MAPPED_READS_NONTEMPLATED" 'BEGIN{printf("%.2f\n",(NONTEMP/MAP)*100)}')

		STATISTICS="=== STATISTICS ===

Total reads processed: 								$(printf "%'.f\n" $NUM_READS_BEFORE_TRIMMING)
Reads remaining after adapter and quality trim:		$(printf "%'.f\n" $NUM_READS_AFTER_TRIMMING) ("$PERCENT_REMAINING_AFTER_TRIM"% of total reads)
Reads that map to "$REF_SEQ_PREFIX" genome:				$(printf "%'.f\n" $NUM_READS_THAT_UNIQUELY_MAP_TO_GENOME) ("$PERCENT_OF_TRIMMED_READS_THAT_MAP"% of trimmed reads)

p"$REPEAT" reads with at least "$MIN_NUM_REPEATS" repeats: 				$(printf "%'.f\n" $NUM_READS_WITH_REPEAT_TAILS)
p"$REPEAT" reads that map to "$REF_SEQ_PREFIX" genome:			$(printf "%'.f\n" $NUM_READS_WITH_REPEAT_TAILS_THAT_MAP_TO_GENOME) ("$PERCENT_OF_REPEAT_READS_THAT_MAP"% of p"$REPEAT" reads)
p"$REPEAT" reads that are nontemplated: 					$(printf "%'.f\n" $NUM_MAPPED_READS_NONTEMPLATED) ("$PERCENT_OF_MAPPED_READS_THAT_ARE_NONTEMPLATED"% of p"$REPEAT" mapping reads)
"

		rm ./results/alignments/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates.sam
		rm ./results/fastq/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates_read_names
		rm ./results/fastq/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates_read_names_with_mate_ids_stripped
		echo "[`date`] FIND NONTEMPLATED "$REPEAT" TAILS COMPLETE" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
		printf -- "%s\n" " " >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt 
		echo "$STATISTICS" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt

		mv "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log_final.txt
		echo "[`date`] FIND NONTEMPLATED "$REPEAT" TAILS COMPLETE. Check "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log_final.txt for summary"
		exit 0
	fi

#single-end read workflow
	if [ $READ_TYPE = "single-end" ]
	then
		#0.make analysis directories
		if [ ! -d "$SAMPLE_PREFIX" ]; then mkdir -p "$SAMPLE_PREFIX"; fi

		#1.remove adapters and low quality bases
		if [ -f "$SAMPLE_PREFIX"/01-TrimReads/"$SAMPLE_PREFIX"_trimmed.fq ]
		then 
			echo "[`date`] Read trimming already complete" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
		else 
			if [ -d "$SAMPLE_PREFIX"/01-TrimReads ]
			then
				echo "[`date`] Read trimming is incomplete... Removing "$SAMPLE_PREFIX"/01-TrimReads directory and starting over." >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
				rm -rf "$SAMPLE_PREFIX"/01-TrimReads
			fi
			echo "[`date`] Trimming reads" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
			mkdir -p "$SAMPLE_PREFIX"/01-TrimReads
			
			# trimming reads with mostly default parameters
			trim_galore -q 20 --dont_gzip --trim-n -o ./"$SAMPLE_PREFIX"/01-TrimReads "$SAMPLE".fastq.gz

			if [ -f "$SAMPLE_PREFIX"/01-TrimReads/"$SAMPLE_PREFIX"_trimmed.fq ]
			then 
				echo "[`date`] Read trimming complete" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
			else 
				echo "[`date`] ERROR: Read trimming failed. Check log for details" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
				mv "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt "$SAMPLE_PREFIX"_log_$(date +"%Y-%m-%d_%H:%M:%S").txt
				echo "[`date`] ERROR: Read trimming failed. Check log for details"
				exit 1
			fi
		fi

		# align trimmed reads to genome
		if [ -f "$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/star/"$SAMPLE_PREFIX"_Log.final.out ]
		then 
			echo "[`date`] Trimmed read alignemnt to "$REF_SEQ_PREFIX" genome already complete" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
		else 
			echo "[`date`] Aligning trimmed reads to "$REF_SEQ_PREFIX" genome" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
			if [ ! -d references/"$REF_SEQ_PREFIX"_STAR_genomeDir ]; then mkdir -p references/"$REF_SEQ_PREFIX"_STAR_genomeDir; fi
			if [ ! -d "$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/star ]; then mkdir -p "$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/star; fi

			if [ ! -f references/"$REF_SEQ_PREFIX"_STAR_genomeDir/SAindex ] &&
			   [ ! -f references/"$REF_SEQ_PREFIX"_STAR_genomeDir/SA ] &&
			   [ ! -f references/"$REF_SEQ_PREFIX"_STAR_genomeDir/genomeParameters.txt ] &&
			   [ ! -f references/"$REF_SEQ_PREFIX"_STAR_genomeDir/Genome ] &&
			   [ ! -f references/"$REF_SEQ_PREFIX"_STAR_genomeDir/chrStart.txt ] &&
			   [ ! -f references/"$REF_SEQ_PREFIX"_STAR_genomeDir/chrNameLength.txt ] &&
			   [ ! -f references/"$REF_SEQ_PREFIX"_STAR_genomeDir/chrName.txt ] &&
			   [ ! -f references/"$REF_SEQ_PREFIX"_STAR_genomeDir/chrLength.txt ]
			then
				REF_SEQ_LENGTH=$("$PATH_TO_SEQTK"/seqtk comp "$REFERENCE_GENOME" | awk -F "\t" '{print $2}' | paste -sd+ - | bc)
				Nbases=$(echo "(l("$REF_SEQ_LENGTH")/l(2))/2 -1" | bc -l)
				printf -v NbasesRounded %.0f "$Nbases"

				if [[ "$NbasesRounded" -lt 14 ]]
				then 
					genomeSAindexNbases="$NbasesRounded"
				else
					genomeSAindexNbases=14
				fi 

				STAR \
				--runMode genomeGenerate \
				--genomeSAindexNbases "$genomeSAindexNbases" \
				--genomeDir ./references/"$REF_SEQ_PREFIX"_STAR_genomeDir/ \
				--genomeFastaFiles "$REFERENCE_GENOME"

				mv Log.out ./references/"$REF_SEQ_PREFIX"_STAR_genomeDir/
			fi

			STAR \
			--runMode alignReads \
			--runThreadN "$MAX_CORES" \
			--genomeDir ./references/"$REF_SEQ_PREFIX"_STAR_genomeDir/ \
			--outSAMtype BAM Unsorted \
			--outFileNamePrefix ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/star/"$SAMPLE_PREFIX"_ \
			--readFilesIn ./"$SAMPLE_PREFIX"/01-TrimReads/"$SAMPLE_PREFIX"_trimmed.fq

			if [ -f "$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/star/"$SAMPLE_PREFIX"_Log.final.out ]
			then 
				echo "[`date`] Trimmed read alignemnt to "$REF_SEQ_PREFIX" genome complete" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
			else 
				echo "[`date`] ERROR: Trimmed read alignemnt to "$REF_SEQ_PREFIX" genome failed. Check log for details" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
				mv "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt "$SAMPLE_PREFIX"_log_$(date +"%Y-%m-%d_%H:%M:%S").txt
				echo "[`date`] ERROR: Trimmed read alignemnt to "$REF_SEQ_PREFIX" genome failed. Check log for details"
				exit 1
			fi
		fi

		#2.Select reads with desired repeat and length
		if [ -f "$SAMPLE_PREFIX"/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT"_tail.fastq ] &&
		   [ -f "$SAMPLE_PREFIX"/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT_REVCOMP"_head.fastq ]
		then 
			echo "[`date`] Read selection already complete" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
		else 
			echo "[`date`] Selecting reads with at least "$MIN_NUM_REPEATS" repeats of "$REPEAT_FRAME_1"/"$REPEAT_FRAME_2" at their 3' end" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
			if [ ! -d "$SAMPLE_PREFIX"/02-GrabReads ]; then mkdir -p "$SAMPLE_PREFIX"/02-GrabReads; fi

			# adding dummy reads to the trimmed reads from the previous step
			DUMMY_READ_SEQ_TAIL=TACATCGAGACGGAGCAGGCCCGTTCCAACATTCGTGAGGACGTTGCAAAGGAGGAAGTTCGTCATATTA"$TAIL_1"
			DUMMY_READ_SEQ_HEAD="$HEAD_1"TAATATGACGAACTTCCTCCTTTGCAACGTCCTCACGAATGTTGGAACGGGCCTGCTCCGTCTCGATGTA
			DUMMY_READ_LEN=${#DUMMY_READ_SEQ_TAIL}

			for i in $(seq 1 "$DUMMY_READ_LEN")
			do
				echo I | tr -d '\n' >> DUMMY_READ_QUALITY
			done
			DUMMY_READ_QUAL=$(cat DUMMY_READ_QUALITY)
			rm DUMMY_READ_QUALITY

			printf -- "%s\n" ""@dummy"$REPEAT".1 "$DUMMY_READ_SEQ_TAIL" + "$DUMMY_READ_QUAL" @dummy"$REPEAT_REVCOMP".1 "$DUMMY_READ_SEQ_HEAD" + "$DUMMY_READ_QUAL""" | cat ./"$SAMPLE_PREFIX"/01-TrimReads/"$SAMPLE_PREFIX"_trimmed.fq - > ./"$SAMPLE_PREFIX"/02-GrabReads/"$SAMPLE_PREFIX"_trimmed_"$REPEAT"_dummied.fq

			# grabbing reads with repeats, frame 1
			grep -i --no-group-separator ${TAIL_1}$ -B 1 -A 2 ./"$SAMPLE_PREFIX"/02-GrabReads/"$SAMPLE_PREFIX"_trimmed_"$REPEAT"_dummied.fq > ./"$SAMPLE_PREFIX"/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT"_tail.fastq &
			grep -i --no-group-separator ^${HEAD_1} -B 1 -A 2 ./"$SAMPLE_PREFIX"/02-GrabReads/"$SAMPLE_PREFIX"_trimmed_"$REPEAT"_dummied.fq > ./"$SAMPLE_PREFIX"/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT_REVCOMP"_head.fastq &
			wait

			# grabbing reads with repeats, frame 2
			grep -i --no-group-separator ${TAIL_2}$ -B 1 -A 2 ./"$SAMPLE_PREFIX"/02-GrabReads/"$SAMPLE_PREFIX"_trimmed_"$REPEAT"_dummied.fq >> ./"$SAMPLE_PREFIX"/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT"_tail.fastq &
			grep -i --no-group-separator ^${HEAD_2} -B 1 -A 2 ./"$SAMPLE_PREFIX"/02-GrabReads/"$SAMPLE_PREFIX"_trimmed_"$REPEAT"_dummied.fq >> ./"$SAMPLE_PREFIX"/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT_REVCOMP"_head.fastq &
			wait

			if [ -f "$SAMPLE_PREFIX"/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT"_tail.fastq ] &&
			   [ -f "$SAMPLE_PREFIX"/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT_REVCOMP"_head.fastq ]
			then 
				echo "[`date`] Read selection complete" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
			else 
				echo "[`date`] ERROR: Read selection failed. Check log for details" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
				mv "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt "$SAMPLE_PREFIX"_log_$(date +"%Y-%m-%d_%H:%M:%S").txt
				echo "[`date`] ERROR: Read selection failed. Check log for details"
				exit 1
			fi
		fi

		#3.Remove repeats from reads
		if [ -f "$SAMPLE_PREFIX"/03-RemoveRepeats/"$MIN_TAIL_LENGTH"_"$REPEAT"_tail.fastq ] &&
		   [ -f "$SAMPLE_PREFIX"/03-RemoveRepeats/"$MIN_TAIL_LENGTH"_"$REPEAT_REVCOMP"_head.fastq ]
		then 
			echo "[`date`] Repeat removal already complete" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
		else 
			echo "[`date`] Removing repeats from reads" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
			if [ ! -d "$SAMPLE_PREFIX"/03-RemoveRepeats ]; then mkdir -p "$SAMPLE_PREFIX"/03-RemoveRepeats; fi

			# remove repeats
			cutadapt -e 0 --minimum-length "$MIN_NONREPETITIVE_LENGTH" -a ${MAX_TAIL_1}$ -a ${MAX_TAIL_2}$ ./"$SAMPLE_PREFIX"/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT"_tail.fastq -o ./"$SAMPLE_PREFIX"/03-RemoveRepeats/${n}_"$REPEAT"_tail.fastq &
			cutadapt -e 0 --minimum-length "$MIN_NONREPETITIVE_LENGTH" -g ^${MAX_HEAD_1} -g ^${MAX_HEAD_2} ./"$SAMPLE_PREFIX"/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT_REVCOMP"_head.fastq -o ./"$SAMPLE_PREFIX"/03-RemoveRepeats/${n}_"$REPEAT_REVCOMP"_head.fastq &
			wait

			for i in $(seq 1 "$REMOVE_UP_TO")
			do
				cutadapt -e 0 --minimum-length "$MIN_NONREPETITIVE_LENGTH" -a ${MAX_TAIL_1:i}$ -a ${MAX_TAIL_2:i}$ ./"$SAMPLE_PREFIX"/03-RemoveRepeats/${n}_"$REPEAT"_tail.fastq -o ./"$SAMPLE_PREFIX"/03-RemoveRepeats/${m}_"$REPEAT"_tail.fastq &
				cutadapt -e 0 --minimum-length "$MIN_NONREPETITIVE_LENGTH" -g ^${MAX_HEAD_1:i} -g ^${MAX_HEAD_2:i} ./"$SAMPLE_PREFIX"/03-RemoveRepeats/${n}_"$REPEAT_REVCOMP"_head.fastq -o ./"$SAMPLE_PREFIX"/03-RemoveRepeats/${m}_"$REPEAT_REVCOMP"_head.fastq &
				wait
				m=$[$m-1] 
				n=$[$n-1] 
			done

			# restore n and m variables
			let i=0	
			let n="$MAX_NUM_REPEATS"*${#REPEAT}
			let m=${n}-1

			# remove temp files
			let start="$MIN_TAIL_LENGTH"+1
			let end=${MAX_NUM_REPEATS}*${#REPEAT}
			for i in $(seq "$start" "$end")
			do
				rm ./"$SAMPLE_PREFIX"/03-RemoveRepeats/${i}_"$REPEAT"_tail.fastq &
				rm ./"$SAMPLE_PREFIX"/03-RemoveRepeats/${i}_"$REPEAT_REVCOMP"_head.fastq &
				wait
			done

			if [ -f "$SAMPLE_PREFIX"/03-RemoveRepeats/"$MIN_TAIL_LENGTH"_"$REPEAT"_tail.fastq ] &&
			   [ -f "$SAMPLE_PREFIX"/03-RemoveRepeats/"$MIN_TAIL_LENGTH"_"$REPEAT_REVCOMP"_head.fastq ]
			then 
				echo "[`date`] Repeat removal complete" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
			else 
				echo "[`date`] ERROR: Repeat removal failed. Check log for details" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
				mv "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt "$SAMPLE_PREFIX"_log_$(date +"%Y-%m-%d_%H:%M:%S").txt
				echo "[`date`] ERROR: Repeat removal failed. Check log for details"
				exit 1
			fi
		fi

		#4.Align reads with repeats removed to genome
		if [ -f "$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT"-tail/accepted_hits.bam ] &&
		   [ -f "$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT_REVCOMP"-head/accepted_hits.bam ]
		then 
			echo "[`date`] Read alignment to "$REF_SEQ_PREFIX" already complete" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
		else 
			echo "[`date`] Aligning reads with repeats removed to "$REF_SEQ_PREFIX" genome" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
			if [ ! -d "$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT"-tail ]; then mkdir -p "$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT"-tail; fi
			if [ ! -d "$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT_REVCOMP"-head ]; then mkdir -p "$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT_REVCOMP"-head; fi

			# read alignment
			if [ $NO_NOVEL_JUNCS = "false" ]
			then
				tophat -p "$MAX_CORES" --read-mismatches 0 --read-gap-length 0 --max-deletion-length 0 --max-insertion-length 0 --read-edit-dist 0 --splice-mismatches 0 --segment-mismatches 0 -o ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT"-tail "$REF_FOR_TOPHAT" ./"$SAMPLE_PREFIX"/03-RemoveRepeats/"$MIN_TAIL_LENGTH"_"$REPEAT"_tail.fastq
				tophat -p "$MAX_CORES" --read-mismatches 0 --read-gap-length 0 --max-deletion-length 0 --max-insertion-length 0 --read-edit-dist 0 --splice-mismatches 0 --segment-mismatches 0 -o ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT_REVCOMP"-head "$REF_FOR_TOPHAT" ./"$SAMPLE_PREFIX"/03-RemoveRepeats/"$MIN_TAIL_LENGTH"_"$REPEAT_REVCOMP"_head.fastq
			else
				if [ $NO_NOVEL_JUNCS = "true" ]
				then
					tophat -p "$MAX_CORES" --no-novel-juncs --read-mismatches 0 --read-gap-length 0 --max-deletion-length 0 --max-insertion-length 0 --read-edit-dist 0 --splice-mismatches 0 --segment-mismatches 0 -o ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT"-tail "$REF_FOR_TOPHAT" ./"$SAMPLE_PREFIX"/03-RemoveRepeats/"$MIN_TAIL_LENGTH"_"$REPEAT"_tail.fastq
					tophat -p "$MAX_CORES" --no-novel-juncs --read-mismatches 0 --read-gap-length 0 --max-deletion-length 0 --max-insertion-length 0 --read-edit-dist 0 --splice-mismatches 0 --segment-mismatches 0 -o ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT_REVCOMP"-head "$REF_FOR_TOPHAT" ./"$SAMPLE_PREFIX"/03-RemoveRepeats/"$MIN_TAIL_LENGTH"_"$REPEAT_REVCOMP"_head.fastq
				fi
			fi

			if [ -f "$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT"-tail/accepted_hits.bam ] &&
			   [ -f "$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT_REVCOMP"-head/accepted_hits.bam ]
			then 
				echo "[`date`] Read alignment to "$REF_SEQ_PREFIX" complete" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
			else 
				echo "[`date`] ERROR: Read alignment to "$REF_SEQ_PREFIX" failed. Check log for details" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
				mv "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt "$SAMPLE_PREFIX"_log_$(date +"%Y-%m-%d_%H:%M:%S").txt
				echo "[`date`] ERROR: Read alignment to "$REF_SEQ_PREFIX" failed. Check log for details"
				exit 1
			fi
		fi

		#5.Add back repeats to aligned reads
		if [ -f "$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped_reads.fastq ] &&
		   [ -f "$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped_reads.fastq ]
		then 
			echo "[`date`] Adding back repeats already complete" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
		else 
			echo "[`date`] Adding back repeats to aligned reads" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
			if [ ! -d "$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/05-AddBackRepeats ]; then mkdir -p "$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/05-AddBackRepeats; fi

			# covert bam to sam
			samtools view ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT"-tail/accepted_hits.bam -o ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped.sam &
			samtools view ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT_REVCOMP"-head/accepted_hits.bam -o ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped.sam &
			wait

			# remove alignments to dummy chr that aren't dummy read
			awk -v rep="$REPEAT" 'BEGIN {FS="\t"; OFS="\t"} ($1 == "dummy"rep".1" && $3 == "dummy_chr_"rep"") || ($1 != "dummy"rep".1" && $3 != "dummy_chr_"rep"") {print $0}' ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped.sam > ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped_cleaned.sam &
			awk -v rep_rev="$REPEAT_REVCOMP" 'BEGIN {FS="\t"; OFS="\t"} ($1 == "dummy"rep_rev".1" && $3 == "dummy_chr_"rep"") || ($1 != "dummy"rep_rev".1" && $3 != "dummy_chr_"rep"") {print $0}' ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped.sam > ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped_cleaned.sam &
			wait

			# grab read names and remove duplicates
			awk '{print $1}' ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped_cleaned.sam | awk '!a[$0]++' > ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped_read_names &
			awk '{print $1}' ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped_cleaned.sam | awk '!a[$0]++' > ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped_read_names &
			wait

			# print fastq for given reads
			"$PATH_TO_SEQTK"/seqtk subseq ./"$SAMPLE_PREFIX"/02-GrabReads/"$SAMPLE_PREFIX"_trimmed_"$REPEAT"_dummied.fq ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped_read_names > ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped_reads.fastq &
			"$PATH_TO_SEQTK"/seqtk subseq ./"$SAMPLE_PREFIX"/02-GrabReads/"$SAMPLE_PREFIX"_trimmed_"$REPEAT"_dummied.fq ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped_read_names > ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped_reads.fastq &
			wait

			if [ -f "$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped_reads.fastq ] &&
			   [ -f "$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped_reads.fastq ]
			then 
				echo "[`date`] Adding back repeats complete" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
			else 
				echo "[`date`] ERROR: Adding back repeats failed. Check log for details" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
				mv "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt "$SAMPLE_PREFIX"_log_$(date +"%Y-%m-%d_%H:%M:%S").txt
				echo "[`date`] ERROR: Adding back repeats failed. Check log for details"
				exit 1
			fi
		fi

		#6.Align reads with repeats added back to genome
		if [ -f "$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT"-tail/accepted_hits.bam ] &&
		   [ -f "$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT_REVCOMP"-head/accepted_hits.bam ]
		then 
			echo "[`date`] Read alignment to "$REF_SEQ_PREFIX" already complete" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
		else 
			echo "[`date`] Aligning reads with repeats added back to "$REF_SEQ_PREFIX" genome" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
			if [ ! -d "$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT"-tail ]; then mkdir -p "$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT"-tail; fi
			if [ ! -d "$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT_REVCOMP"-head ]; then mkdir -p "$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT_REVCOMP"-head; fi

			# read alignment
			if [ $NO_NOVEL_JUNCS = "false" ]
			then
				tophat -p "$MAX_CORES" --read-mismatches 1 --read-gap-length 1 --max-deletion-length 1 --max-insertion-length 1 --read-edit-dist 1 --splice-mismatches 1 --segment-mismatches 1 -o ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT"-tail $REF_FOR_TOPHAT ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped_reads.fastq
				tophat -p "$MAX_CORES" --read-mismatches 1 --read-gap-length 1 --max-deletion-length 1 --max-insertion-length 1 --read-edit-dist 1 --splice-mismatches 1 --segment-mismatches 1 -o ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT_REVCOMP"-head $REF_FOR_TOPHAT ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped_reads.fastq
			else
				if [ $NO_NOVEL_JUNCS = "true" ]
				then
					tophat -p "$MAX_CORES" --no-novel-juncs --read-mismatches 1 --read-gap-length 1 --max-deletion-length 1 --max-insertion-length 1 --read-edit-dist 1 --splice-mismatches 1 --segment-mismatches 1 -o ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT"-tail $REF_FOR_TOPHAT ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped_reads.fastq
					tophat -p "$MAX_CORES" --no-novel-juncs --read-mismatches 1 --read-gap-length 1 --max-deletion-length 1 --max-insertion-length 1 --read-edit-dist 1 --splice-mismatches 1 --segment-mismatches 1 -o ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT_REVCOMP"-head $REF_FOR_TOPHAT ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped_reads.fastq
				fi
			fi

			if [ -f "$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT"-tail/accepted_hits.bam ] &&
			   [ -f "$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT_REVCOMP"-head/accepted_hits.bam ]
			then 
				echo "[`date`] Read alignment to "$REF_SEQ_PREFIX" complete" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
			else 
				echo "[`date`] ERROR: Read alignment to "$REF_SEQ_PREFIX" failed. Check log for details" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
				mv "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt "$SAMPLE_PREFIX"_log_$(date +"%Y-%m-%d_%H:%M:%S").txt
				echo "[`date`] ERROR: Read alignment to "$REF_SEQ_PREFIX" failed. Check log for details"
				exit 1
			fi
		fi

		#7.Print reads that failed to map and get mates
		if [ -f ./results/fastq/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates_read_names ] &&
		   [ -f ./results/fastq/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates.fastq ]
		then 
			echo "[`date`] Printing nontemplated p"$REPEAT" tail reads already complete" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
		else 
			echo "[`date`] Printing nontemplated p"$REPEAT" tail reads and their mates" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
			if [ ! -d ./results/fastq ]; then mkdir -p ./results/fastq; fi
			if [ ! -d "$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/07-PrintReads ]; then mkdir -p "$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/07-PrintReads; fi

			# covert bam to sam
			samtools view ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT"-tail/unmapped.bam -o ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/07-PrintReads/"$REPEAT"_tail_unmapped.sam &
			samtools view ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT_REVCOMP"-head/unmapped.bam -o ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/07-PrintReads/"$REPEAT_REVCOMP"_head_unmapped.sam &
			wait

			# grab read names and remove duplicates
			awk '{print $1}' ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/07-PrintReads/"$REPEAT"_tail_unmapped.sam | awk '!a[$0]++' > ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/07-PrintReads/"$REPEAT"_tail_unmapped_read_names &
			awk '{print $1}' ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/07-PrintReads/"$REPEAT_REVCOMP"_head_unmapped.sam | awk '!a[$0]++' > ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/07-PrintReads/"$REPEAT_REVCOMP"_head_unmapped_read_names &
			wait

			# combine read names and remove duplicates
			cat ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/07-PrintReads/"$REPEAT"_tail_unmapped_read_names ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/07-PrintReads/"$REPEAT_REVCOMP"_head_unmapped_read_names | sort | uniq > ./results/fastq/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates_read_names

			# if read name files are empty (i.e. no candidates were identified), then exit script.
			if [[ $(wc ./results/fastq/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates_read_names | awk -F " " '{print $1}') = 0 ]]
			then
				# STATISTICS
				NUM_READS_BEFORE_TRIMMING=$(grep "Total reads processed" ./"$SAMPLE_PREFIX"/01-TrimReads/"$SAMPLE_PREFIX".fastq.gz_trimming_report.txt | awk '{print $4}' | sed 's/,//g')

				NUM_READS_AFTER_TRIMMING=$(grep "Number of input reads" ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/star/"$SAMPLE_PREFIX"_Log.final.out | awk -F "\t" '{print $NF}')

				PERCENT_REMAINING_AFTER_TRIM=$(awk -v TOTAL="$NUM_READS_BEFORE_TRIMMING" -v TRIM="$NUM_READS_AFTER_TRIMMING" 'BEGIN{printf("%.2f\n",(TRIM/TOTAL)*100)}')

				NUM_READS_THAT_UNIQUELY_MAP_TO_GENOME=$(grep "Uniquely mapped reads number" ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/star/"$SAMPLE_PREFIX"_Log.final.out | awk -F "\t" '{print $NF}')

				PERCENT_OF_TRIMMED_READS_THAT_MAP=$(awk -v TRIM="$NUM_READS_AFTER_TRIMMING" -v MAP="$NUM_READS_THAT_UNIQUELY_MAP_TO_GENOME" 'BEGIN{printf("%.2f\n",(MAP/TRIM)*100)}')

				NUM_READS_WITH_REPEAT_TAILS=$(cat \
				    ./"$SAMPLE_PREFIX"/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT"_tail.fastq \
				    ./"$SAMPLE_PREFIX"/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT_REVCOMP"_head.fastq \
				    | awk 'NR % 4 == 1' | grep -v dummy | sort | uniq | wc -l)

				NUM_READS_WITH_REPEAT_TAILS_THAT_MAP_TO_GENOME=$(cat \
					./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped_read_names \
					./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped_read_names \
					| grep -v dummy | wc -l)
					#read names from SAM files do not have read 1 or read 2 identifiers, so don't uniqify

				PERCENT_OF_REPEAT_READS_THAT_MAP=$(awk -v REPEATS="$NUM_READS_WITH_REPEAT_TAILS" -v MAP="$NUM_READS_WITH_REPEAT_TAILS_THAT_MAP_TO_GENOME" 'BEGIN{printf("%.2f\n",(MAP/REPEATS)*100)}')

				NUM_MAPPED_READS_NONTEMPLATED=$(wc -l ./results/fastq/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates_read_names | awk '{print $1}')

				PERCENT_OF_MAPPED_READS_THAT_ARE_NONTEMPLATED=$(awk -v MAP="$NUM_READS_WITH_REPEAT_TAILS_THAT_MAP_TO_GENOME" -v NONTEMP="$NUM_MAPPED_READS_NONTEMPLATED" 'BEGIN{printf("%.2f\n",(NONTEMP/MAP)*100)}')

				STATISTICS="=== STATISTICS ===

Total reads processed: 								$(printf "%'.f\n" $NUM_READS_BEFORE_TRIMMING)
Reads remaining after adapter and quality trim:		$(printf "%'.f\n" $NUM_READS_AFTER_TRIMMING) ("$PERCENT_REMAINING_AFTER_TRIM"% of total reads)
Reads that map to "$REF_SEQ_PREFIX" genome:				$(printf "%'.f\n" $NUM_READS_THAT_UNIQUELY_MAP_TO_GENOME) ("$PERCENT_OF_TRIMMED_READS_THAT_MAP"% of trimmed reads)

p"$REPEAT" reads with at least "$MIN_NUM_REPEATS" repeats: 				$(printf "%'.f\n" $NUM_READS_WITH_REPEAT_TAILS)
p"$REPEAT" reads that map to "$REF_SEQ_PREFIX" genome:			$(printf "%'.f\n" $NUM_READS_WITH_REPEAT_TAILS_THAT_MAP_TO_GENOME) ("$PERCENT_OF_REPEAT_READS_THAT_MAP"% of p"$REPEAT" reads)
p"$REPEAT" reads that are nontemplated: 					$(printf "%'.f\n" $NUM_MAPPED_READS_NONTEMPLATED) ("$PERCENT_OF_MAPPED_READS_THAT_ARE_NONTEMPLATED"% of p"$REPEAT" mapping reads)
"
				
				echo "[`date`] Failed to find any p"$REPEAT" reads that are nontemplated" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
				echo "[`date`] FIND NONTEMPLATED "$REPEAT" TAILS COMPLETE" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
				printf -- "%s\n" " " >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt 
				echo "$STATISTICS" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt

				mv "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log_final.txt
				echo "[`date`] FIND NONTEMPLATED "$REPEAT" TAILS COMPLETE. Check "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log_final.txt for summary"
				exit 0
			fi

			# grab reads
			echo | tr -d '\n' > ./results/fastq/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates.fastq

			for i in $(cat ./results/fastq/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates_read_names)
			do
				grep -A 3 @"$i" ./"$SAMPLE_PREFIX"/02-GrabReads/"$SAMPLE_PREFIX"_trimmed_"$REPEAT"_dummied.fq >> ./results/fastq/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates.fastq
			done

			# copy 
			if [ ! -d ./../candidates ]; then mkdir -p ./../candidates; fi
			cp ./results/fastq/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates.fastq ./../candidates

			if [ -f ./results/fastq/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates_read_names ] &&
			   [ -f ./results/fastq/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates.fastq ]
			then 
				echo "[`date`] Printing nontemplated p"$REPEAT" tail reads complete" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
			else 
				echo "[`date`] ERROR: Printing nontemplated p"$REPEAT" tail reads failed. Check log for details" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
				mv "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt "$SAMPLE_PREFIX"_log_$(date +"%Y-%m-%d_%H:%M:%S").txt
				echo "[`date`] ERROR: Printing nontemplated p"$REPEAT" tail reads failed. Check log for details"
				exit 1
			fi
		fi
	
		#8.Align candidates to genome
		if [ -f ./results/alignments/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates.bam ] &&
		   [ -f ./results/alignments/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates.bam.bai ]
		then 
			echo "[`date`] p"$REPEAT" tail candidate read alignment to "$REF_SEQ_PREFIX" genome already complete" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
		else 
			echo "[`date`] Aligning nontemplated p"$REPEAT" tail candidates to "$REF_SEQ_PREFIX" genome" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
			if [ ! -d ./results/alignments ]; then mkdir -p ./results/alignments; fi
			if [ ! -d "$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/08-AlignCandidates ]; then mkdir -p "$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/08-AlignCandidates; fi

			# covert bam to sam
			samtools view ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT"-tail/accepted_hits.bam -o ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/08-AlignCandidates/"$REPEAT"_tail_mapped.sam &
			samtools view ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/06-AlignReadsWithRepeatsAddedBack-"$REPEAT_REVCOMP"-head/accepted_hits.bam -o ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/08-AlignCandidates/"$REPEAT_REVCOMP"_head_mapped.sam &
			wait

			# grab read names and remove duplicates
			awk '{print $1}' ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/08-AlignCandidates/"$REPEAT"_tail_mapped.sam | awk '!a[$0]++' > ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/08-AlignCandidates/"$REPEAT"_tail_mapped_read_names &
			awk '{print $1}' ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/08-AlignCandidates/"$REPEAT_REVCOMP"_head_mapped.sam | awk '!a[$0]++' > ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/08-AlignCandidates/"$REPEAT_REVCOMP"_head_mapped_read_names &
			wait

			# remove reads that aligned with repeats added back from the first alignment with reads with repeats removed
			for i in $(cat ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/08-AlignCandidates/"$REPEAT"_tail_mapped_read_names); do grep -v "$i" ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped_cleaned.sam > ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/08-AlignCandidates/"$REPEAT"_tail_mapped_cleaned_geno_encode_removed.sam; done &
			for i in $(cat ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/08-AlignCandidates/"$REPEAT_REVCOMP"_head_mapped_read_names); do grep -v "$i" ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped_cleaned.sam > ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/08-AlignCandidates/"$REPEAT_REVCOMP"_head_mapped_cleaned_geno_encode_removed.sam; done &
			wait

			# combine alignments and add header
			samtools view -H ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/04-AlignReadsWithRepeatsRemoved-"$REPEAT"-tail/accepted_hits.bam | cat - ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/08-AlignCandidates/"$REPEAT"_tail_mapped_cleaned_geno_encode_removed.sam ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/08-AlignCandidates/"$REPEAT_REVCOMP"_head_mapped_cleaned_geno_encode_removed.sam > ./results/alignments/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates.sam

			# sort sam and output as bam
			samtools sort ./results/alignments/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates.sam -O BAM -o ./results/alignments/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates.bam
			
#			if [ $FEATURE_COUNTS = "true" ]
#			then
#				"$PATH_TO_SUBREAD"/featureCounts \
#				-a "$ANNOTATION" \
#				-o ./results/alignments/"$SAMPLE_PREFIX"_p"$REPEAT"_counts.txt \
#				-O \
#				./results/alignments/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates.bam
#			fi

			# index bam
			samtools index ./results/alignments/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates.bam

			# copy
			if [ ! -d ./../candidates ]; then mkdir -p ./../candidates; fi
			cp ./results/alignments/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates.bam ./../candidates &
			cp ./results/alignments/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates.bam.bai ./../candidates &
			wait

			if [ -f ./results/alignments/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates.bam ] &&
			   [ -f ./results/alignments/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates.bam.bai ]
			then 
				echo "[`date`] Aligning nontemplated p"$REPEAT" tail candidates to "$REF_SEQ_PREFIX" genome complete" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
			else 
				echo "[`date`] ERROR: Aligning nontemplated p"$REPEAT" tail candidates to "$REF_SEQ_PREFIX" genome failed. Check log for details" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
				mv "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt "$SAMPLE_PREFIX"_log_$(date +"%Y-%m-%d_%H:%M:%S").txt
				echo "[`date`] ERROR: Aligning nontemplated p"$REPEAT" tail candidates to "$REF_SEQ_PREFIX" genome failed. Check log for details"
				exit 1
			fi
		fi
	
		# STATISTICS
		NUM_READS_BEFORE_TRIMMING=$(grep "Total reads processed" ./"$SAMPLE_PREFIX"/01-TrimReads/"$SAMPLE_PREFIX".fastq.gz_trimming_report.txt | awk '{print $4}' | sed 's/,//g')

		NUM_READS_AFTER_TRIMMING=$(grep "Number of input reads" ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/star/"$SAMPLE_PREFIX"_Log.final.out | awk -F "\t" '{print $NF}')

		PERCENT_REMAINING_AFTER_TRIM=$(awk -v TOTAL="$NUM_READS_BEFORE_TRIMMING" -v TRIM="$NUM_READS_AFTER_TRIMMING" 'BEGIN{printf("%.2f\n",(TRIM/TOTAL)*100)}')

		NUM_READS_THAT_UNIQUELY_MAP_TO_GENOME=$(grep "Uniquely mapped reads number" ./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/star/"$SAMPLE_PREFIX"_Log.final.out | awk -F "\t" '{print $NF}')

		PERCENT_OF_TRIMMED_READS_THAT_MAP=$(awk -v TRIM="$NUM_READS_AFTER_TRIMMING" -v MAP="$NUM_READS_THAT_UNIQUELY_MAP_TO_GENOME" 'BEGIN{printf("%.2f\n",(MAP/TRIM)*100)}')

		NUM_READS_WITH_REPEAT_TAILS=$(cat \
		    ./"$SAMPLE_PREFIX"/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT"_tail.fastq \
		    ./"$SAMPLE_PREFIX"/02-GrabReads/"$MIN_NUM_REPEATS"_"$REPEAT_REVCOMP"_head.fastq \
		    | awk 'NR % 4 == 1' | grep -v dummy | sort | uniq | wc -l)

		NUM_READS_WITH_REPEAT_TAILS_THAT_MAP_TO_GENOME=$(cat \
			./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT"_tail_mapped_read_names \
			./"$SAMPLE_PREFIX"/"$REF_SEQ_PREFIX"/05-AddBackRepeats/"$REPEAT_REVCOMP"_head_mapped_read_names \
			| grep -v dummy | wc -l)
			#read names from SAM files do not have read 1 or read 2 identifiers, so don't uniqify

		PERCENT_OF_REPEAT_READS_THAT_MAP=$(awk -v REPEATS="$NUM_READS_WITH_REPEAT_TAILS" -v MAP="$NUM_READS_WITH_REPEAT_TAILS_THAT_MAP_TO_GENOME" 'BEGIN{printf("%.2f\n",(MAP/REPEATS)*100)}')

		NUM_MAPPED_READS_NONTEMPLATED=$(wc -l ./results/fastq/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates_read_names | awk '{print $1}')

		PERCENT_OF_MAPPED_READS_THAT_ARE_NONTEMPLATED=$(awk -v MAP="$NUM_READS_WITH_REPEAT_TAILS_THAT_MAP_TO_GENOME" -v NONTEMP="$NUM_MAPPED_READS_NONTEMPLATED" 'BEGIN{printf("%.2f\n",(NONTEMP/MAP)*100)}')

		STATISTICS="=== STATISTICS ===

Total reads processed: 								$(printf "%'.f\n" $NUM_READS_BEFORE_TRIMMING)
Reads remaining after adapter and quality trim:		$(printf "%'.f\n" $NUM_READS_AFTER_TRIMMING) ("$PERCENT_REMAINING_AFTER_TRIM"% of total reads)
Reads that map to "$REF_SEQ_PREFIX" genome:				$(printf "%'.f\n" $NUM_READS_THAT_UNIQUELY_MAP_TO_GENOME) ("$PERCENT_OF_TRIMMED_READS_THAT_MAP"% of trimmed reads)

p"$REPEAT" reads with at least "$MIN_NUM_REPEATS" repeats: 				$(printf "%'.f\n" $NUM_READS_WITH_REPEAT_TAILS)
p"$REPEAT" reads that map to "$REF_SEQ_PREFIX" genome:			$(printf "%'.f\n" $NUM_READS_WITH_REPEAT_TAILS_THAT_MAP_TO_GENOME) ("$PERCENT_OF_REPEAT_READS_THAT_MAP"% of p"$REPEAT" reads)
p"$REPEAT" reads that are nontemplated: 					$(printf "%'.f\n" $NUM_MAPPED_READS_NONTEMPLATED) ("$PERCENT_OF_MAPPED_READS_THAT_ARE_NONTEMPLATED"% of p"$REPEAT" mapping reads)
"

		rm ./results/alignments/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates.sam
		rm ./results/fastq/"$SAMPLE_PREFIX"_p"$REPEAT"_tail_candidates_read_names
		echo "[`date`] FIND NONTEMPLATED "$REPEAT" TAILS COMPLETE" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt
		printf -- "%s\n" " " >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt 
		echo "$STATISTICS" >> "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt

		mv "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log.txt "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log_final.txt
		echo "[`date`] FIND NONTEMPLATED "$REPEAT" TAILS COMPLETE. Check "$SAMPLE_PREFIX"_p"$REPEAT"-"$MIN_NUM_REPEATS"x_log_final.txt for summary"
		exit 0
	fi

STATUS=$?

exit $STATUS
