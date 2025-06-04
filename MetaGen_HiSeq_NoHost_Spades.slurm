#!/usr/bin/env bash

# - View information about jobs using:  squeue
# Find more instructions at https://slurm.schedmd.com/quickstart.html
# And at http://www.brightcomputing.com/blog/bid/174099/slurm-101-basic-slurm-usage-for-linux-clusters
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --time=96:00:00
#SBATCH --job-name=DuskMatter
#SBATCH --output=DuskMatter_%j.out

## Metagenomics pipeline for Illumina paired-end reads ##
## RNA-Viral Spades and Diamond to Viral Discovery in SRA Public Data ###

echo ; echo "### STARTING AT ###"
date

### Name of the compute node
HOST=$(hostname)
echo $HOST

### use $1, $2, $3, ... to receive arguments via command line. For this script it will be necessary:

### $1 = COMPLETE PATH FOR LIBRARY DIR
### $2 = COMPLETE PATH FOR KRONA TAXONOMY FOLDER
### $3 = COMPLETE PATH FOR DIAMOND DB
### $4 = PATTERN OF LIBRARY DIRECTORIES
### $5 = MEMORY
### $6 = NAME OF RESULTS DIRECTORY

###EXAMPLE:

### date; sbatch /data04/projects04/MarceloSoares/viromapn/bin/Scripts/Others/MetaGen_HiSeq_NoHost_Spades.slurm /data04/projects04/MarceloSoares/viromapn/data/Henrique/MOREY2017/Henrique_DuskMatter_215987/ /data04/projects04/MarceloSoares/viromapn/bin/KronaTools-2.8.1/taxonomy/ /data03/databases/sources/ncbi/blast/dbV5/nr SRR 90 DuskMatter_RNAspades_MOREY2017

LIBDIR=$1
KRONA=$2
DIAMOND_DB=$3
PATTERN=$4
MEMO=$5
RESULTS=$6

OUTDIR=/data04/projects04/MarceloSoares/viromapn/results/DuskMatter/${RESULTS}_${SLURM_JOB_ID}

### Define a name for the temporary directory in the scratch disc of the compute node - MUST be in /scr  !!!

TMPDIR=/scr/${USER}-${SLURM_JOB_ID}

### threads

NTs=$SLURM_CPUS_PER_TASK


################################################################################################

### functions (Nicole)

checkfile() { if [[ -f $1 ]] ; then echo "$1 OK!";  else echo "- File $1 not found!"; exit; fi; }
checkdir() { if [[ -d $1 ]] ; then echo "$1 OK!";  else echo "- Directory $1 not found!"; exit; fi; }
chkmkdir() { if [[ -e $1 ]] ; then echo "$1 already exists.";  else  mkdir -p $1 ; checkdir $1 ; fi; }
checkNonEmpty()  { if ! [[ -s $1 ]] ; then echo "- File $1 not found or zero sized!"; exit; fi; }

################################################################################################

echo ; echo "### Check scratch directory"

ls -rtl /scr

### Create temporary directories

chkmkdir $TMPDIR
chkmkdir $TMPDIR/DIAMOND_OUTPUT_CONTIGS
chkmkdir $TMPDIR/ASSEMBLIES
chkmkdir $TMPDIR/FILTERED_DATA

echo ; echo "### Temporary files will be writen in $TMPDIR execution host $HOST"

### Create output directories

chkmkdir $OUTDIR
chkmkdir $OUTDIR/DIAMOND_OUTPUT_CONTIGS
chkmkdir $OUTDIR/ASSEMBLIES
chkmkdir $OUTDIR/FILTERED_DATA

echo ; echo "### Final files will be writen in $OUTDIR"


###--------------------------------------------------------------------------------------###

###--------------------------------------------------------------------------------------###

### store your current working directory in a variable

OLDWD=$(pwd -P)


###--------------------------------------------------------------------------------------###

###--------------------------------------------------------------------------------------###

### Go to your scratch directory before start the script

echo ; echo "### Change your current working directory to $TMPDIR"

cd $TMPDIR/

### check current working directory

echo
pwd
ls -rtl


###--------------------------------------------------------------------------------------###

###--------------------------------------------------------------------------------------###

### Activate LOOP modules

# Activate SAMTOOLS
module load samtools/1.16.1

# Activate SPADES module
module load spades/3.15.5

# Activate DIAMOND module
module load diamond/2.1.9


###--------------------------------------------------------------------------------------###

###--------------------------------------------------------------------------------------###

### LOOP

# For each sample

for fullpath in $LIBDIR/${PATTERN}*; do

    # Extract just the directory name (without full path)

    i=$(basename "$fullpath")
    
    # Copy directory to tempdir

    echo ; echo -e "\nCopying $i to $TMPDIR..."

    scp -r "$fullpath" "$TMPDIR/"
    
    # Verify it's a valid directory

    if [ -d "$TMPDIR/$i" ]; then

	echo -e "\n#############################"
        echo "Processing directory: ${i}"
        echo "#############################"
        echo ; date ; echo

        # Enter directory for processing

	cd "$TMPDIR/$i" || { echo "Error entering directory"; exit 1; }

	# Get name of raw data 

        R1="$(ls *.fastq* | head -1)"
        echo "Paired-end reads file 1: $R1"
        R2="$(ls *.fastq* | tail -1)"
        echo "Paired-end reads file 2: $R2"
        NAME=$(echo $R1 | sed -E 's/_1.+//g')

	
	echo ; echo "########## STEP 01: QUALITY CONTROL (FASTP) ##########" ; date

	# Filter data with fastp

	echo ; echo "Removing sequences with: <50bp Lenght, <30 QScore and redundancy (duplication)"

	echo ; echo "/data04/projects04/MarceloSoares/viromapn/bin/fastp --dedup --detect_adapter_for_pe -l 50 -q 30 --thread $NTs -i $R1 -I $R2 -o $NAME.trim_p_R1.fastq -O $NAME.trim_p_R2.fastq"

	/data04/projects04/MarceloSoares/viromapn/bin/fastp --dedup --detect_adapter_for_pe -l 50 -q 30 --thread $NTs -i $R1 -I $R2 -o $NAME.trim_p_R1.fastq -O $NAME.trim_p_R2.fastq
	
        echo ; echo "Checking FILTER files in the directory..."
        pwd 
        ls -rtl

        # Move filtered data to output diretory                                                                        

        rsync -tvr $NAME.trim_p_R1.fastq $TMPDIR/FILTERED_DATA/
        rsync -tvr $NAME.trim_p_R2.fastq $TMPDIR/FILTERED_DATA/  
	

        echo; echo "########## STEP 02: DE NOVO ASSEMBLY (SPADES AND SAMTOOLS) ##########" ; date

        # De novo assembly (contigs)

        echo ; echo "Performing 'de novo' assembling (contigs) using Spades for library $NAME"
	
	echo ; echo "spades.py --rnaviral -1 $NAME.trim_p_R1.fastq -2 $NAME.trim_p_R2.fastq -t $NTs -m $MEMO --only-assemble -o $NAME.SPADES_ASSEMBLY"

	spades.py --rnaviral -1 $NAME.trim_p_R1.fastq -2 $NAME.trim_p_R2.fastq -t $NTs -m $MEMO --only-assemble -o $NAME.SPADES_ASSEMBLY
	
        # Go to assembly directory

        cd $NAME.SPADES_ASSEMBLY/
         
        #May add an if here

        if [ -f contigs.fasta ]; then

        mv contigs.fasta $NAME.contigs.fasta

        fi

	### SAMTOOLS
	
	# Index the FASTA file using samtools faidx
	samtools faidx "$NAME.contigs.fasta"

	# Get all sequence names first

	grep "^>" "$NAME.contigs.fasta" | cut -d' ' -f1 | sed 's/>//' > contig_names.txt
	
	# Calculate splitting (more robust)

	total_seqs=$(wc -l < contig_names.txt)

	echo "Number of sequences: $total_seqs"

	# Split into N equal parts (10 files)

	N=10

	seqs_per_file=$(( ("$total_seqs" + N - 1) / N ))
	
	# Split names file

	split -l $seqs_per_file -d contig_names.txt contig_names_part_
	
	# Extract using samtools

	for p in {00..09}; do

	    samtools faidx "$NAME.contigs.fasta" $(cat "contig_names_part_$p") > "${NAME}_part${p#0}_contigs.fasta"

	done

        echo ; echo "Checking DE NOVO files in the directory..."
        pwd
        ls -rtl

	
        echo ; echo "########## STEP 03: TAXONOMIC ASSIGNMENT (DIAMOND AND KRONA) ##########" ; date

        # Use diamond to perform taxonomic assignment - contigs

        echo ; echo "Performing DIAMOND taxonomic assignment (contigs) of library $NAME -> output DIAMOND (.txt)"

	for split_file in *part*.fasta; do

	    # Safely handle cases where no files match the pattern

	    [ -e "$split_file" ] || continue
    
	    # Get base name without extension

	    base_name="${split_file%.fasta}"
    
	    echo -e "\nProcessing $split_file → ${base_name}_diamond.txt"

	    echo ; echo "diamond blastx --threads "$NTs" --query "$split_file" --db "$DIAMOND_DB" --out "${base_name}_diamond.txt" --outfmt 6 --max-target-seqs 1 --evalue 0.00001"
	    
	    diamond blastx --threads "$NTs" --query "$split_file" --db "$DIAMOND_DB" --out "${base_name}_diamond.txt" --outfmt 6 --max-target-seqs 1 --evalue 0.00001
    
	    # Verify output

	    if [[ ! -f "${base_name}_diamond.txt" ]]; then

		echo "ERROR: DIAMOND failed for $split_file" >&2

		exit 1

	    fi

	done

	echo ; echo "Checking DIAMOND files (contigs) in the directory..."
        pwd
        ls -rtl

	# Create unified file

	echo ; echo -e "\nCreating Unified DIAMOND file for $NAME"
	
	if ls *part*_diamond.txt 1> /dev/null 2>&1; then
	    
	    echo -e "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore" > "${NAME}_unified_diamond.txt"
	    
	    cat *part*_diamond.txt >> "${NAME}_unified_diamond.txt"
	    
	    echo "Results combined in ${NAME}_unified_diamond.txt"
	else
	    echo "No DIAMOND result files (*part*_diamond.txt) found to combine."
	fi
		    
        # Move only important files to output directories
	
        rsync -tvr ${NAME}_unified_diamond.txt $TMPDIR/DIAMOND_OUTPUT_CONTIGS/
	rsync -tvr ${NAME}.contigs.fasta $TMPDIR/ASSEMBLIES/

	# Go back to sample directory and delete SPADES directory

	cd ../                                                                
        rm -rf ${NAME}.SPADES_ASSEMBLY
	
	# Go to DIAMOND directory                                                                                          

        cd $TMPDIR/DIAMOND_OUTPUT_CONTIGS/
	
        echo ; echo "Creating DIAMOND diversity plot (contigs) for library $NAME -> output Krona (.html)"

        echo ; echo "perl /data04/projects04/MarceloSoares/viromapn/bin/KronaTools-2.8.1/scripts/ImportBLAST.pl "${NAME}_unified_diamond.txt" -tax "$KRONA" -o "${NAME}_unified_krona.html""
	
	perl /data04/projects04/MarceloSoares/viromapn/bin/KronaTools-2.8.1/scripts/ImportBLAST.pl "${NAME}_unified_diamond.txt" -tax "$KRONA" -o "${NAME}_unified_krona.html"

        echo ; echo "Checking DIAMOND/Krona files in the directory..."
	pwd
	ls -rtl
	
        # Go back to the basal directory

        cd $TMPDIR/

	# Move final directories to output directories                                                                                      
	mv DIAMOND_OUTPUT_CONTIGS/* ${OUTDIR}/DIAMOND_OUTPUT_CONTIGS/.
	mv ASSEMBLIES/* ${OUTDIR}/ASSEMBLIES/.
	mv FILTERED_DATA/* ${OUTDIR}/FILTERED_DATA/.
	
        # Remove processed sapmle to free memory

	echo ; echo "Delete $i directory from tempdir"

	echo ; echo "rm -rf $TMPDIR/$i"
        rm -rf $TMPDIR/$i

        echo ; echo "Final check in the tempdir before moving to next sample..."
	pwd
	ls -rtl

    else
        echo ; echo "Warning: $i is not a valid directory or wasn't copied correctly" ; date ; echo
    fi

    echo ; date ; echo
    
done


###--------------------------------------------------------------------------------------###

###--------------------------------------------------------------------------------------###

### Deactivate LOOP modules

# Deactivate SAMTOOLS module
module unload samtools/1.16.1

# Deactivate SPADES module
module unload spades/3.15.5

# Deactivating DIAMOND module 
module unload diamond/2.1.9


###--------------------------------------------------------------------------------------###

### check the size of your directories -- just to keep track

###--------------------------------------------------------------------------------------###

echo ; echo "### Check the size of $TMPDIR"
du -hs $TMPDIR/

echo ; echo "### Check the size of $OUTDIR"
du -hs $OUTDIR/


###--------------------------------------------------------------------------------------###

### clear the scratch for your colleagues

###--------------------------------------------------------------------------------------###

echo ; echo "## Remove temporary directory"

if   [[ -d $TMPDIR ]] ; then     

    echo ; echo "$TMPDIR found. Removing."; 
    /bin/rm -rf $TMPDIR/ ; 

else 

    echo ; echo "- Directory $TMPDIR not found!";   

fi 

### just check -- debug

echo ; echo "### Check again scratch directory"
ls -rtl /scr


###--------------------------------------------------------------------------------------### 

### get out of /scr

###--------------------------------------------------------------------------------------###

echo ; echo "### Navegate back to your initial working directory $OLDWD"
cd $OLDWD
pwd
echo ; date ; echo

###--------------------------------------------------------------------------------------###

echo "#######################"
echo "####### The end #######"
echo "#######################"

exit


