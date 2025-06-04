#!/usr/bin/env bash

# - View information about jobs using:  squeue
# Find more instructions at https://slurm.schedmd.com/quickstart.html
# And at http://www.brightcomputing.com/blog/bid/174099/slurm-101-basic-slurm-usage-for-linux-clusters
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=56
#SBATCH --time=96:00:00
#SBATCH --job-name=SRA_DuskMatter
#SBATCH --output=SRA_DuskMatter_%j.out
#SBATCH --mem=60G

 ##  pipeline for fastq Download to explore DuskMatter ##
 echo ; echo "### STARTING AT ###"
 date

### Name of the compute node
HOST=$(hostname)
echo $HOST


### use $1, $2, $3, ... to receive arguments via command line. For this script it will be necessary:

###EXAMPLE:

### date; sbatch ~/viromapn/bin/Scripts/download_SRA.slurm ~/viromapn/data/Henrique/Henrique_DuskMatter_accesions.txt Henrique_DuskMatter

### $1 = TXT FILE WITH ONE SRA ACESSION PER LINE OF INTEREST
### $2 = OUTPUT NAME OF INTEREST

# Variables
ACC=$1
RESULTS=$2

# Create variables of dirs of interest
OUTDIR=/data04/projects04/MarceloSoares/viromapn/data/SRA/${RESULTS}_${SLURM_JOB_ID}


# Activate sratoolkit module
module load sratoolkit/3.1.0


### functions (Nicole)
checkfile() { if [[ -f $1 ]] ; then echo "$1 OK!";  else echo "- File $1 not found!"; exit; fi; }
checkdir() { if [[ -d $1 ]] ; then echo "$1 OK!";  else echo "- Directory $1 not found!"; exit; fi; }
chkmkdir() { if [[ -e $1 ]] ; then echo "$1 already exists.";  else  mkdir -p $1 ; checkdir $1 ; fi; }
checkNonEmpty()  { if ! [[ -s $1 ]] ; then echo "- File $1 not found or zero sized!"; exit; fi; }
################################################################################################

### Create output dir
chkmkdir $OUTDIR
echo ; echo "### Final files will be writen in $OUTDIR"

###--------------------------------------------------------------------------------------###
### store your current working directory in a variable
OLDWD=$(pwd -P)
###--------------------------------------------------------------------------------------###
	while read -r p; 
		do
		[[ -z $p ]] && continue  # Skip empty lines
		echo "$p"
		chkmkdir $TMPDIR/"$p"
		fasterq-dump --outdir $OUTDIR/"$p" "$p" --verbose
		sleep 240
	done < "$ACC"
	
###--------------------------------------------------------------------------------------###

###--------------------------------------------------------------------------------------###
### Deactivate LOOP modules

# Deactivate sratoolkits module
module unload sratoolkit/3.1.0

###--------------------------------------------------------------------------------------###
### check the size of your temporary directory -- just to keep track
###--------------------------------------------------------------------------------------###

du -hs $TMPDIR/
echo ; echo "### Final results will be found in $OUTDIR." ; ls -rtl $OUTDIR/

###--------------------------------------------------------------------------------------###
### clear the scratch for your colleagues
###--------------------------------------------------------------------------------------###

echo ; echo "## Remove temporary directory"
if   [[ -d $TMPDIR ]] ; then 
    
    echo ; echo "$TMPDIR found. Removing."; 
    rm -rf $TMPDIR/ ; 

else 
    echo ; echo "- Directory $TMPDIR not found!";   
fi 

### just check -- debug
echo ; echo "### Check again scratch directory" ; ls -rtl /scr

###--------------------------------------------------------------------------------------### 
### get out of /scr
###--------------------------------------------------------------------------------------###

echo ; echo "### Navegate back to your initial working directory $OLDWD"
cd $OLDWD ; pwd ; echo ; date ; echo

###--------------------------------------------------------------------------------------###

echo "#######################"
echo "########The end########"
echo "#######################"

exit
