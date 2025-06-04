#!/usr/bin/env bash



# - View information about jobs using:  squeue

# Find more instructions at https://slurm.schedmd.com/quickstart.html

# And at http://www.brightcomputing.com/blog/bid/174099/slurm-101-basic-slurm-usage-for-linux-clusters

#SBATCH --ntasks=1

#SBATCH --cpus-per-task=56

#SBATCH --time=48:00:00

#SBATCH --job-name=RdRpScan

#SBATCH --output=RdRpScan_%j.out



## DuskMatter pipeline for Viral Identification on Diamond Unnasigned Contigs ##

## Developed by Cosentino & D'arc ####



echo ; echo "### STARTING AT ###"

date



### Name of the compute node

HOST=$(hostname)

echo $HOST





echo "##############################"

echo "##Step 0. Veryfing Arguments##"

echo "##############################"



if [[ -d $INPUT ]] 

then 

	echo "Root directory for samples indentified ->  $INPUT"

else 

	echo "Samples directory not found. Please specify."

	exit

fi



# directory for output

if [[ -d $RESULT ]] 

then 

	echo "Output directory indentified ->  $RESULT"

	echo "Please specify another to avoid overwriting previous analysis."

	exit

else 

	echo "Samples directory created -> $RESULT"

fi



# directory for rdpr_scan database

if [[ -d $RDRP ]] 

then 

	echo "Rdpr_scan directory indentified ->  $RDRP"

else 

	echo "Rdpr_scan directory not indentified -> $RDRP"

	exit

fi



# directory for DIAMOND files previously done 

if [[ -d $DIAMOND ]] 

then 

	echo "Diamond directory indentified ->  $DIAMOND"

else 

	echo "Diamond directory not indentified  -> $DIAMOND"

	exit

fi





###### Define a name for the Output Directory - MUST be in RdRpScan

OUTPUT=/data04/projects04/MarceloSoares/viromapn/results/RdRpScan/${RESULT}_${SLURM_JOB_ID}



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



### CREATE this temporary directory (if not already available):

chkmkdir $TMPDIR

echo ; echo "### Temporary files will be writen in $TMPDIR execution host $HOST"



### Create output dir

chkmkdir $OUTPUT

echo ; echo "### Final files will be writen in $OUTPUT"



###--------------------------------------------------------------------------------------###

### store your current working directory in a variable

OLDWD=$(pwd -P)

###--------------------------------------------------------------------------------------###



###--------------------------------------------------------------------------------------###

### you can / should copy your input files to your scratch directory and save this location for the rest of the script



scp -r $INPUT/* $TMPDIR/

scp -r $DIAMOND/* $TMPDIR/



echo ; echo "### Change your current working directory to $TMPDIR"

cd $TMPDIR/





### check current working directory

echo ; pwd ; ls -rtl



###--------------------------------------------------------------------------------------###



###--------------------------------------------------------------------------------------###

### Activate LOOP modules



#Activate getorf module
module load emboss/6.6.0



#Activate diamond module
module load diamond/2.1.9



#Activate hmmer module
module load hmmer/3.3.2


#active module of seqtk to subset fastas
module load seqtk


#Activate conda module
module load miniconda3



#Activate conda module of cd-hit
eval "$(conda shell.bash hook)"

conda activate cdhit
echo "Modules Activated"

###--------------------------------------------------------------------------------------###

### LOOP



for fasta in *".fasta"

do

    #working on fasta

    echo "Working on -> $fasta"

    

    # Print sample directory NAME

    NAME=$(echo $fasta |  sed -E 's/\.contigs\.fasta$//')

    echo "Working on -> $NAME"    



    echo "#################################################################################"

    echo "####################  1. Geting Contigs without Diamond hit #####################"

    echo "#################################################################################"

    

    echo "$NAME.unmapped_contigs_diamond.txt"

    if [ -e "$NAME.unmapped_contigs_diamond.txt" ]; then

        echo "File $fasta exists in $TMPDIR"

    else

        echo "File $fasta does not exist in $TMPDIR"

    fi



    # retrieve information of interest

    awk '{print $1}' $NAME.unmapped_contigs_diamond.txt | sort -u > $NAME.contigs_hit.txt |

    grep -Fvf $NAME.contigs_hit.txt $fasta > $NAME.unmatch.fasta 

        

    # Filter by size the fasta    

    perl $RDRP/filterfastarange.fasta.pl 600 $NAME.unmatch.fasta > $NAME.unmatch.600nt.fasta



    # Clean the directory with temp files  

            

    echo "###############################################################################################"

    echo "####################  2. Translate ORFs (cut-off 200 aa) using GetORF #########################"

    echo "##############################################################################################"



    # Define the list of genetic code tables

    genetic_code_tables="1 3 4 5 6 11 16"



    # Loop through each genetic code table and run getorf

    for table in $genetic_code_tables; do

        echo "################## Finding ORFS for "$table" genetic code ##########################"

        output_file="ORFs_Frame${table}_"$NAME"_min600.fasta"

        getorf $NAME.unmatch.600nt.fasta -minsize 600 -table "$table" -find 0 -outseq "$output_file"

        echo "################### ORFS found for "$table" genetic code ##########################"
	
        # Add genetic table info to ORF headers
        sed -i "s/^>/>table${table}_/g" "$output_file"
	
	echo "################### Genetic Code "$table" added to fasta  ##########################"

    done



    # Concatenating multiple files in only one fasta and removing originals from getorf

    cat ORFs_Frame* > $NAME.Orfs_600.fasta

    rm ORFs_Frame*



    echo "###############################################################################################"
    echo "####################### 3. Remove redundant sequences using CD-Hit ############################"
    echo "##############################################################################################"

    #remove duplicated orfs
    cd-hit -i $NAME.Orfs_600.fasta -o $NAME.CD_hit_98.fasta -c 0.98 -d 1

    #substitute spaces per _
    sed -i.bak '/^>/ s/ /_/g' $NAME.CD_hit_98.fasta


    echo "###############################################################################################"
    echo "####################################### 4. HMM-RdRp hits ######################################"
    echo "##############################################################################################"

    # use hmmer to idd putative rdrp-scan
    hmmscan -E 1e-6 --cpu $NTs --domtblout $NAME.hmmscan.txt $RDRP/Profile_db_and_alignments/RdRp_HMM_profile_CLUSTALO.db.h3m $NAME.CD_hit_98.fasta

    #get fasta 
    #get a list with sequences matching
    grep -v "^#" $NAME.hmmscan.txt | awk '{print $4}' | sort -u > $NAME.hmmscan_contigs_hits.txt

    #subset from original fasta
	    
    seqtk subseq $NAME.CD_hit_98.fasta $NAME.hmmscan_contigs_hits.txt > RdRp_$NAME.fasta
    #python ~/viromapn/bin/Scripts/RdRp-scan/seq_extractor.py $NAME.hmmscan_contigs_hits.txt $NAME.CD_hit_98.fasta RdRp_$NAME.fasta
	
    echo "######################################################################################################"
    echo "####################################### 5. Organizing Results ######################################"
    echo "####################################################################################################"
    
    #### Make Directorys for Lib
    #chkmkdir $OUTPUT/$NAME
    chkmkdir $TMPDIR/$NAME

    ####Move Results of Interest to $NAME Directorys
    mv $NAME.CD_hit_98.fas* $TMPDIR/$NAME
    mv $NAME.unmatch.fasta $TMPDIR/$NAME
    mv RdRp_$NAME.fasta $TMPDIR/$NAME
    mv $NAME.hmmscan.txt $TMPDIR/$NAME
   # mv putativeRdrp_$NAME.fasta $TMPDIR/$NAME

    ####RSYNC TO OUTPUT
    rsync -tvr $TMPDIR/$NAME $OUTPUT/

    ####Remove All Files in $TMPDIR
    cd $TMPDIR
    rm -rf *$NAME*    

           
done



# put data in final file

#scp -r $TMPDIR/ $OUTDIR/
#rsync -tvr $TMPDIR/ $OUTPUT/


#Unload getorf module
module unload emboss/6.6.0


#Unload seqtk
module unload seqtk


#Unload diamond module
module unload diamond/2.1.9


#Unload hmmer module
module unload hmmer/3.3.2


#Deactivate Conda
conda deactivate


#Unload module of miniconda3
module unload miniconda3



###--------------------------------------------------------------------------------------###

### check the size of your temporary directory -- just to keep track

###--------------------------------------------------------------------------------------###



echo ; du -hs $TMPDIR

echo ; echo "### Final results will be found in $OUTDIR." ; ls -rtl $OUTDIR/



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
