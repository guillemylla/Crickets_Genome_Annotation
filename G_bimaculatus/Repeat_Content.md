# __Gryllus bimaculatus__ Repetitive Content Masking 

Here, I describe how I constructed a repeat library for *G. bimaculatus* and used to mask the genome.

## Requirments
 
* Repeat Library
	* MITE-tracker
	* LTRharvest + LTRDigest
	* RepeatModeler (de novo predicts repeats)
	* Transposons PSI
	* SINE database
	* RepBase

* RepeatMasker
 
 
## Repetitive elements library
 
Based on: http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/Repeat_Library_Construction-Advanced

### 0. prepare directories
mkdir TransposonPSI
mkdir MITEs
mkdir LTRs
mkdir RepeatModeler
mkdir Custom_Library


### 1.  LTRs (long terminal repeat) retrotransposons 

 - 1.1. LTRharvest: to collect LTRtransposons (http://genometools.org/documents/ltrharvest.pdf)
 - 1.2. LTRdigest: to filter LTRharvest results (http://genometools.org/documents/ltrdigest.pdf)
	
Prepare LTRharvest and LTRdigest

 ```bash
#srun -p test --pty --mem 10gb -t 0-01:00 /bin/bash 
cd ~/LTRs
wget http://genometools.org/pub/genometools-1.5.9.tar.gz
tar -zxvf  genometools-1.5.9.tar.gz
cd genometools-1.5.9
make with-hmmer=yes threads=yes
#tRNA database for LTRdifgest
 cd ~/LTRs
 wget http://gtrnadb2009.ucsc.edu/download/tRNAs/eukaryotic-tRNAs.fa.gz
 gunzip eukaryotic-tRNAs.fa.gz
```

Run as sbatch script


```bash
#!/bin/sh
#SBATCH --job-name=Run_LTRs    # Job name
#SBATCH --mail-type=ALL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=guillem_ylla@fas.harvard.edu     # Where to send mail	
#SBATCH -n 6                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -p shared   # Partition to submit toshared
#SBATCH --mem=15gb                     # Job memory request
#SBATCH --time=7-00:00               # Time limit days-hrs:min
#SBATCH -o out_LTRs_%j.log   # Standard output and error log
#SBATCH -e error_LTRs_%j.err   # Standard output and error log

pwd; hostname; date

echo "Set vars"

cd ~/LTRs/genometools-1.5.9

GenomeFasta="Gbimaculatus_Gap_filled.fasta"

echo "Run suffixerator"

./bin/gt suffixerator -db $GenomeFasta -indexname GbiV3index -dna -tis -suf -lcp -des -ssp

echo "suffixerator done"

echo "Run LtrHarvest with default parameters"

./bin/gt ltrharvest -index GbiV3index -v -out GbiV3LtrHarvest.out -outinner GbiV3LtrHarvest.outinner -gff3 GbiV3LtrHarvest.gff

echo "sort gff3"

./bin/gt gff3 -sort GbiV3LtrHarvest.gff > GbiV3LtrHarvest.gff.sort

echo "Run LTRdigest"

./bin/gt ltrdigest -trnas ~/LTRs/eukaryotic-tRNAs.fa \
	-outfileprefix GbiV3_LtrDigest \
	GbiV3LtrHarvest.gff.sort GbiV3index > GbiV3_LtrDigest.gff.dgt

echo "Done!"

pwd; hostname; date

```



### 2. MITE tracker

The installation of MITE is described for the annotation of V2.
 
 ```bash
#!/bin/sh
#SBATCH --job-name=MiteTracker    # Job name
#SBATCH --mail-type=ALL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=guillem_ylla@fas.harvard.edu     # Where to send mail	
#SBATCH -n 10                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -p general   # Partition to submit toshared
#SBATCH --mem=15gb                     # Job memory request
#SBATCH --time=7-00:00               # Time limit days-hrs:min
#SBATCH -o out_MiteTracker_%j.log   # Standard output and error log
#SBATCH -e error_MiteTracker_%j.err   # Standard output and error log

pwd; hostname; date

echo "Load Modules"
module load Anaconda/5.0.1-fasrc02
module load blast/2.6.0+-fasrc01
module load python/3.6.3-fasrc02


echo "Load python-conda environment"
source activate MITEtrackerenvironment
echo "Activated"

cd MITE-Tracker



echo "Run MITE tracker"

GenomeFasta="~/Gbimaculatus_Gap_filled.fasta"

python -m MITETracker.py -g $GenomeFasta -w 10 -j GbiV3_Mitetracker

echo "result should be at: /n/home09/gylla/Software/MITE-Tracker/results"

source deactivate

echo "Done!"

pwd; hostname; date

```
 
Collect the results:

```bash
cp /n/home09/gylla/Software/MITE-Tracker/results/GbiV3_Mitetracker/all.fasta .
```
 
### 3. RepeatModeler




 ```bash
#!/bin/sh
#SBATCH --job-name=runClass    # Job name
#SBATCH --mail-type=ALL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH -n 30               # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH --mem=70gb                     # Job memory request
#SBATCH --time=30-00:00               # Time limit days-hrs:min
#SBATCH -o out_class_%j.log   # Standard output and error log
#SBATCH -e error_class_%j.err   # Standard output and error log
pwd; hostname; date
module load repeatmodeler/1.0.8

GenomeFasta="Gbimaculatus_Gap_filled.fasta"
BuildDatabase -name GbiV3database -engine ncbi $GenomeFasta

echo "Run repeatModeler"
RepeatModeler -pa $SLURM_NTASKS -database GbiV3database 
echo "Done!"

pwd; hostname; date
```






### 4. TransposonPSI

Installed in V2.1 I will run it as a job array (to make it faster) afetr splitting the genome.

*~/TransposonPSI*

Split genome
 
 ```bash
srun -p test --pty --mem 5gb -n 2 -N 1 -t 0-06:00 /bin/bash
module load ucsc/20150820-fasrc01

mkdir GbiV3_splitted 

cd GbiV3_splitted

faSplit sequence ~/Gbimaculatus_Gap_filled.fasta  10 GbiV3_part

```

Create the Job array

```bash
#!/bin/sh
#SBATCH --job-name=GbiTransposonPSI    # Job name
#SBATCH --mail-type=ALL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH -n 10                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -p shared   # Partition to submit toshared
#SBATCH --mem=10gb                     # Job memory request
#SBATCH --time=7-00:00               # Time limit days-hrs:min
#SBATCH -o out_TransposonPSI_%A_%a.out   # Standard output and error log
#SBATCH -e error_TransposonPSI_%A_%a.out   # Standard output and error log

module load Anaconda/5.0.1-fasrc02

source activate MITEtrackerenvironment

cd ~/TransposonPSI

##from 00 to 09
GenomeFasta=~/TransposonPSI/GbiV3_splitted/GbiV3_part0"${SLURM_ARRAY_TASK_ID}".fa

perl transposonPSI.pl $GenomeFasta nuc

source deactivate


```

And later run the Job Array:

```bash
##from 00 to 09
sbatch --array=0-9  run_TransposonPSI.sbatch 
```


The file *GbiV3_part09.fa* run out of time..  I can't continuie from where it got interrupted... I split and re-run...

 ```bash
srun -p test --pty --mem 5gb -n 2 -N 1 -t 0-06:00 /bin/bash
module load ucsc/20150820-fasrc01

cd ~/TransposonPSI/GbiV3_splitted
faSplit sequence GbiV3_part09.fa  10 GbiV3_part09_repart

```
Re run the GbiV3_part09.fa parted in 10...

```bash
##from 00 to 09
 sbatch --array=0-9 run_TransposonPSI_part09reparted.sbatch
 ```




 When it finishes, I collect results, transform to gff3 and get fasta
 
 ```bash
ls -lah  *.fa.TPSI.allHits
ls -lah  *.fa.TPSI.allHits | wc -l
#19 files : 00 to 08 (9 files) + 00 to 09 (10 files)

cat *.fa.TPSI.allHits > GbiV3_joinparts.fa.TPSI.allHits
wc -l  GbiV3_joinparts.fa.TPSI.allHits
#60814

## I created the script https://github.com/guillemylla/Gryllus_genome_annotation/blob/master/Annotation_v1/Repetitive_content/Scripts/TransposonPSI2gff3.py
python TransposonPSI2gff3.py GbiV3_joinparts.fa.TPSI.allHits > GbiV3_joinparts.fa.TPSI.allHits.gff3  

wc -l GbiV3_joinparts.fa.TPSI.allHits.gff3  
#60814

module load bedtools2/2.26.0-fasrc01
GenomeFasta="~/Gbimaculatus_Gap_filled.fasta"
bedtools getfasta -name -s -bed GbiV3_joinparts.fa.TPSI.allHits.gff3   -fi $GenomeFasta -fo GbiV3_joinparts.TPSI.allHits.fa

 ```




### 5. Join libraries

I copy the ourput of each softare to a new directory

```bash
cd ~/Custom_Library

# Extract taxon-specific repeats from RepeatMasker database (RepBase)
module load Anaconda/5.0.1-fasrc02
source activate CondarepeatMasker
perl queryRepeatDatabase.pl -clade insecta  > MakeRepeatBaseInsecta.fa
## RepeatMaskerDB
tail -n+2 MakeRepeatBaseInsecta.fa > RepeatMaskerInsecta.lib
rm MakeRepeatBaseInsecta.fa
source deactivate 
 
# SINE databse
 wget http://sines.eimb.ru/banks/SINEs.bnk 
 
# ## LTRdigest (after LTRharvest)
cp ../LTRs/genometools-1.5.9/GbiV3_LtrDigest_complete.fas .
## MITE Tracker
cp ../MITEs/all.fasta  MiteTracker.lib
##TransposonPSI
cp  ~/TransposonPSI/GbiV3_joinparts.TPSI.allHits.fa TransposonPSI_GBIv3.lib
## RepeatModeler
cp ~/RepeatModeler/consensi.fa.classified RepeatModeler_GBIv3consensi.fa.classified.lib

```
 
 And now I join the libraries:
 
 ```bash
cat  RepeatMaskerInsecta.lib \
  SINEs.bnk \
  GbiV3_LtrDigest_complete.fas \
  MiteTracker.lib \
  TransposonPSI_GBIv3.lib \
  RepeatModeler_GBIv3consensi.fa.classified.lib > CombinedLibrary.lib
  
grep -c ">"  CombinedLibrary.lib 
#126437

 ```
 
  
### 6. Filter, remove redundancies and classify

```bash
srun -p test --pty --mem 5gb -n 2 -N 1 -t 0-06:00 /bin/bash

##Removing short seqs (>50nts)
module load seqtk/1.2-fasrc01
seqtk seq -L 50 CombinedLibrary.lib > CombinedLibrary.lib.minlen50
grep -c ">"  CombinedLibrary.lib.minlen50
#126410 Only removes 27

##Removing redundancies (like in https://www.nature.com/protocolexchange/protocols/6761#/procedure)
module load usearch/9.2.64-fasrc01
usearch -cluster_fast CombinedLibrary.lib.minlen50 -id 0.8 -consout CombinedLibrary.lib.minlen50.nr

##output
#   Seqs  124167 (124.2k)
#   Clusters  87301 (87.3k)
#   Max size  850
#   Avg size  1.4
#   Min size  1
# Singletons  80211 (80.2k), 64.6% of seqs, 91.9% of clusters

grep -c ">" CombinedLibrary.lib.minlen50.nr
#87301


```
Run Repeat Classifier as Job array.

```bash
srun -p test --pty --mem 5gb -n 2 -N 1 -t 0-06:00 /bin/bash
module load ucsc/20150820-fasrc01

mkdir Split_Combinedlibrary
cp CombinedLibrary.lib.minlen50.nr Split_Combinedlibrary/
cd Split_Combinedlibrary

faSplit sequence CombinedLibrary.lib.minlen50.nr 20 CombinedLibrary.nr

mkdir directory_{00..19}


for i in `seq 0 9`;   do echo "mv CombinedLibrary0"$i".fa directory_0"$i'/' ;  done
for i in `seq 10 19`;    do echo "mv CombinedLibrary"$i".fa directory_"$i'/' ;  done
# Do it, move
for i in `seq 0 9`;    do mv "CombinedLibrary0"$i".fa"  "directory_0"$i'/' ;  done
for i in `seq 10 19`;  do mv "CombinedLibrary"$i".fa"  "directory_"$i'/' ;  done

### Odyssey 
 sbatch --array=0-4  run_RepeatClassifier_jobarray_Zero.sbatch
 sbatch --array=5-9  run_RepeatClassifier_jobarray_Zero.sbatch
 sbatch --array=10-14  run_RepeatClassifier_jobarray.sbatch
 sbatch --array=15-19  run_RepeatClassifier_jobarray.sbatch

 find . -type f -name "*.fa.classified"  | wc -l
 find . -type f -name "*.fa.classified"  
 find . -type f -name "*.fa.classified"  -exec cat {} \; > ../CombinedLibrary.lib.minlen50.nr.classified.fa
 
 
 grep -c ">" CombinedLibrary.lib.minlen50.nr.classified.fa
 #87296 (87301 minus some empty seqd)
 ```


## 7. Remove possible  protein-coding genes from non-transposable element 

In order to remove posible  protein-coding genes from non-transposable element , we filter out all thos elements classified as "Unknown" by the RepeatClassifier
and with a Blastx hit (evalue'<'1e-10) to Insect proteins from the reviewed Swiss-Prot database.
 
I download 9.229 REVIEWED protein seqeunces of Insects from uniprot/Swiss-Prot (February 13 2019)
 
 (https://www.uniprot.org/uniprot/?query=reviewed:yes%20taxonomy:50557)
 
 
 File: uniprot-reviewed%3Ayes+taxonomy%3A50557.fasta

 
```bash

## Get sequences classified as Unknown
srun -p test --pty --mem 20gb -n 4 -N 1 -t 0-06:00 /bin/bash

module load ucsc/20150820-fasrc01
module load blast/2.6.0+-fasrc01

cd ~/Custom_Library/Filter_CutomLibBlast 


##copy here the classified custom library
cp ../CombinedLibrary.lib.minlen50.nr.classified.fa .
grep -c ">" CombinedLibrary.lib.minlen50.nr.classified.fa
#87296

grep "Unknown" CombinedLibrary.lib.minlen50.nr.classified.fa > UnknownIdslist.txt
wc -l UnknownIdslist.txt
#15537
sed -i 's/>//g' UnknownIdslist.txt

#get unknowns fasta
faSomeRecords CombinedLibrary.lib.minlen50.nr.classified.fa UnknownIdslist.txt UnknownRepeats.fasta
 
#prepare uniprot database
# I previously downlaoded the database (v1)
cp -r /n/regal/extavour_lab/gylla/Annotation_v1/Repetitive_content/Custom_library_v1/FilterDb/UniprotSeqs/ .
#makeblastdb -in UniprotSeqs/uniprot-reviewed%3Ayes+taxonomy%3A50557.fasta  -dbtype prot
 
 
blastx -query UnknownRepeats.fasta -db UniprotSeqs/uniprot-reviewed%3Ayes+taxonomy%3A50557.fasta -evalue 1e-10 \
 	-num_threads 4 -max_target_seqs 1 -outfmt '6 qseqid sseqid evalue bitscore sgi sacc stitle' -out Blast_out.txt
 
 
awk '{print $1}' Blast_out.txt  | sort | uniq | wc -l 

awk -F "\t" '{print $1,$7}' Blast_out.txt  | sort | uniq
# 835
awk -F "\t" '{print $1,$7}' Blast_out.txt  | sort | uniq | grep -i -v "transposon" | grep -i -v "Copia protein" | grep -i -v "mobile element" | grep -i -v "transposable"  | grep -i -v "transposase" | wc -l
# 68
awk -F "\t" '{print $1,$7}' Blast_out.txt  | sort | uniq | grep -i -v "transposon" | grep -i -v "Copia protein" | grep -i -v "mobile element" | grep -i -v "transposable"  | grep -i -v "transposase" | awk '{print $1}' > Unknowns_with_Port_hit.txt
#68

## Filter Cutsom library
faSomeRecords -exclude CombinedLibrary.lib.minlen50.nr.classified.fa Unknowns_with_Port_hit.txt Gbi_Rep_CombinedLibrary.lib.minlen50.nr.classified.filtered.fa

grep -c ">" Gbi_Rep_CombinedLibrary.lib.minlen50.nr.classified.filtered.fa
# 87228 (87296-68=767)


## check...
wc -l Unknowns_with_Port_hit.txt
#68
grep -c ">" CombinedLibrary.lib.minlen50.nr.classified.fa
#87296
grep -c ">" Gbi_Rep_CombinedLibrary.lib.minlen50.nr.classified.filtered.fa
#87228

cp  Gbi_Rep_CombinedLibrary.lib.minlen50.nr.classified.filtered.fa ../
```

 The final custom repeat library is:  Gbi_Rep_CombinedLibrary.lib.minlen50.nr.classified.filtered.fa

Path: "*~/Custom_Library/Gbi_Rep_CombinedLibrary.lib.minlen50.nr.classified.filtered.fa*"

Backup: *~/Gbi_Rep_CombinedLibrary.lib.minlen50.nr.classified.filtered.fa*




### 8. RepeatMasker

MAKERE was not producing the repeat content report... I run RepeatMasker first and provide it to maker later.


```bash
GenomeFasta="~/Gbimaculatus_Gap_filled.fasta"
RepLib="~/Custom_Library/Gbi_Rep_CombinedLibrary.lib.minlen50.nr.classified.filtered.fa"

cd ~/RepeatMasker

RepeatMasker  -pa $SLURM_NTASKS -lib $RepLib $GenomeFasta -gff -dir ~/RepeatMasker

```
 
 I ran RepeatMasker with 30 cpus, and took ~2days (46h).
 

### Backup
 
I will collect all scripts and essential files, donwload in local and upload to github

```bash 
cd ~

## List fo sbatch files
./Custom_Library/Split_Combinedlibrary/run_RepeatClassifier_jobarray_Zero.sbatch
./Custom_Library/Split_Combinedlibrary/run_RepeatClassifier_jobarray.sbatch
./Custom_Library/Split_Combinedlibrary/run_RepeatClassifier_repeat11.sbatch
./MITEs/run_MITEtracker.sbatch
./RepeatMasker/run_RepeatMasker.sbatch
./TransposonPSI/run_TransposonPSI_part09reparted.sbatch
./TransposonPSI/run_TransposonPSI.sbatch
./LTRs/run_LTR_harvest_digest.sbatch


find . -type f -name "*.sbatch" -exec cp {} Sbatch_cripts_bkup/ \;



scp -r gylla@odyssey.rc.fas.harvard.edu:~/Sbatch_cripts_bkup .

##

```

Download custom repeat library in local and later upload it to github using "git lfs":

```bash
scp -r gylla@odyssey.rc.fas.harvard.edu:~/Custom_Library/Gbi_Rep_CombinedLibrary.lib.minlen50.nr.classified.filtered.fa .

# git lfs instructions
##Local
git lfs track "*.fa"
git add .gitattributes
git commit -m "track *.fa files using Git LFS"
git add Gbi_Rep_CombinedLibrary.lib.minlen50.nr.classified.filtered.fa
git commit -m "GBI Custom Repeat library"
git push origin master
```

