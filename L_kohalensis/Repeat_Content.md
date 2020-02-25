# __Laupala kohalensis__ Repetitive Content Masking

This document contains all the commands used for masking the repetitive content of *Laupala_kohalensis* followin the same pipeline applied to __Gryllus bimaculatus__ genome. [See here](https://github.com/guillemylla/Gryllus_genome_annotation/blob/master/GBI_Genome_v3/Repetitive_content/README.md).


For more details about the methodology please see tha [*G. bimaculatus* pipline](https://github.com/guillemylla/Gryllus_genome_annotation/blob/master/GBI_Genome_v3/Repetitive_content/README.md). 


## Requirments

* Repeat Library
	* MITE-tracker
	* LTRharvest + LTRDigest
	* RepeatModeler (de novo predicts repeats)
	* Transposons PSI
	* SINE database
	* RepeatMasker (RepBase)


### 1.  LTRs (long terminal repeat) retrotransposons 

 - 1.1. LTRharvest: to collect LTRtransposons (http://genometools.org/documents/ltrharvest.pdf)
 - 1.2. LTRdigest: to filter LTRharvest results (http://genometools.org/documents/ltrdigest.pdf)


```bash
#srun -p test --pty --mem 10gb -t 0-01:00 /bin/bash#Install
cd ~/Laupala_kohalensis_annotation/Repetitive_content
mkdir LTRs
wget http://genometools.org/pub/genometools-1.5.9.tar.gz
tar -zxvf  genometools-1.5.9.tar.gz
cd genometools-1.5.9
make with-hmmer=yes threads=yes
#tRNA database for LTRdifgest
cd ~/Laupala_kohalensis_annotation/Repetitive_content/LTRs
wget http://gtrnadb2009.ucsc.edu/download/tRNAs/eukaryotic-tRNAs.fa.gz
gunzip eukaryotic-tRNAs.fa.gz

## Run in several sbtach scripts (~/Annotation_v1/sbatch_scripts_v2.1)
cd ~/Annotation_v1/Repetitive_content/LTRs/genometools-1.5.9

GenomeFasta="~/Laupala_kohalensis/GCA_002313205.1_ASM231320v1_genomic.fna"

#Get index
./bin/gt suffixerator -db $GenomeFasta -indexname Lkoindex -dna -tis -suf -lcp -des -ssp

#Run LtrHarvest with default parameters
./bin/gt ltrharvest -index Lkoindex -v -out LkoLtrHarvest.out -outinner LkoLtrHarvest.outinner -gff3 LkoLtrHarvest.gff

#Run LTRdigest with default parameters
./bin/gt gff3 -sort LkoLtrHarvest.gff > LkoLtrHarvest.gff.sort

./bin/gt ltrdigest -trnas ~/Annotation_v1/Repetitive_content/LTRs/eukaryotic-tRNAs.fa \
	-outfileprefix Lko_LtrDigest \
	LkoLtrHarvest.gff.sort Lkoindex > Lko_LtrDigest.gff.dgt


## cp ~/Annotation_v1/Repetitive_content/LTRs/genometools-1.5.9/Lko_LtrDigest*  ~/Laupala_kohalensis_annotation/Repetitive_content/LTRs/


```
### 2. MITE tracker

```bash
## RUN as sbatch
module load Anaconda/5.0.1-fasrc02
module load blast/2.6.0+-fasrc01
module load python/3.6.3-fasrc02


source activate ~/.conda/envs/MITEtrackerenvironment

cd ~/Software/MITE-Tracker

GenomeFasta="~/Laupala_kohalensis/GCA_002313205.1_ASM231320v1_genomic.fna"
python -m MITETracker.py -g $GenomeFasta -w 10 -j Lko_Mitetracker

source deactivate

```
### 3. RepeatModeler


```bash
GenomeFasta="GCA_002313205.1_ASM231320v1_genomic.fna"

BuildDatabase -name Lkodatabase -engine ncbi $GenomeFasta

RepeatModeler -pa 30 -database Lkodatabase 

```



### 4. TransposonPSI


I split the genome and run it in fragments as job array to speed it up.

```bash
srun -p test --pty --mem 5gb -n 2 -N 1 -t 0-06:00 /bin/bash
module load ucsc/20150820-fasrc01

mkdir Lko_genome_splitted 

cd Lko_genome_splitted

faSplit sequence ~/Laupala_kohalensis/GCA_002313205.1_ASM231320v1_genomic.fna 9 Lko_genome_part

```

Create sbatch file for running as Job Array:

```bash
#!/bin/sh
#SBATCH --job-name=TransposonPSI    # Job name
#SBATCH --mail-type=ALL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH -n 9                # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -p shared   # Partition to submit toshared
#SBATCH --mem=10gb                     # Job memory request
#SBATCH --time=7-00:00               # Time limit days-hrs:min
#SBATCH -o out_TransposonPSI_%A_%a.out   # Standard output and error log
#SBATCH -e error_TransposonPSI_%A_%a.out   # Standard output and error log

module load Anaconda/5.0.1-fasrc02

source activate MITEtrackerenvironment

cd ~/Laupala_kohalensis_annotation/Repetitive_content/Transposonpsi

GenomeFasta=~/Laupala_kohalensis_annotation/Repetitive_content/Transposonpsi/Lko_genome_splitted/Lko_genome_part"${SLURM_ARRAY_TASK_ID}".fa

perl ~/Software/TransposonPSI_08222010/transposonPSI.pl $GenomeFasta nuc

source deactivate

echo "outputs at: ~/Software/TransposonPSI_08222010/"

```

And later run the Job Array:

```bash

sbatch --array=0-8  4.1-TransposonPSI_ARRAY.sbatch

```

ERROR in Lko_genome_part0.fa Contain too many sequences and after 7 days didn't finish...

I need to re-split this fragment

```bash
srun -p test --pty --mem 5gb -n 2 -N 1 -t 0-06:00 /bin/bash
module load ucsc/20150820-fasrc01

cd ~/Laupala_kohalensis_annotation/Repetitive_content/Transposonpsi/Lko_genome_splitted

mkdir Lko_genome_part0_resplitted 

cd ~/Laupala_kohalensis_annotation/Repetitive_content/Transposonpsi/Lko_genome_splitted/Lko_genome_part0_resplitted

faSplit sequence ../Lko_genome_part0.fa 20 GbiV3_part0_resplited_

sbatch --array=0-9  run_TransposonPSI_ARRAY_0_9.sbatch
sbatch --array=0-9  run_TransposonPSI_ARRAY_10_19.sbatch

```

Still, 2 files GbiV3_part0_resplited_00 and GbiV3_part0_resplited_01  didn't finish..


```bash

#lets see where did they get stack:
tail GbiV3_part0_resplited_00.fa.TPSI.allHits
#last NNCF01032388.1
tail GbiV3_part0_resplited_01.fa.TPSI.allHits
#last NNCF01073432.1

## get the line in fasta where got stcked
grep -n "NNCF01032388.1" GbiV3_part0_resplited_00.fa
#192833
grep -n "NNCF01073432.1" GbiV3_part0_resplited_01.fa
#208935

# put the reminding sequences in rnew files
tail -n +192233 GbiV3_part0_resplited_00.fa | head
tail -n +208935 GbiV3_part0_resplited_01.fa | head

tail -n +192233 GbiV3_part0_resplited_00.fa > GbiV3_part0_resplited_00_relaunchlastlines.fa
tail -n +208935 GbiV3_part0_resplited_01.fa > GbiV3_part0_resplited_01_relaunchlastlines.fa

grep -c ">" GbiV3_part0_resplited_00_relaunchlastlines.fa
#9972
grep -c ">" GbiV3_part0_resplited_01_relaunchlastlines.fa
#3662

sbatch run_TransposonPSI_00rerunlastlines.sbatch
sbatch run_TransposonPSI_01rerunlastlines.sbatch
```

When it finishes, I collect results, transform to gff3 and get fasta
 
 
 ```bash
 ## the output of the file 0 resplited in 20 (00 to 19):
cd ~/Laupala_kohalensis_annotation/Repetitive_content/Transposonpsi

## the output from 01 to 08  (00 crashed) + the 00 after joining
ls -lah *.fa.TPSI.allHits
## Part 00 resplitted and, the re-run of part0_00 and part0_01 (GbiV3_part0_resplited_00_relaunchlastlines GbiV3_part0_resplited_01_relaunchlastlines)
ls -lah Lko_genome_splitted/Lko_genome_part0_resplitted/*.fa.TPSI.allHits
# join all the parts 0
cat Lko_genome_splitted/Lko_genome_part0_resplitted/*.fa.TPSI.allHits > Lko_genome_part0_JOINED.fa.TPSI.allHits

# join all the parts
cat *.fa.TPSI.allHits > Lko_joinparts.fa.TPSI.allHits

wc -l Lko_joinparts.fa.TPSI.allHits
#119423 Lko_joinparts.fa.TPSI.allHits

## I use the script https://github.com/guillemylla/Gryllus_genome_annotation/blob/master/Annotation_v1/Repetitive_content/Scripts/TransposonPSI2gff3.py
python TransposonPSI2gff3.py Lko_joinparts.fa.TPSI.allHits > Lko_joinparts.fa.TPSI.allHits.gff3  

wc -l Lko_joinparts.fa.TPSI.allHits.gff3
#119423

module load bedtools2/2.26.0-fasrc01
GenomeFasta="~/Laupala_kohalensis/GCA_002313205.1_ASM231320v1_genomic.fna"
bedtools getfasta -name -s -bed Lko_joinparts.fa.TPSI.allHits.gff3  -fi $GenomeFasta -fo LKO_TPSI.allHits.fa

grep -c ">" LKO_TPSI.allHits.fa
#119423

 ```


### 5. Join libraries

All the softwares outputs in the same directory

```bash
srun -p test --pty --mem 5gb -n 2 -N 1 -t 0-06:00 /bin/bash

cd ~/Laupala_kohalensis_annotation/Repetitive_content/Custom_library_lko

 module load RepeatMasker/4.0.5-fasrc05
# Extract taxon-specific repeats from RepeatMasker database (RepBase)
#perl /n/helmod/apps/centos7/Core/RepeatMasker/4.0.5-fasrc05/util/queryRepeatDatabase.pl -clade insecta  > MakeRepeatBaseInsecta.fa
cp ~/Annotation_v1/Repetitive_content/Custom_library_v1/MakeRepeatBaseInsecta_2.fa  .
## LTRdigest (after LTRharvest)
cp ../LTRs/Lko_LtrDigest_complete.fas LtrDigest_lko.lib
## MITE Tracker
cp ~/Software/MITE-Tracker/results/Lko_Mitetracker/all.fasta MiteTracker_lko.lib
## SINE databse
cp ~/Annotation_v1/Repetitive_content/Custom_library_v1/SINEs.bnk .
## RepeatModeler
cp ../RM_12173.SatFeb231247572019/consensi.fa.classified RepeatModeler_consensi_classified_Lko.lib
###TransposonPSI
cp ../Transposonpsi/LKO_TPSI.allHits.fa .
```

Join Libraries

```bash
cat LtrDigest_lko.lib \
	MiteTracker_lko.lib \
	MakeRepeatBaseInsecta_2.fa \
	RepeatModeler_consensi_classified_Lko.lib \
	SINEs.bnk  \
	LKO_TPSI.allHits.fa > CombinedLibrary_Lko.lib

grep -c ">" CombinedLibrary_Lko.lib
#185834

```

### 6. Filter, remove redundancies and classify


```bash
srun -p test --pty --mem 5gb -n 2 -N 1 -t 0-06:00 /bin/bash
##Removing short seqs (>50nts)
module load seqtk/1.2-fasrc01
seqtk seq -L 50 CombinedLibrary_Lko.lib > CombinedLibrary_Lko.lib.minlen50
grep -c ">"  CombinedLibrary_Lko.lib.minlen50
#185767
```

### 6. Remove redundancies and classify

```bash
##Removing redundancies (like in https://www.nature.com/protocolexchange/protocols/6761#/procedure)
module load usearch/9.2.64-fasrc01
usearch -cluster_fast CombinedLibrary_Lko.lib.minlen50 -id 0.8 -consout  CombinedLibrary_Lko.lib.minlen50.nr

#     Seqs  180646 (180.6k)
# Clusters  112078 (112.1k)


grep -c ">"  CombinedLibrary_Lko.lib.minlen50.nr
#112078

```

Run Repeat Classifier as Job array.

```bash
#####
#To make it fast, I do fasta chunk and run it by chunk
module load ucsc/20150820-fasrc01
mkdir Split_Combinedlibrary 
cp CombinedLibrary_Lko.lib.minlen50.nr Split_Combinedlibrary/
cd Split_Combinedlibrary

faSplit sequence CombinedLibrary_Lko.lib.minlen50.nr 20 LkoCombinedLibrary.nr

mkdir directory_{00..19}

#try
for i in `seq 0 9`;   do echo "mv LkoCombinedLibrary0"$i".fa directory_0"$i'/' ;  done  
for i in `seq 10 19`;    do echo "mv LkoCombinedLibrary"$i".fa directory_"$i'/' ;  done  
# Do it, move
for i in `seq 0 9`;    do mv "LkoCombinedLibrary0"$i".fa"  "directory_0"$i'/' ;  done  
for i in `seq 10 19`;  do mv "LkoCombinedLibrary"$i".fa"  "directory_"$i'/' ;  done  


sbatch --array=0-4  run_RepeatClassifier_jobarray_ZERO.sbatch
sbatch --array=5-9  run_RepeatClassifier_jobarray_ZERO.sbatch
sbatch --array=10-14  run_RepeatClassifier_jobarray.sbatch
sbatch --array=15-19  run_RepeatClassifier_jobarray.sbatch




find . -type f -name "*.fa.classified"  | wc -l
# 20 files
find . -type f -name "*.fa.classified"  
find . -type f -name "*.fa.classified"  -exec cat {} \; > CombinedLibrary_Lko.lib.minlen50.nr.classified.fa


#check
grep -c ">" CombinedLibrary_Lko.lib.minlen50.nr
#112078
grep -c ">" CombinedLibrary_Lko.lib.minlen50.nr.classified.fa
#112041 # there were a few empty seqs...


```


### 7. Remove possible  protein-coding genes from non-transposable element 

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

cd ~/Laupala_kohalensis_annotation/Repetitive_content/Custom_library_lko/Filter_CutomLibBlast


##copy here the classified custom library
cp ../CombinedLibrary_Lko.lib.minlen50.nr.classified.fa .
grep -c ">" CombinedLibrary_Lko.lib.minlen50.nr.classified.fa
#112041

grep "Unknown" CombinedLibrary_Lko.lib.minlen50.nr.classified.fa > UnknownIdslist.txt
wc -l UnknownIdslist.txt
#10737
sed -i 's/>//g' UnknownIdslist.txt

#get unknowns fasta
faSomeRecords CombinedLibrary_Lko.lib.minlen50.nr.classified.fa UnknownIdslist.txt UnknownRepeats.fasta

grep -c ">" UnknownRepeats.fasta
#10737

#prepare uniprot database
# I previously downlaoded the database (v1)
cp -r ~/Annotation_v1/Repetitive_content/Custom_library_v1/FilterDb/UniprotSeqs/ .
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
faSomeRecords -exclude CombinedLibrary.lib.minlen50.nr.classified.fa Unknowns_with_Port_hit.txt CombinedLibrary.lib.minlen50.nr.classified.filtered.fa

grep -c ">" CombinedLibrary.lib.minlen50.nr.classified.filtered.fa
# 87228 (87296-68=767)


## check...
wc -l Unknowns_with_Port_hit.txt
#68
grep -c ">" CombinedLibrary.lib.minlen50.nr.classified.fa
#87296
grep -c ">" CombinedLibrary.lib.minlen50.nr.classified.filtered.fa
#87228

cp  CombinedLibrary.lib.minlen50.nr.classified.filtered.fa ../
```

 The final custom repeat library is:  CombinedLibrary.lib.minlen50.nr.classified.filtered.fa

Path: "*~/GBI_Genome_v3/Repetitive_content/Custom_Library *"


The final custom repeat library is: CombinedLibrary.nr.classified.filtered.lib


### 8. RepeatMasker

MAKERE was not producing the repeat content report... I run RepeatMasker first and provide it to maker later.


```bash
GenomeFasta="~/Laupala_kohalensis/GCA_002313205.1_ASM231320v1_genomic.fna"
RepLib="~/Laupala_kohalensis_annotation/Repetitive_content/Custom_library_lko/CombinedLibrary_Lko.lib.minlen50.nr.classified.fa"

cd ~/Laupala_kohalensis_annotation/Repetitive_content/RepeatMasker

RepeatMasker  -pa $SLURM_NTASKS -lib $RepLib $GenomeFasta -gff -dir ~/Laupala_kohalensis_annotation/Repetitive_content/RepeatMasker

```


RepeatMasker output

```bash

==================================================
file name: GCA_002313205.1_ASM231320v1_genomic.fna
sequences:        148784
total length: 1595214429 bp  (1563778341 bp excl N/X-runs)
GC level:         35.58 %
bases masked:  566518287 bp ( 35.51 %)
==================================================
               number of      length   percentage
               elements*    occupied  of sequence
--------------------------------------------------
SINEs:            29510      7083717 bp    0.44 %
      ALUs          304       101257 bp    0.01 %
      MIRs         1248       430584 bp    0.03 %

LINEs:           1035151    322470849 bp   20.21 %
      LINE1         941       367057 bp    0.02 %
      LINE2      584526    167380843 bp   10.49 %
      L3/CR1      10257      4624100 bp    0.29 %

LTR elements:     57347     29690552 bp    1.86 %
      ERVL          231        43500 bp    0.00 %
      ERVL-MaLRs      0            0 bp    0.00 %
      ERV_classI   1821       585650 bp    0.04 %
      ERV_classII   389       125302 bp    0.01 %

DNA elements:    189815     62384975 bp    3.91 %
     hAT-Charlie  15008      5154516 bp    0.32 %
     TcMar-Tigger  8896      2459752 bp    0.15 %

Unclassified:    409303    128822550 bp    8.08 %

Total interspersed repeats:550452643 bp   34.51 %


Small RNA:        13816      3005585 bp    0.19 %

Satellites:        2088       882748 bp    0.06 %
Simple repeats:  307925     19782955 bp    1.24 %
Low complexity:   48386      2381730 bp    0.15 %
==================================================

* most repeats fragmented by insertions or deletions
  have been counted as one element
                                                      

The query species was assumed to be homo          
RepeatMasker version open-4.0.5 , default mode
                                   
run with rmblastn version 2.2.27+
The query was compared to classified sequences in ".../CombinedLibrary_Lko.lib.minlen50.nr.classified.fa"
RepBase Update 20160829, RM database version 20160829


```




