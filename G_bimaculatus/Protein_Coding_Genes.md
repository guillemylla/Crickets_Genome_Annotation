# Annotation of *Gryllus bimaculatus* protein-coding genes V3

Protein coding genes of the *Gryllus bimaculatus* genome version V3. 

## Ab-initio gene predictors

### BUSCO - Augustus
 
 I ran Busco at arthtopod and Insect levels.
 
 Busco creates a "augustus_output" that can be passed to MAKER. We need the option --long, which takes longer time but optimize the HMM search model to train Augustus and produce a trained HMM for MAKER. (https://gist.github.com/darencard/bb1001ac1532dd4225b030cf0cd61ce2)
 
 Later, the files must be renamed for maker (https://wiki.cyverse.org/wiki/display/TUT/Training+ab+initio+Gene+Predictors+for+MAKER+genome+annotation & https://gist.github.com/darencard/bb1001ac1532dd4225b030cf0cd61ce2)
 
 
 
```bash
 Run as sbatch jobs

#Insect
Genome="Gbimaculatus_Gap_filled.fasta"
outfilename="BuscoInsecta_GBIv3"
python ~/busco/scripts/run_BUSCO.py --cpu $SLURM_NTASKS --long -i $Genome  -o $outfilename -l ~/busco/datasets/insecta_odb9 -m geno

#Arthropod 
Genome="Gbimaculatus_Gap_filled.fasta"
outfilename="BuscoArthropoda_GBIv3"
python ~/busco/scripts/run_BUSCO.py --cpu $SLURM_NTASKS --long -i $Genome  -o $outfilename -l ~/busco/datasets/arthropoda_odb9 -m geno


# ############ Get Augustus files trained by Busco
# Before running maker I need to put the training files under the Augustus   config/species directory
```

 
### GeneMark-ES
 
 A *de-novo* gene predictor in self-training mode. (full GeneMark installation on Annotation_v1)
 

Run Gmark as sbatch script (*run_GeneMArk.sbatch*)

```bash
#!/bin/sh
#SBATCH --job-name=runGeneMark   # Job name
#SBATCH --mail-type=ALL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH -n 20               # Number of cores
#SBATCH -N 1                # Ensure that all cores are on one machine
#SBATCH -p shared   # Partition to submit toshared
#SBATCH --mem=20gb                     # Job memory request
#SBATCH --time=7-00:00               # Time limit days-hrs:min
#SBATCH -o out_GeneMark_%j.log   # Standard output and error log
#SBATCH -e err_GeneMark_%j.err   # Standard output and error log

pwd; hostname; date

echo "Load modules"
module load perl-modules/5.26.1-fasrc01
echo "Modules loaded"

cd "GeneMark"

echo "Run GeneMark"

Genome="Gbimaculatus_Gap_filled.fasta"


perl gmes_petap/gmes_petap.pl --ES \
        --cores  $SLURM_NTASKS \
        --sequence $Genome

echo "Genemark finished!"

echo "Done!"

date
 
```


## Repteat  Content

The Repeat content library was produced following the pipline described here [Repeat_Content.md](Repeat_Content.md).



## Prepare RNA-seq Hisat2


Run Hisat2 with --dta option : *" --dta== --downstream-transcriptome-assembly. HISAT2 provides options for transcript assemblers (e.g., StringTie and Cufflinks) to work better with the alignment from HISAT2 (see options such as --dta and --dta-cufflinks"*.
 

As SBATCH scripts at: *RNA_seq_Gbi*

```bash
#run_hisat_joinedf.sbatch
indexname="GbV3Hisat"
output="JoinedRNAseqEggHisatV3.sam"
R1="All the fastq files"
R2="All fastqfiles same order R1"
hisat2 -x $indexname  \
        --rna-strandness FR \
        -1 $R1 \
        -2 $R2 \
        -p  $SLURM_NTASKS \
        --dta \
        -S $output

# 805155792 reads; of these:
#   805155792 (100.00%) were paired; of these:
#     346261054 (43.01%) aligned concordantly 0 times
#     372076986 (46.21%) aligned concordantly exactly 1 time
#     86817752 (10.78%) aligned concordantly >1 times
#     ----
#     346261054 pairs aligned concordantly 0 times; of these:
#       14501406 (4.19%) aligned discordantly 1 time
#     ----
#     331759648 pairs aligned 0 times concordantly or discordantly; of these:
#       663519296 mates make up the pairs; of these:
#         543662182 (81.94%) aligned 0 times
#         72201682 (10.88%) aligned exactly 1 time
#         47655432 (7.18%) aligned >1 times
# 66.24% overall alignment rate

#run_hisat_joinedEmbryo.sbatch
indexname="GbV3Hisat"
output="RNA_seq_Gbi/JoinedEmbryoHisat.sam"
R1="All the fastq files"
R2="All fastqfiles same order R1"
hisat2 -x $indexname  \
        -1 $R1 \
        -2 $R2 \
        -p  $SLURM_NTASKS \
        --dta \
        -S $output

# 25653266 reads; of these:
#   25653266 (100.00%) were paired; of these:
#     22835861 (89.02%) aligned concordantly 0 times
#     1901181 (7.41%) aligned concordantly exactly 1 time
#     916224 (3.57%) aligned concordantly >1 times
#     ----
#     22835861 pairs aligned concordantly 0 times; of these:
#       44484 (0.19%) aligned discordantly 1 time
#     ----
#     22791377 pairs aligned 0 times concordantly or discordantly; of these:
#       45582754 mates make up the pairs; of these:
#         44791354 (98.26%) aligned 0 times
#         601101 (1.32%) aligned exactly 1 time
#         190299 (0.42%) aligned >1 times
# 12.70% overall alignment rate



```

 
## ESTs / assembled transcriptomes

  * 1- I will used the previously assembled transcriptome in the lab and located at */Asgard_version/Gryllus.fa*

  * 2- ESTS from NCBI. I retrieevd all the "Nucleotide" sequences form *G. bimaculatus* from ncbi. Search select "Nucelotide" and search  "Gryllus bimaculatus"[porgn:__txid6999] ". I got 43,595 sequences. They include the contigs from a Noji's trancriptome of Nymphs... I downloaded the 43,595 sequences through the NCBI webpage. And Locted at: *"Gbi_NCBI_ESTs/ncbi_ests_2019_noambiguous.fasta"*.

  * 3- I have found a recently published transcriptome assembly of prothoracic ganglion (https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0199070). And accessible from:PRJNA376023. *Gbi_prothoracic_gangalion_transcriptome/GFMG01.1.fsa_nt*
  

* Details on data downloading and processing at *Annotation_v1*
 
## ESTs *Laupala kohalensis* 

Maker accepts ESTs from closely related species that are searched within the genome by tblastx. To that end, I downloaded 14,391  ESTs from *Laupala kohalensis* available at NCBI-nucleotide. Most of these EST come from this paper https://bmcgenomics.biomedcentral.com/articles/10.1186/1471-2164-8-109. *Laupala_kohalensis_annotation/Protein_coding_genes/Ncbi_ESTs/ncbi_Lko_jan2019_noambiguous.fasta *

* Details about processing *Annotation_v1* 

## Proteins

Maker will use protein database to BlastX against the genome and identify protein coding genes. To that end, I will use the Insect rvised uniprot/Swiss-Prot database.

*Insect_swissprot/uniprot-reviewed%3Ayes+taxonomy%3A50557.fasta*

## Prepare First MAKER cycle (v3)

Details about installation at *Annotation_v1* 

```bash
cp -r Annotation_v1/Protein_coding_genes/Maker/Maker_v3/maker/ .

cp -r Annotation_v1/Protein_coding_genes/Maker/Maker_v3/config/ .

## I need to add the AUgustus models created by Busco for Genome V3

cp -r GBI_Genome_v3/Busco_Analysis/Busco_Genome/run_BuscoInsecta_GBIv3/augustus_output/retraining_parameters/ config/species/GBimaculatusV3
# rename files
cd config/species/GBimaculatusV3

 rename "BUSCO_BuscoInsecta_GBIv3_3289539766" "GBimaculatusV3"  *
# #We also need to rename the files cited within certain HMM configuration files.
sed -i 's/BUSCO_BuscoInsecta_GBIv3_3289539766/GBimaculatusV3/g' GBimaculatusV3_parameters.cfg
#sed -i 's/BUSCO_BuscoInsecta_GBIv3_3289539766/GBimaculatusV3/g' GBimaculatusV3_parameters.cfg.orig1
# 
```


On the maker_opts.ctl we put the Augustus species *"GBimaculatusV3"* and we indicate Augustus where to find this organism adding the following line to our sbatch script:

```bash
## Must remember to put the follwoing line in the sbatch file:
export AUGUSTUS_CONFIG_PATH=Maker/config

````
 
Create the Maker config files:
  
```bash
srun -p test --pty --mem 5gb -n 2 -N 1 -t 0-06:00 /bin/bash

echo "load modules"
module load centos6/0.0.1-fasrc01
module load tRNAscan-SE/1.23-fasrc01
module load perl/5.10.1-fasrc01
module load perl-modules/5.10.1-fasrc12
module load GeneMark-ES/4.30-fasrc01
 
echo "load Maker modules on CentOS6"
module load gcc/5.2.0-fasrc01
module load openmpi/2.0.1-fasrc01
module load maker/2.31.8-fasrc01

maker/bin/maker -CTL

```
 
This creates the files:
 
   * maker_exe.ctl - contains the path information for the underlying executables.
   * maker_bopt.ctl - contains filtering statistics for BLAST and Exonerate
   * maker_opts.ctl - contains all other information for MAKER, including the location of the input genome file. 

 For the 1st round of maker I edited the maker_opts.ctl  lines:
 
* Genome
 * genome=Gbimaculatus_Gap_filled.fasta
* ESTEvidence
   * est=/Asgard_version/Gryllus.fa,Gbi_prothoracic_gangalion_transcriptome/GFMG01.1.fsa_nt,Gbi_NCBI_ESTs/ncbi_ests_2019_noambiguous.fasta
  * altest=Laupala_kohalensis_annotation/Protein_coding_genes/Ncbi_ESTs/ncbi_Lko_jan2019_noambiguous.fasta #EST/cDNA sequence file in fasta format from an alternate organism
  * est_gff =RNA_seq_Gbi/stringtie_JoinedRNAseqEggHisat_sorted.gff3,RNA_seq_Gbi/stringtieout_JoinedEmbryoHisat_sorted.gff3
* Protein Homology Evidence
  * protein=Insect_swissprot/uniprot-reviewed%3Ayes+taxonomy%3A50557.fasta
* Repeat Masking section 
  * rmlib=GBI_Genome_v3/Repetitive_content/Custom_Library/CombinedLibrary.lib.minlen50.nr.classified.filtered.fa  #provide an organism specific repeat library in fasta format for RepeatMasker
* Gene Prediction
  * gmhmm=GeneMark/output/gmhmm.mod
  * augustus_species=GBimaculatusV3 #config/species/GBimaculatusV3/
  * est2genome=1 !!!!!! *This option allows you to make gene models directly from the transcript evidence. This option is useful when you don't have a gene predictor trained on your organism and ere is not a training file available for a closely related organism. The gene models from this option are going to be fragmented and incomplete because of the nature of transcript data, especial mRNA-Seq. These gene models are most useful for first round training of gene finders. One you have a trained gene predictor turn this option off. *
  * protein2genome=1 !!!!!! *Similar to est2genome this option will make gene models out of protein data. Like est2genome this option is most useful for training gene predictors. *
  * trna=1 #find tRNAs with tRNAscan, 1 = yes, 0 = no
* External Application Behavior Options
  * cpus = 1 ## USING MPI ELAVE 1
* MAKER Behavior Options
  * single_exon=1 
  * split_hit=25000  # In Bger (30000 bc was large genome), defauklt is 100000
  * tries=6 #number of times to try a contig if there is a failure for some reason

 
The  maker_opts.ctl  file looks like this:
 
```bash
#-----Genome (these are always required)
genome=Gbimaculatus_Gap_filled.fasta #genome sequence (fasta file or fasta embeded in GFF3 file)
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

#-----Re-annotation Using MAKER Derived GFF3
maker_gff= #MAKER derived GFF3 file
est_pass=0 #use ESTs in maker_gff: 1 = yes, 0 = no
altest_pass=0 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
protein_pass=0 #use protein alignments in maker_gff: 1 = yes, 0 = no
rm_pass=0 #use repeats in maker_gff: 1 = yes, 0 = no
model_pass=0 #use gene models in maker_gff: 1 = yes, 0 = no
pred_pass=0 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
other_pass=0 #passthrough anyything else in maker_gff: 1 = yes, 0 = no

#-----EST Evidence (for best results provide a file for at least one)
est=/Asgard_version/Gryllus.fa,Gbi_prothoracic_gangalion_transcriptome/GFMG01.1.fsa_nt,Gbi_NCBI_ESTs/ncbi_ests_2019_noambiguous.fasta #set of ESTs or assembled mRNA-seq in fasta format
altest=Laupala_kohalensis_annotation/Protein_coding_genes/Ncbi_ESTs/ncbi_Lko_jan2019_noambiguous.fasta #EST/cDNA sequence file in fasta format from an alternate organism
est_gff=RNA_seq_Gbi/stringtie_JoinedRNAseqEggHisat_sorted.gff3,RNA_seq_Gbi/stringtieout_JoinedEmbryoHisat_sorted.gff3 #aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=Insect_swissprot/uniprot-reviewed%3Ayes+taxonomy%3A50557.fasta #protein sequence file in fasta format (i.e. from mutiple oransisms)
protein_gff=  #aligned protein homology evidence from an external GFF3 file

#-----Repeat Masking (leave values blank to skip repeat masking)
model_org= #select a model organism for RepBase masking in RepeatMasker
rmlib=GBI_Genome_v3/Repetitive_content/Custom_Library/CombinedLibrary.lib.minlen50.nr.classified.filtered.fa #provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein= #provide a fasta file of transposable element proteins for RepeatRunner
rm_gff= #pre-identified repeat elements from an external GFF3 file
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

#-----Gene Prediction
snaphmm= #SNAP HMM file
gmhmm=GeneMark/output/gmhmm.mod #GeneMark HMM file
augustus_species=GBimaculatusV3 #Augustus gene prediction species model
fgenesh_par_file= #FGENESH parameter file
pred_gff= #ab-initio predictions from an external GFF3 file
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
est2genome=1 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=1 #infer predictions from protein homology, 1 = yes, 0 = no
trna=1 #find tRNAs with tRNAscan, 1 = yes, 0 = no
snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no

#-----Other Annotation Feature Types (features MAKER doesn't recognize)
other_gff= #extra features to pass-through to final MAKER generated GFF3 file

#-----External Application Behavior Options
alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
cpus=1 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)

#-----MAKER Behavior Options
max_dna_len=100000 #length for dividing up contigs into chunks (increases/decreases memory usage)
min_contig=1 #skip genome contigs below this length (under 10kb are often useless)

pred_flank=200 #flank for extending evidence clusters sent to gene predictors
pred_stats=0 #report AED and QI statistics for all predictions as well as models
AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
min_protein=0 #require at least this many amino acids in predicted proteins
alt_splice=0 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
always_complete=0 #extra steps to force start and stop codons, 1 = yes, 0 = no
map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)

split_hit=25000 #length for the splitting of hits (expected max intron size for evidence alignments)
single_exon=1 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
correct_est_fusion=0 #limits use of ESTs in annotation to avoid fusion genes

tries=6 #number of times to try a contig if there is a failure for some reason
clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
TMP= #specify a directory other than the system default temporary directory for temporary files

```


## RUN First MAKER cycle

I will use MPI to make it faster:
 
 
 
```bash
#!/bin/bash
#SBATCH --job-name=Maker_r1
#SBATCH --mail-type=ALL 
#SBATCH -n 200 # Number of cores
#SBATCH -t 7-00:00:00 # Runtime in minutes
#SBATCH -p shared # Partition to submit to
#SBATCH --mem-per-cpu=5Gb # Memory per cpu in MB (see also --mem)
#SBATCH -o out_Maker_r1_%j.out
#SBATCH -e err_Maker_r1_%j.err
pwd; hostname; date

echo "load modules"
module load centos6/0.0.1-fasrc01
module load tRNAscan-SE/1.23-fasrc01
module load perl/5.10.1-fasrc01
module load perl-modules/5.10.1-fasrc12
module load GeneMark-ES/4.30-fasrc01

echo "load Maker modules on CentOS6"
module load gcc/5.2.0-fasrc01
module load openmpi/2.0.1-fasrc01
module load maker/2.31.8-fasrc01

echo "export Augustus config dir"

export AUGUSTUS_CONFIG_PATH=Maker/config

echo "Run Maker"
srun -n $SLURM_NTASKS --mpi=pmi2 maker/bin/maker  
echo "Done!"
pwd; hostname; date

```
With 200 cpus and 5Gb of RAM each CPU, took ~2days. 


To collect the maker output we need to run:

```bash
srun -p test --pty --mem 20gb -n 4 -N 1 -t 0-06:00 /bin/bash

module load centos6/0.0.1-fasrc01
module load gcc/5.2.0-fasrc01
module load openmpi/2.0.1-fasrc01
module load maker/2.31.8-fasrc01


# backup outputs in a higher directory
cd Maker/Maker_GbiV3_v1

mkdir Maker_GbiV3_v1_round1
cd Maker_GbiV3_v1_round1

../maker/bin/fasta_merge -d ../Gbimaculatus_Gap_filled.maker.output/Gbimaculatus_Gap_filled_master_datastore_index.log
../maker/bin/gff3_merge -d  ../Gbimaculatus_Gap_filled.maker.output/Gbimaculatus_Gap_filled_master_datastore_index.log 
../maker/bin/gff3_merge -n -g -o Maker_GbiV3_v1_r1.onlyGenes.gff -d ../Gbimaculatus_Gap_filled.maker.output/Gbimaculatus_Gap_filled_master_datastore_index.log # use -n to not include fasta  + If you want the MAKER gene models only use the -g flag


##count number of protein codying genes 
grep -c -P "\tmaker\tgene"  Gbimaculatus_Gap_filled.all.gff

# ## from Gb_v2.1.all.gff get the different evidences:
# awk '{ if ($2 == "est2genome") print $0 }' Gbimaculatus_Gap_filled.all.gff > Gbimaculatus_Gap_filled.maker.est2genome.gff
# # protein alignments
# awk '{ if ($2 == "protein2genome") print $0 }' Gbimaculatus_Gap_filled.all.gff  >Gbimaculatus_Gap_filled.maker.protein2genome.gff
# # repeat alignments
# awk '{ if ($2 ~ "repeat") print $0 }'  Gbimaculatus_Gap_filled.all.gff > Gbimaculatus_Gap_filled.maker.repeats.gff


## Transform to SNAP format
../maker/bin/maker2zff Gbimaculatus_Gap_filled.all.gff
## creates genome.ann and genome.dna

```


### SNAP round 1

After 1st MAKER cycle, we can train SNAP follwoing MAKER guidleines

 * http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/MAKER_Tutorial_for_WGS_Assembly_and_Annotation_Winter_School_2018#Running_MAKER_with_example_data

Another usefule tutorial:

 * https://gist.github.com/darencard/bb1001ac1532dd4225b030cf0cd61ce2

**MAKER gives the user the option to produce gene annotations directly from the EST evidence. You can then use these imperfect gene models to train gene predictor program. Once you have re-run MAKER with the newly trained gene predictor, you can use the second set of gene annotations to train the gene predictors yet again. This boot-strap process allows you to iteratively improve the performance of ab initio gene predictors**

```bash
srun -p test --pty --mem 20gb -n 4 -N 1 -t 0-06:00 /bin/bash
module load centos6/0.0.1-fasrc01
module load gcc/5.2.0-fasrc01
module load openmpi/2.0.1-fasrc01
module load maker/2.31.8-fasrc01
module load snap/2013.11.29-fasrc01


cd SNAP_training/After_Maker_roun1

## copy the 2 files for Snap
cp Maker/Maker_GbiV3_v1_round1/genome.ann .
cp Maker/Maker_GbiV3_v1_round1/genome.dna .

#The basic steps for training SNAP are first to filter the input gene models, then capture genomic sequence immediately surrounding each model locus, and finally uses those captured segments to produce the HMM. You can explore the internal SNAP documentation for more details if you wish.
 fathom -categorize 1000 genome.ann genome.dna
 fathom -export 1000 -plus uni.ann uni.dna
 forge export.ann export.dna
 hmm-assembler.pl Gbimaculatus_Gap_filled . > Gbimaculatus_Gap_filled.hmm

```

## Run 2nd MAKER cycle

I backup the *maker_opts.ctl* and I edit this one for the 2nd maker round.

```bash
cp maker_opts.ctl maker_optsctl_round1

```

After havintg SNAP trained with the 1st Maker iteration models, I run maker again in the same directory, *only changeing snaphmm, est2genome and protein2genome *.

*I have previusly done some checks, and running in differetnt directory, takes much longer and produced worse results.*


* Gene Prediction
  * snaphmm=SNAP_training/After_Maker_roun1/Gbimaculatus_Gap_filled.hmm
  * est2genome=0
  * protein2genome=0
  * alt_splice=1


And the maker_opts.ctl file should look like:

```bash
#-----Genome (these are always required)
genome=Gbimaculatus_Gap_filled.fasta #genome sequence (fasta file or fasta embeded in GFF3 file)
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

#-----Re-annotation Using MAKER Derived GFF3
maker_gff= #MAKER derived GFF3 file
est_pass=0 #use ESTs in maker_gff: 1 = yes, 0 = no
altest_pass=0 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
protein_pass=0 #use protein alignments in maker_gff: 1 = yes, 0 = no
rm_pass=0 #use repeats in maker_gff: 1 = yes, 0 = no
model_pass=0 #use gene models in maker_gff: 1 = yes, 0 = no
pred_pass=0 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
other_pass=0 #passthrough anyything else in maker_gff: 1 = yes, 0 = no

#-----EST Evidence (for best results provide a file for at least one)
est=/Asgard_version/Gryllus.fa,Gbi_prothoracic_gangalion_transcriptome/GFMG01.1.fsa_nt,Gbi_NCBI_ESTs/ncbi_ests_2019_noambiguous.fasta #set of ESTs or assembled mRNA-seq in fasta format
altest=Laupala_kohalensis_annotation/Protein_coding_genes/Ncbi_ESTs/ncbi_Lko_jan2019_noambiguous.fasta #EST/cDNA sequence file in fasta format from an alternate organism
est_gff=RNA_seq_Gbi/stringtie_JoinedRNAseqEggHisat_sorted.gff3,RNA_seq_Gbi/stringtieout_JoinedEmbryoHisat_sorted.gff3 #aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=Insect_swissprot/uniprot-reviewed%3Ayes+taxonomy%3A50557.fasta #protein sequence file in fasta format (i.e. from mutiple oransisms)
protein_gff=  #aligned protein homology evidence from an external GFF3 file

#-----Repeat Masking (leave values blank to skip repeat masking)
model_org= #select a model organism for RepBase masking in RepeatMasker
rmlib=GBI_Genome_v3/Repetitive_content/Custom_Library/CombinedLibrary.lib.minlen50.nr.classified.filtered.fa #provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein= #provide a fasta file of transposable element proteins for RepeatRunner
rm_gff= #pre-identified repeat elements from an external GFF3 file
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

#-----Gene Prediction
snaphmm=SNAP_training/After_Maker_roun1/Gbimaculatus_Gap_filled.hmm #SNAP HMM file
gmhmm=GeneMark/output/gmhmm.mod #GeneMark HMM file
augustus_species=GBimaculatusV3 #Augustus gene prediction species model
fgenesh_par_file= #FGENESH parameter file
pred_gff= #ab-initio predictions from an external GFF3 file
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
est2genome=0 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=0 #infer predictions from protein homology, 1 = yes, 0 = no
trna=1 #find tRNAs with tRNAscan, 1 = yes, 0 = no
snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no

#-----Other Annotation Feature Types (features MAKER doesn't recognize)
other_gff= #extra features to pass-through to final MAKER generated GFF3 file

#-----External Application Behavior Options
alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
cpus=1 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)

#-----MAKER Behavior Options
max_dna_len=100000 #length for dividing up contigs into chunks (increases/decreases memory usage)
min_contig=1 #skip genome contigs below this length (under 10kb are often useless)

pred_flank=200 #flank for extending evidence clusters sent to gene predictors
pred_stats=0 #report AED and QI statistics for all predictions as well as models
AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
min_protein=0 #require at least this many amino acids in predicted proteins
alt_splice=1 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
always_complete=0 #extra steps to force start and stop codons, 1 = yes, 0 = no
map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)

split_hit=25000 #length for the splitting of hits (expected max intron size for evidence alignments)
single_exon=1 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
correct_est_fusion=0 #limits use of ESTs in annotation to avoid fusion genes

tries=10 #number of times to try a contig if there is a failure for some reason
clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
TMP= #specify a directory other than the system default temporary directory for temporary files
```


 Now, we can launch a sbatch job for running Maker (run_maker_r2.sbatch)
 
```bash
#!/bin/bash
#SBATCH --job-name=Maker_r2
#SBATCH --mail-type=ALL 
#SBATCH -n 50 # Number of cores
#SBATCH -t 7-00:00:00 # Runtime in minutes
#SBATCH -p general # Partition to submit to
#SBATCH --mem-per-cpu=5Gb # Memory per cpu in MB (see also --mem)
#SBATCH -o out_Maker_r2_%j.out
#SBATCH -e err_Maker_r2_%j.err
pwd; hostname; date

echo "load modules"
module load centos6/0.0.1-fasrc01
module load tRNAscan-SE/1.23-fasrc01
module load perl/5.10.1-fasrc01
module load perl-modules/5.10.1-fasrc12
module load GeneMark-ES/4.30-fasrc01

echo "load Maker modules on CentOS6"
module load gcc/5.2.0-fasrc01
module load openmpi/2.0.1-fasrc01
module load maker/2.31.8-fasrc01

echo "export Augustus config dir"

export AUGUSTUS_CONFIG_PATH=Maker/config

echo "Run Maker"
srun -n $SLURM_NTASKS --mpi=pmi2 maker/bin/maker  
echo "Done!"
pwd; hostname; date


```
Collect maker round 2 output  as sbatch script.

```bash

# backup outputs in a higher directory
cd Maker/Maker_GbiV3_v1

#mkdir Maker_GbiV3_v1_round2
cd Maker_GbiV3_v1_round2

../maker/bin/fasta_merge -d ../Gbimaculatus_Gap_filled.maker.output/Gbimaculatus_Gap_filled_master_datastore_index.log
../maker/bin/gff3_merge -d  ../Gbimaculatus_Gap_filled.maker.output/Gbimaculatus_Gap_filled_master_datastore_index.log 
../maker/bin/gff3_merge -n -g -o Maker_GbiV3_v1_r2.onlyGenes.gff -d ../Gbimaculatus_Gap_filled.maker.output/Gbimaculatus_Gap_filled_master_datastore_index.log # use -n to not include fasta  + If you want the MAKER gene models only use the -g flag

## Transform to SNAP format
../maker/bin/maker2zff Gbimaculatus_Gap_filled.all.gff
## creates genome.ann and genome.dna

```

### SNAP round 2

After 2nd MAKER cycle, we can train SNAP

```bash
srun -p test --pty --mem 20gb -n 4 -N 1 -t 0-06:00 /bin/bash
module load centos6/0.0.1-fasrc01
module load gcc/5.2.0-fasrc01
module load openmpi/2.0.1-fasrc01
module load maker/2.31.8-fasrc01
module load snap/2013.11.29-fasrc01


cd SNAP_training/After_Maker_round2

## copy the 2 files for Snap
# cp Maker/Maker_GbiV3_v1_round2/genome.ann .
# cp Maker/Maker_GbiV3_v1_round2/genome.dna .

#The basic steps for training SNAP are first to filter the input gene models, then capture genomic sequence immediately surrounding each model locus, and finally uses those captured segments to produce the HMM. You can explore the internal SNAP documentation for more details if you wish.
 fathom -categorize 1000 genome.ann genome.dna
 fathom -export 1000 -plus uni.ann uni.dna
 forge export.ann export.dna
 hmm-assembler.pl Gbimaculatus_Gap_filled . > Gbimaculatus_Gap_filled.hmm

```


## Run 3rd MAKER cycle

I backup the *maker_opts.ctl* and I edit this one for the 2nd maker round.

```bash
cp maker_opts.ctl maker_optsctl_round2

```

After havintg SNAP trained with the 2nd Maker iteration models, I run maker again in the same directory, only changing:

* Gene Prediction
  * snaphmm=SNAP_training/After_Maker_round2/Gbimaculatus_Gap_filled.hmm



And the maker_opts.ctl file should look like:

```bash
#-----Genome (these are always required)
genome=Gbimaculatus_Gap_filled.fasta #genome sequence (fasta file or fasta embeded in GFF3 file)
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

#-----Re-annotation Using MAKER Derived GFF3
maker_gff= #MAKER derived GFF3 file
est_pass=0 #use ESTs in maker_gff: 1 = yes, 0 = no
altest_pass=0 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
protein_pass=0 #use protein alignments in maker_gff: 1 = yes, 0 = no
rm_pass=0 #use repeats in maker_gff: 1 = yes, 0 = no
model_pass=0 #use gene models in maker_gff: 1 = yes, 0 = no
pred_pass=0 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
other_pass=0 #passthrough anyything else in maker_gff: 1 = yes, 0 = no

#-----EST Evidence (for best results provide a file for at least one)
est=/Asgard_version/Gryllus.fa,Gbi_prothoracic_gangalion_transcriptome/GFMG01.1.fsa_nt,Gbi_NCBI_ESTs/ncbi_ests_2019_noambiguous.fasta #set of ESTs or assembled mRNA-seq in fasta format
altest=Laupala_kohalensis_annotation/Protein_coding_genes/Ncbi_ESTs/ncbi_Lko_jan2019_noambiguous.fasta #EST/cDNA sequence file in fasta format from an alternate organism
est_gff=RNA_seq_Gbi/stringtie_JoinedRNAseqEggHisat_sorted.gff3,RNA_seq_Gbi/stringtieout_JoinedEmbryoHisat_sorted.gff3 #aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=Insect_swissprot/uniprot-reviewed%3Ayes+taxonomy%3A50557.fasta #protein sequence file in fasta format (i.e. from mutiple oransisms)
protein_gff=  #aligned protein homology evidence from an external GFF3 file

#-----Repeat Masking (leave values blank to skip repeat masking)
model_org= #select a model organism for RepBase masking in RepeatMasker
rmlib=GBI_Genome_v3/Repetitive_content/Custom_Library/CombinedLibrary.lib.minlen50.nr.classified.filtered.fa #provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein= #provide a fasta file of transposable element proteins for RepeatRunner
rm_gff= #pre-identified repeat elements from an external GFF3 file
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

#-----Gene Prediction
snaphmm=SNAP_training/After_Maker_round2/Gbimaculatus_Gap_filled.hmm #SNAP HMM file
gmhmm=GeneMark/output/gmhmm.mod #GeneMark HMM file
augustus_species=GBimaculatusV3 #Augustus gene prediction species model
fgenesh_par_file= #FGENESH parameter file
pred_gff= #ab-initio predictions from an external GFF3 file
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
est2genome=0 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=0 #infer predictions from protein homology, 1 = yes, 0 = no
trna=1 #find tRNAs with tRNAscan, 1 = yes, 0 = no
snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no

#-----Other Annotation Feature Types (features MAKER doesn't recognize)
other_gff= #extra features to pass-through to final MAKER generated GFF3 file

#-----External Application Behavior Options
alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
cpus=1 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)

#-----MAKER Behavior Options
max_dna_len=100000 #length for dividing up contigs into chunks (increases/decreases memory usage)
min_contig=1 #skip genome contigs below this length (under 10kb are often useless)

pred_flank=200 #flank for extending evidence clusters sent to gene predictors
pred_stats=0 #report AED and QI statistics for all predictions as well as models
AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
min_protein=0 #require at least this many amino acids in predicted proteins
alt_splice=1 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
always_complete=0 #extra steps to force start and stop codons, 1 = yes, 0 = no
map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)

split_hit=25000 #length for the splitting of hits (expected max intron size for evidence alignments)
single_exon=1 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
correct_est_fusion=0 #limits use of ESTs in annotation to avoid fusion genes

tries=10 #number of times to try a contig if there is a failure for some reason
clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
TMP= #specify a directory other than the system default temporary directory for temporary files

```


Collect maker round 3 output:

```bash
## No scaffold skipped by maker!
grep  "DIED_SKIPPED_PERMANENT" Gbimaculatus_Gap_filled.maker.output/Gbimaculatus_Gap_filled_master_datastore_index.log

# backup outputs in a higher directory
cd Maker/Maker_GbiV3_v1

#mkdir Maker_GbiV3_v1_round3
cd Maker_GbiV3_v1_round3

../maker/bin/fasta_merge -d ../Gbimaculatus_Gap_filled.maker.output/Gbimaculatus_Gap_filled_master_datastore_index.log
../maker/bin/gff3_merge -d  ../Gbimaculatus_Gap_filled.maker.output/Gbimaculatus_Gap_filled_master_datastore_index.log 
../maker/bin/gff3_merge -n -g -o Maker_GbiV3_v1_r2.onlyGenes.gff -d ../Gbimaculatus_Gap_filled.maker.output/Gbimaculatus_Gap_filled_master_datastore_index.log # use -n to not include fasta  + If you want the MAKER gene models only use the -g flag
```


## Annotations quality after each round:

I can use BUSCO to check how complete are my annotations.

Quality of the transcriptome after each Maker round (version 3).


```bash
#### Round 1
cd Annotations_Quality/Maker_GbiV3_v1_round1_completed

Fasta="Maker/Maker_GbiV3_v1_round1/Gbimaculatus_Gap_filled.all.maker.transcripts.fasta"
outfilename="Busco_maker_v1r1"
python ~/busco/scripts/run_BUSCO.py --cpu $SLURM_NTASKS -i $Fasta  -o $outfilename -l ~/busco/datasets/insecta_odb9 -m tran

#### Round 2
cd Annotations_Quality/Maker_GbiV3_v1_round2_completed

Fasta="Maker/Maker_GbiV3_v1_round2/Gbimaculatus_Gap_filled.all.maker.transcripts.fasta"
outfilename="Busco_maker_v1r2"
python ~/busco/scripts/run_BUSCO.py --cpu $SLURM_NTASKS -i $Fasta  -o $outfilename -l ~/busco/datasets/insecta_odb9 -m tran

#### Round 3
cd Annotations_Quality/Maker_GbiV3_v1_round3_completed

Fasta="Maker/Maker_GbiV3_v1_round3/Gbimaculatus_Gap_filled.all.maker.transcripts.fasta"
outfilename="Busco_maker_v1r2"
python ~/busco/scripts/run_BUSCO.py --cpu $SLURM_NTASKS -i $Fasta  -o $outfilename -l ~/busco/datasets/insecta_odb9 -m tran

### comapre

## Buscos missing in Round 2 but present in Round 1
 grep -v -f ../Maker_GbiV3_v1_round1_completed/run_Busco_maker_v1r1/missing_busco_list_Busco_maker_v1r1.tsv  run_Busco_maker_v1r2/missing_busco_list_Busco_maker_v1r2.tsv

```
