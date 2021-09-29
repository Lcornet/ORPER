# ORganismPlacER (ORPER)

This repository contains singularity definition file and Nextflow script for ORPER.  

## Install

### singularity
Please first install singularity, follow:
https://stackoverflow.com/questions/63865962/singularity-3-6-2-installation  

### ORPER

Please build the singularity container using ORPER.def and make sure to place RNAmmer in the ORPER directory. Please visit https://services.healthtech.dtu.dk/software.php to download RNAmmer v1.2.

    git clone https://github.com/Lcornet/ORPER
    sudo singularity build ORPER.sif ORPER.def

## Usage

### Quick Usage

Classical command line to work on a phylum, here Cyanobacteria.

    nextflow ORPER.nf --reftaxolevel=phylum --refgroup=Cyanobacteria --outtaxolevel=phylum --outgroup=Melainabacteria --outgenbank=yes --cpu=30 --SSU=SequencesULC4Luc.fasta -with-singularity ORPER.sif

This command should give a nextflow report like this:

    executor >  local (23)
    [3f/8bd288] process > RiboDBSetUp (1)             [100%] 1 of 1 ✔
    [1b/d33ae5] process > Taxonomy (1)                [100%] 1 of 1 ✔
    [0b/e20645] process > RefSeq (1)                  [100%] 1 of 1 ✔
    [ff/04aef1] process > GenBank (1)                 [100%] 1 of 1 ✔
    [b8/c9dc0a] process > GetRefGenomesRefseq (1)     [100%] 1 of 1 ✔
    [71/94add5] process > GetRefGenomesGenbank (1)    [100%] 1 of 1 ✔
    [40/965983] process > RefGenomesCheckm (1)        [100%] 1 of 1 ✔
    [aa/5f3aa2] process > RefGenomesBarnap (1)        [100%] 1 of 1 ✔
    [23/4dcb16] process > RefGenomesFilter (1)        [100%] 1 of 1 ✔
    [1a/059ff6] process > RefGenomesDereplication (1) [100%] 1 of 1 ✔
    [65/b930c7] process > GetRefRelProteomes (1)      [100%] 1 of 1 ✔
    [d9/2f9308] process > GetOutGenomesRefSeq (1)     [100%] 1 of 1 ✔
    [2b/666443] process > GetOutGenomesGenbank (1)    [100%] 1 of 1 ✔
    [40/269a15] process > OutGenomesCheckm (1)        [100%] 1 of 1 ✔
    [4d/cb7dbd] process > OutGenomesBarnap (1)        [100%] 1 of 1 ✔
    [c9/1792af] process > OutGenomesFilter (1)        [100%] 1 of 1 ✔
    [f5/a63730] process > GetOutRelProteomes (1)      [100%] 1 of 1 ✔
    [01/69d4d9] process > RiboDBFortytwo (1)          [100%] 1 of 1 ✔
    [91/8f3231] process > AlignmentMUSCLE (1)         [100%] 1 of 1 ✔
    [79/153441] process > ConcatScafos (1)            [100%] 1 of 1 ✔
    [d5/6d85e1] process > ReferenceTreeRaxml (1)      [100%] 1 of 1 ✔
    [84/6240bd] process > SSUDereplication (1)        [100%] 1 of 1 ✔
    [2a/23433b] process > ConstrainTreeRaxml (1)      [100%] 1 of 1 ✔
    [1e/1edc61] process > PublicationResults (1)      [100%] 1 of 1 ✔
    Completed at: 07-Sep-2021 09:22:27
    Duration    : 20h 47m 27s
    CPU hours   : 25.6
    Succeeded   : 23

### Advance usage

Please follow the help of ORPER:

    Mandatory arguments:
    --refgroup          Group of interest
    --outgroup          Outgroup
    --reftaxolevel      Taxonomic level of reference group - Choice between four taxa levels: phylum, class, order, family
    --outtaxolevel      Taxonomic level of outgroup - Choice between four taxa levels: phylum, class, order, family
    --SSU               Path to fasta file containing SSU sequences                     

    Optional arguments:
    --ribodb            Path to directory containing ribodb fasta files - automatic setup by default
    --companion         Path to ORPER-companion.py 
    --taxdump 		Path to taxdump directory - automatic setup by default
    --cpu               Number of cpus, default = 1
    --genbank           Download GenBank metadata, activated by default - yes or no. Needed if refgenbank or outgenbank is activated
    --refgenbank        add GenBank for reference group, Deactivated by default - yes or no
    --outgenbank        add GenBank for outgroup, Deactivated by default - yes or no
    --dRep              dRep dereplication for group of interest, activated by default - yes or no
    --cdhit             cdhit dereplication of provided SSU sequences, activated by default - yes or no
    --shrink            TreeShrink cutoff value, 0.1 by default

### HPC usage

For HPC usage, we recommend to use Nextflow configuration file: nextflow.config.  
This file, which as to be in the folder where you ORPER should look like this:

	process.container = '/scratch/ulg/bioec/lcornet/ORPER/ORPER.sif'
	singularity.enabled = true
	singularity.cacheDir = "$PWD"
	singularity.autoMounts = false
	singularity.runOptions = '-B <Path-to-ORPER-working-dir> -B <Path-to-symlink-dir>'

Please note that if your HPC system use symbolic link to mount partition, both links (hard link and symlink) has to be added
in the bind (-B) linl of the config file.  
