#!/usr/bin/env nextflow
nextflow.enable.dsl=2
/*
========================================================================================
                         Nextflow-ORganismPlacER(ORPER)
========================================================================================
GIT url : https://github.com/Lcornet/ORPER
----------------------------------------------------------------------------------------
*/

// Show help message and exit the workflow (--help option)
params.help = null
if (params.help){
    helpMessage()
    exit 0
}

// mandatory arguments
params.reftaxon = null
if (params.reftaxon == null) {
    exit 1, "Missing mandatory --reftaxon argument."
}

params.outtaxon = null
if (params.outtaxon == null) {
    exit 1, "Missing mandatory --outtaxon argument."
}

params.reflevel = null
if (params.reflevel == null) {
    exit 1, "Missing mandatory --reflevel argument."
}

params.outlevel = null
if (params.outlevel == null) {
    exit 1, "Missing mandatory --outlevel argument."
}

params.rRnaFile = null
if (params.rRnaFile == null) {
    exit 1, "Missing mandatory --rRnaFile argument."
}

params.kingdom == null
if (params.kingdom == null) {
    exit 1, "Missing mandatory --kingdom argument."
}

// optional arguments
params.ribodb  = 'local'
params.taxdump = 'local'
params.eukccdb = 'local'
params.refgenbank = 'off'
params.outgenbank = 'off'
params.cpu = 1
params.rrna = 'auto'
params.dRep = 'off'
params.cdhit = 'off'
params.checkrrna = 'off'
params.shrink = '0.1'
params.outdir = "$workflow.projectDir/ORPER-results"

// load modules
include { SetupTaxonomy } from './modules/ncbi.nf'
include { SetupRiboDB } from './modules/ribodb.nf'
include { SetUpEukccDB } from './modules/eukccdb.nf'
include { DerepRefGenomes } from './modules/utils.nf'
include { ProteomeDownloader as RefProteomeDownloader; ProteomeDownloader as OutProteomeDownloader} from './modules/ncbi.nf'
include { RiboDBFortytwo } from './modules/fortytwo.nf'
include { Mafft } from './modules/mafft.nf'
include { Scafos } from './modules/scafos.nf'
include { CdHit } from './modules/cdhit.nf'
include { RaxmlReferenceTree; RaxmlConstrainTree } from './modules/phylogeny.nf'

// load workflows
include { GetGenomes as GetRefGenomes; GetGenomes as GetOutGenomes } from './workflows/genome.nf'

workflow {

    // set up ribo database
    SetupRiboDB(params.ribodb)

    // set up taxdump
    SetupTaxonomy(params.taxdump)

    // set up eukcc database
    SetUpEukccDB(params.eukccdb)

    // launch main pipeline, which download genomes
    // assess contamination and predict rRNA
    // for reference genomes
    GetRefGenomes(SetupTaxonomy.out, SetUpEukccDB.out, params.reftaxon, params.reflevel, params.refgenbank)
    // ... and outgroup genomes
    GetOutGenomes(SetupTaxonomy.out, SetUpEukccDB.out, params.outtaxon, params.outlevel, params.outgenbank)

    // combien rRNA prediction for reference and out genomes
    all_rRna = GetRefGenomes.out.rrna_files.combine(GetOutGenomes.out.rrna_files)

    // dereplication of reference genomes if specify by user
    if (params.dRep == 'on') {
        DerepRefGenomes(GetRefGenomes.out.genome_files, GetRefGenomes.out.filtered_gca)
        // download corresponding proteome for reference genomes
        RefProteomeDownloader(GetRefGenomes.out.ftp_files, DerepRefGenomes.out.derep_gca, 'ref')
    } else {
        RefProteomeDownloader(GetRefGenomes.out.ftp_files, GetRefGenomes.out.filtered_gca, 'ref')
    }

    // ... and outgroup genomes
    OutProteomeDownloader(GetOutGenomes.out.ftp_files, GetOutGenomes.out.filtered_gca, 'out')

    // combien reference and outgroup proteome
    all_proteome = RefProteomeDownloader.out.abbr_files.combine(OutProteomeDownloader.out.abbr_files)

    // Reference phylogenomic tree
    //RiboDBFortytwo(RefProteomeDownloader.out.abbr_files, SetupRiboDB.out, SetupTaxonomy.out)
    RiboDBFortytwo(all_proteome, SetupRiboDB.out, SetupTaxonomy.out)
    Mafft(RiboDBFortytwo.out.enrich_files)
    // TODO: Manage empty files in Scafos
    Scafos(Mafft.out.AlignFastaFiles)
    RaxmlReferenceTree(Scafos.out.ConcatFastaFile)

    if (params.cdhit == 'on') {
        CdHit(Channel.fromPath(params.rRnaFile), 1.0)
        RaxmlConstrainTree(RaxmlReferenceTree.out.Tree, RaxmlReferenceTree.out.IdList, all_rRna,
        CdHit.out.ClstFasta, GetOutGenomes.out.filtered_gca, SetupTaxonomy.out)
    } else {
        RaxmlConstrainTree(RaxmlReferenceTree.out.Tree, RaxmlReferenceTree.out.IdList, all_rRna,
        Channel.fromPath(params.rRnaFile), GetOutGenomes.out.filtered_gca, SetupTaxonomy.out)
    }

}

def helpMessage() {
    log.info """

    Description:
    ORPER performs a phylogenetic placement of rRNA sequences in a tree,
    composed of RefSeq genomes, and constrained by a ribosomal phylogenomic tree.

    Citation:
    Please cite:
    Cornet, L., Ahn, A.-C., Wilmotte, A., & Baurain, D. (2021). ORPER: A Workflow for Constrained
    SSU rRNA Phylogenies. Genes, 12 (11). doi:10.3390/genes12111741

    Usage:
    The typical command for running the pipeline is as follows:
    nextflow ORPER.nf --reftaxon=Cyanobacteriota --reflevel=phylum --outtaxon=Melainabacteria
    --outlevel=phylum --cpu=45 --kingdom=Bacteria --refgenbank=off --outgenbank=off --dRep=on
    --cdhit=on --rRnaFile=SequencesULC4Luc-NOS.fasta

    Mandatory arguments:
    --reftaxon          NCBI's taxonomic name of reference group.
    --outtaxon          NCBI's taxonomic name of outgroup.
    --reflevel          Taxonomic level of reference group (e.g., phylum, class, order, family).
    --outlevel          Taxonomic level of outgroup (e.g., phylum, class, order, family).
    --rRnaFile          Path to FASTA file containing rRNA sequences to be placed in reference tree.
    --kingdom           NCBI's kingdom taxon of studied organisms. The two available kingdom are Bacteria and Eukaryota.

    Optional arguments:
    --cpu               Number of threads to run in parallel [default: 1].
    --taxdump           Path to taxdump directory [default: automatic set-up].
    --ribodb            Path to directory containing ribodb fasta files [default: automatic set-up].
    --eukccdb           Path to EukCC database directory [default: automatic set-up].
    --refgenbank        GenBank genomes of reference group are considered for analyses [default: off].
    --outgenbank        GenBank genomes of outgroup are considered for analyses [default: off].
    --rrna              rRNA type to be considered for analyses [default: autodetection]. Avaible rRNA are:
                        16S (default), 23S and 5S for Bacteria and 18S, 28S (default) and 5.8S for Eukaryota.
                        ITS region (ITS1-5.8S-ITS2) is only available for Fungi.
    --checkrrna         Run a validation check when rRNA 28S predictions [default: off].
    --dRep              Run dRep to dereplicate reference group genomes [default: off].
    --cdhit             Run cd-hit on provided rRNA files [default: off].
    --shrink            TreeShrink cutoff value [default: 0.1].
    --outdir            Specify the output directory [default: ./ORPER-results].

    """.stripIndent()
    
}
