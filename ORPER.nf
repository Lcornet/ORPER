#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// parameters
params.genbank = 'off'
params.kingdom = 'Bacteria'
params.ribodb  = 'local'
params.taxdump = 'local'
params.eukccdb = 'local'
params.cpu = 1
params.rrna = 'auto'
params.dRep = 'off'
params.checkrrna = 'off'
params.shrink = '0.1'

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
    GetRefGenomes(SetupTaxonomy.out, SetUpEukccDB.out, params.reftaxon, params.reflevel)
    // ... and outgroup genomes
    GetOutGenomes(SetupTaxonomy.out, SetUpEukccDB.out, params.outtaxon, params.outlevel)

    // combien rRNA prediction for ea
    rRnaFiles = GetRefGenomes.out.rrna_files.combine(GetOutGenomes.out.rrna_files)

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
    }

    RaxmlConstrainTree(RaxmlReferenceTree.out.Tree, RaxmlReferenceTree.out.IdList, rRnaFiles,
        CdHit.out.ClstFasta, GetOutGenomes.out.filtered_gca, SetupTaxonomy.out)
    //RaxmlConstrainTree(RaxmlReferenceTree.out.Tree, RaxmlReferenceTree.out.IdList, GetRefGenomes.out.rrna_files,
    //    CdHit.out.ClstFasta, GetOutGenomes.out.filtered_gca, SetupTaxonomy.out)
    
}
