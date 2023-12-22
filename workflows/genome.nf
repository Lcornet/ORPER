// load modules
include { GetURL; GetURL as GbGetURL; GenomeDownloader;
    GenomeDownloader as GbGenomeDownloader } from '../modules/ncbi.nf'
include { CheckM; EukCC } from '../modules/contamination.nf'
include { CheckRefSeq; CheckGenBank; ReduceGCA;
    GetUncontaminated; SeqSelector; ChangeIdsAli } from '../modules/utils.nf'
include { Barrnap; FilterBarrnap } from '../modules/barrnap.nf'
include { rRNAFortyTwo; ITSFortyTwo } from '../modules/fortytwo.nf'
include { CdHit } from '../modules/cdhit.nf'

workflow GetGenomes {

    take:
        Taxonomy
        EukccDB
        Taxon
        Level
        GenBank

    main:
        // Download RefSeq Genomes
        RefSeqftp = file('ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt')
        GetURL(Taxonomy, RefSeqftp,Taxon, Level)
        GenomeDownloader(GetURL.out.ftp_links, GetURL.out.gca_list)

        // Check if RefSeq Genomes are downloaded
        CheckRefSeq(GenomeDownloader.out.abbr_files)


        // Download GenBank Genomes
        //if (params.genbank == "on") {
        if (GenBank == 'on') {
            GenBankftp = file('ftp://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt')
            GbGetURL(Taxonomy, GenBankftp, Taxon, Level)
            ReduceGCA(GetURL.out.gca_list, GbGetURL.out.gca_list)
            GbGenomeDownloader(GbGetURL.out.ftp_links, ReduceGCA.out.reduced_gca_list)
            CheckGenBank(GenomeDownloader.out.abbr_files, GbGenomeDownloader.out.abbr_files, Taxon)
        }

        // Combine FTP links files for proteomes downloading
        //all_ftp = params.genbank == 'on' ? GetURL.out.ftp_links.combine(GbGetURL.out.ftp_links)
        all_ftp = GenBank == 'on' ? GetURL.out.ftp_links.combine(GbGetURL.out.ftp_links)
            : GetURL.out.ftp_links

        // Combine RefSeq and GenBank genomes (if users specify --genbank=on)
        //all_genomes = params.genbank == 'on' ? GenomeDownloader.out.abbr_files.combine(GbGenomeDownloader.out.abbr_files)
        all_genomes = GenBank == 'on' ? GenomeDownloader.out.abbr_files.combine(GbGenomeDownloader.out.abbr_files)
            : GenomeDownloader.out.abbr_files

        // Assess genome contamination for procaryotes
        if (params.kingdom == 'Bacteria') {
            CheckM(all_genomes)
        }
        // ... and eukaryotes
        if (params.kingdom == 'Eukaryota') {
            EukCC(all_genomes, EukccDB)
        }

        // If user specify ITS...
        if (params.rrna == 'ITS' && params.kingdom == 'Eukaryota') {
            ITSFortyTwo(all_genomes, Taxonomy)
        }
        else if (params.rrna == 'ITS' && params.kingdom != 'Eukaryota') {
            exit 1, "[ORPER-WARNING] ITS analysis is only available for Eukaryota - Fungi -"
        } else {
            // ... otherwise, it predict rRNA
            Barrnap(all_genomes)
            // ... and get rRNA type specified by user
            rrna = ( params.rrna == 'auto' && params.kingdom == 'Bacteria' ) ? '16S'
                : ( ( params.rrna == 'auto' && params.kingdom == 'Eukaryota' ) ? '28S'
                : params.rrna )
            FilterBarrnap(Barrnap.out.result, rrna)
        }

        // user can assess barrnap's rRNA predictions
        if (params.checkrrna == 'on' && params.kingdom == 'Eukaryota' && params.rrna != 'ITS') {
            rRNAFortyTwo(FilterBarrnap.out.filter_files, Taxonomy)
            CdHit(rRNAFortyTwo.out.fasta_files, 0.995)
            // Perl script that select only 1 rRNA sequence per organism
            SeqSelector(CdHit.out.ClstFasta)
        } else {
            // Same here...
            // TODO: Check if no bug when change-ids-ali with params.rrna = ITS
            params.rrna == 'ITS' ? CdHit(ITSFortyTwo.out.fasta_files, 0.995) : CdHit(FilterBarrnap.out.filter_files, 0.995)
            // adding org name to IDS (because no 42 step here)
            ChangeIdsAli(CdHit.out.ClstFasta, Taxonomy)
            SeqSelector(ChangeIdsAli.out.Outfiles)
        }

        // filter uncontaminated genomes with predicted rRNA
        contam_report = params.kingdom == 'Bacteria' ? CheckM.out.result : EukCC.out.result
        GetUncontaminated(contam_report, SeqSelector.out.gca_list)

    emit:
        ftp_files = all_ftp
        genome_files = all_genomes
        rrna_files = SeqSelector.out.UniqFasta
        filtered_gca = GetUncontaminated.out.reliable

}
