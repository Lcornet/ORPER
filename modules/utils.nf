process ReduceGCA {
    input:
	file gcf_list
        file gca_list

    output:
        path 'reduced_gca.list', emit: reduced_gca_list


    script:
        """
        if [ ! -s $gcf_list ]; then
            cp $gca_list reduced_gca.list
        elif [ ! -s $gca_list ]; then
            touch reduced_gca.list
        elif [ -s $gcf_list ] && [ -s $gca_list ]; then
            perl -ple 's/GCF/GCA/' $gcf_list > refseq2genbank.gca
            if cmp -s refseq2genbank.gca $gca_list; then
                touch reduced_gca.list
            else
                grep -vf refseq2genbank.gca $gca_list > reduced_gca.list
            fi
        fi
        """
}

process CheckRefSeq {
    input:
        file refseq_genomes

    script:
        files = refseq_genomes.name
        if ('refseq_empty-abbr.fna' in files && (params.refgenbank == 'off' || params.outgenbank == 'off')) {
            exit 1, "[ORPER-WARNING] No RefSeq genomes have been downloaded. Rerun ORPER with `--genbank` option enable"
        }

        else if ('refseq_empty-abbr.fna' in files && (params.refgenbank == 'on' || params.outgenbank == 'on')) {
            println("[ORPER-INFO] No RefSeq genomes have been downloaded. Continuing with GenBank genomes only")
            """
            touch checkrefseq.log
            """
        }

        else {
            """
            touch checkrefseq.log
            """
        }
}

process CheckGenBank {
    input:
        file refseq_genomes
        file genbank_genomes
        val taxon

    script:
        rs_files = refseq_genomes.name
        gb_files = genbank_genomes.name
        if ( 'refseq_empty-abbr.fna' in rs_files && 'genbank_empty-abbr.fna' in gb_files ) {
            exit 1, "[ORPER-WARNING] No genomes downloaded with taxon: ${taxon}. Try with another taxon name"
        }

	else if ( 'refseq_empty-abbr.fna' in rs_files && !('genbank_empty-abbr.fna' in gb_files) ) {
            println("[ORPER-INFO] No genomes with taxon: ${taxon} in RefSeq database. Using GenBank only for the following analyses")
            """
            touch checkgenbank.log
            """
        }

        else if ( 'genbank_empty-abbr.fna' in gb_files && !('refseq_empty-abbr.fna' in rs_files) ) {
            println("[ORPER-INFO] Redundancy between GenBank and RefSeq database. Only RefSeq genomes are kept for the following analyses")
            """
            touch checkgenbank.log
            """
        }

	else {
            """
            touch checkgenbank.log
            """
        }
}

process GetUncontaminated {
    input:
	file contam
        file barrnap

    output:
        path 'reliable_gca.list', emit: reliable

    script:
        compl = contam.name == 'checkm_result.csv' ? 90 : 85

        """
        perl -F"," -anle 'm/(GC[AF]_\\d{9}\\.\\d+)/ ; print \$1 if \$F[1] > $compl && \$F[2] < 5' $contam > noncontam_gca.list
        comm -12 <(sort noncontam_gca.list) <(sort $barrnap) > reliable_gca.list
        """
}

process DerepRefGenomes {
    input:
        file genomes
        file reliable

    output:
        path 'derep_gca.list', emit: derep_gca

    script:
        """
        mkdir genomes
        for f in `cat $reliable`; do mv \$f*.fna genomes; done
        rm -f *.fna
        dRep dereplicate DREP -g genomes/*.fna -p ${params.cpu} 2>log
        ls DREP/dereplicated_genomes/*.fna | perl -nle 'm/(GC[AF]_\\d{9}\\.\\d+)/ ; print \$1' > derep_gca.list
        """
}

process SeqSelector {
    input:
        path FastaFiles

    output:
        path '*-uniq.fasta', emit: UniqFasta
        path 'gca.list', emit: gca_list

    script:
        """
        #!/usr/bin/env perl
        use Modern::Perl '2011';
        use autodie;
        use Smart::Comments '###';
        use Getopt::Euclid qw(:vars);
        use List::AllUtils qw( :all );
        use Bio::MUST::Core;
        use aliased 'Bio::MUST::Core::Ali';
        use Bio::MUST::Core::Utils qw(:filenames secure_outfile);

        my \$outfile = 'gca.list';
        open my \$out, '>', \$outfile;
        my \$regex = qr/(GC[AF]_\\d{9}\\.\\d+)/;
        my @gca = map { /\$regex/ ? \$1 : () } split ' ', "${FastaFiles}";
        say {\$out} \$_ for @gca; 

        FILE:
        for my \$infile (split ' ', "${FastaFiles}") {
            my \$ali = Ali->load(\$infile);
            if (Ali->instant_count(\$infile) == 1) {
                \$ali->store_fasta(change_suffix(insert_suffix(\$infile, '-uniq'), '.fasta'));
                next FILE;
            }
            elsif (Ali->instant_count(\$infile) > 1) {
                warn "[WARN] Dereplicated file \$infile have more than 1 sequence. Selecting the longest sequence.";
                my @longest_seqs = ( max_by { \$_->nomiss_seq_len } \$ali->all_seqs )[0];
                my \$ali2 = Ali->new( seqs => \\@longest_seqs );
                \$ali2->store_fasta(change_suffix(insert_suffix(\$infile, '-uniq'), '.fasta'));
            }
        }
        """
}

process ChangeIdsAli {
    input:
        path Infiles
        val taxdir

    output:
        path '*-long.fasta', emit: Outfiles

    script:
        """        
        perl -nle 'm/(GC[AF]_\\d{9}\\.\\d+)/ ; print \$1' *.fasta | sort -u > infiles.gca
        fetch-tax.pl infiles.gca --taxdir=${taxdir} --org-mapper --item-type=taxid
        change-ids-ali.pl *.fasta --org-mapper=infiles.org-idm --mode=abbr2long --out=-long
        """
}

