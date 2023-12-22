process RiboDBFortytwo {
    input:
        file FastaFiles
        path ribodir
        val taxdir

    output:
        path '*enrich.fasta', emit: enrich_files

    script:
        path_ref = params.kingdom == 'Bacteria' ? '42-ribo-msas/ref_orgs/prokaryotes/*.fasta'
            : '42-ribo-msas/ref_orgs/eukaryotes/*.faa'
        path_idm = params.kingdom == 'Bacteria' ? '42-ribo-msas/ref_orgs/prokaryotes/mapper-prokaryotes.idm'
            : '42-ribo-msas/ref_orgs/eukaryotes/mapper-eukaryotes.idm'
        suffix = params.kingdom == 'Bacteria' ? '.fasta' : '.faa'
        """
        #Reference organisms
        mkdir ref-bank && cd ref-bank/
        git clone https://bitbucket.org/phylogeno/42-ribo-msas/
        mv ${path_ref} ./
        for REFORG in *${suffix}; do makeblastdb -in \$REFORG -dbtype prot -out `basename \$REFORG ${suffix}` -parse_seqids; done
        mv ${path_idm} ./ref-mapper.idm
        cd ..

        #Query organisms
        grep -h '>' ${ribodir}/*.fasta | cut -f1 -d'@' | sort | uniq -c | sort -rh \
        | tail -n+2 | head -50 | cut -f2 -d'>' | sed 's/_/ /' > queries.idl

        #Candidate organisms
        mkdir candidate_orgs && mv *abbr.faa candidate_orgs
        cd candidate_orgs/
        for ORG in *-abbr.faa; do makeblastdb -in \$ORG -dbtype prot -out `basename \$ORG .faa` -parse_seqids; done
        ls *-abbr.faa | perl -nle 'm/(GC[AF]_\\d{9}\\.\\d+)/ ; print join "\t", \$1, \$_' > gca_file.idm

        fetch-tax.pl gca_file.idm --taxdir=${taxdir} --col=1 --sep='\t' --item-type=taxid --org-mapper
        paste gca_file.org-idm gca_file.idm | cut -f1,4 | sed 's/.faa//' > bank-mapper.idm
        #ls *-abbr.faa | sed 's/.faa//' | perl -nle 'm/(GC[AF]_\\d{9}\\.\\d+)/ ; print join "\t", \$1, \$_' > bank-mapper.idm
        cd ..

        #Part for alignements
        mkdir ribodb
        cp ${ribodir}/*.fasta ribodb/

        #YAML generator
        yaml-generator-42.pl --run_mode=phylogenomic --out_suffix=-ORPER --queries queries.idl --evalue=1e-05 --homologues_seg=yes \
        --max_target_seqs=1000 --templates_seg=no --bank_dir candidate_orgs --bank_suffix=.psq --bank_mapper candidate_orgs/bank-mapper.idm --ref_brh=on \
        --ref_bank_dir ref-bank --ref_bank_suffix=.psq --ref_bank_mapper ref-bank/ref-mapper.idm --ref_org_mul=0.15 --ref_score_mul=0.99 \
        --trim_homologues=off --ali_keep_lengthened_seqs=keep --aligner_mode=off --tax_reports=off --tax_dir ${taxdir} \
        --megan_like --tol_check=off

        #Run forty-two
        forty-two.pl ribodb/*.fasta --config=config-ORPER.yaml --verbosity=1 --threads=${params.cpu}

        #extract new sequences from enriched ribodb
        cd ribodb/
        fasta2ali.pl *-ORPER.fasta
        grep -c "#NEW" *ORPER.fasta > count-enrich.list
        /opt/ORPER/ORPER-companion.py count-enrich.list --mode=forty
        for f in `cat enrich.list`; do grep -A1 "#NEW#" \$f.ali > \$f-enrich.ali; done
        ali2fasta.pl --degap --noguessing *enrich.ali
        mv *enrich.fasta ../
        cd ..
        """
}

process rRNAFortyTwo {
    input:
        path FastaFiles
        path taxdir

    output:
        path '*-42-split.fasta', emit: fasta_files

    script:
        """
        # Download NCBI rRNA DB
        wget https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Fungi/fungi.18SrRNA.fna.gz
        wget https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Fungi/fungi.28SrRNA.fna.gz
        wget https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Fungi/fungi.ITS.fna.gz
        #wget https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Bacteria/bacteria.16SrRNA.fna.gz
        gunzip *.gz
        perl -nle 'if (m/^>/) { m/^>(\\S+)\\s(\\S+)\\s(\\S+)\\s/ ; print ">\$2_\$3\\@\$1" } else { print }' fungi.28SrRNA.fna > fungi.28SrRNA-renamed.fna
        perl -i.bak -ple 's/\\]//g ; s/\\[//g' fungi.28SrRNA-renamed.fna
        grep '>' fungi.28SrRNA-renamed.fna | cut -f1 -d'@' | cut -c2- | sort | uniq -c | sort -rn | head -n20 | perl -ple 's/\\s+\\d+\\s//; s/_/ /' > queries.txt
        #perl -nle 'if (m/^>/) { m/^>(\\S+)\\s(\\S+)\\s(\\S+)\\s/ ; print ">\$2_\$3\\@\$1" } else { print }' bacteria.16SrRNA.fna >  bacteria.16SrRNA-renamed.fna
        #perl -i.bak -ple 's/\\]//g ; s/\\[//g' bacteria.16SrRNA-renamed.fna
        #grep '>'  bacteria.16SrRNA-renamed.fna | cut -f1 -d'@' | cut -c2- | sort | uniq -c | sort -rn | head -n20 | perl -ple 's/\\s+\\d+\\s//; s/_/ /' > queries.txt

        # Download mitochondrial rRNA sequence as `.para` files
        #for f in L01493.1 EU546103.1 KR704917.1 JQ303212.1 KR045776.1; do \$HOME/edirect/efetch -id \$f -db nuccore -format fasta; done > mitochondrial_16S.fasta
        for f in AM495067.1 AM495066.1 AM495076.1 AM495068.1 AM495062.1 ; do \$HOME/edirect/efetch -id \$f -db nuccore -format fasta; done > mitochondrial_23S.fasta
        ln -s mitochondrial_23S.fasta fungi.28SrRNA-renamed.para
        #ln -s mitochondrial_23S.fasta bacteria.16SrRNA-renamed.para

        # Bank
        for REFORG in `echo "${FastaFiles}" | tr ' ' '\n'`; do makeblastdb -in \$REFORG -dbtype nucl -out `basename \$REFORG .fna` -parse_seqids; done
        echo "${FastaFiles}" | tr ' ' '\n' | perl -nle 'm/(GC[AF]_\\d{9}\\.\\d+)/ ; print join "\t", \$1, \$_' > gca_file.idm
        fetch-tax.pl gca_file.idm --taxdir=${taxdir} --col=1 --sep='\t' --item-type=taxid --org-mapper
        paste gca_file.org-idm gca_file.idm | cut -f1,4 | sed 's/.fna//' > bank.idm

        # YAML generator
        yaml-generator-42.pl --run_mode=phylogenomic --SSUrRNA --out_suffix=-42 --queries queries.txt \
        --evalue=1e-05 --max_target_seqs=1000 \
        --bank_dir \$PWD --bank_suffix=.fna --bank_mapper bank.idm --code=1 \
        --ref_brh=off --trim_homologues=off --merge_orthologues=off \
        --aligner_mode=blast --ali_skip_self=off --ali_cover_mul=1.1 --ali_keep_old_new_tags=off --ali_keep_lengthened_seqs=keep \
        --tax_reports=on --tax_min_score=0 --tax_score_mul=0 --tax_min_ident=0 --tax_min_len=0 --tol_check=off

        # Run 42
        forty-two.pl fungi.28SrRNA-renamed.fna --config=config-42.yaml --verbosity=1
        #forty-two.pl bacteria.16SrRNA-renamed.fna --config=config-42.yaml --verbosity=1

        # Split results on individual FASTA
        fasta2ali.pl fungi.28SrRNA-renamed-42.fna --degap
        grep NEW fungi.28SrRNA-renamed-42.ali | cut -c2- > fungi.28SrRNA-renamed-42.idl
        prune-ali.pl fungi.28SrRNA-renamed-42.idl --out=-new
        grep '>' ./fungi.28SrRNA-renamed-42-new.ali | perl -nle 'm/(GC[AF]_\\d{9}\\.\\d+)/ ; print \$1' | sort -u > org_gca.lis
        for f in `cat org_gca.lis`; do grep -A1 \$f fungi.28SrRNA-renamed-42-new.ali > \$f-42-split.ali; done
        ali2fasta.pl *-42-split.ali
        perl -i.bak -ple 's/#NEW#//' *-42-split.fasta
        """
}

process ITSFortyTwo {
    input:
	path FastaFiles
        path taxdir

    output:
        path '*-ITS-split.fasta', emit: fasta_files

    script:
        """
        # Download fungal ITS file from NCBI rRNA DB
        wget https://ftp.ncbi.nlm.nih.gov/refseq/TargetedLoci/Fungi/fungi.ITS.fna.gz
        gunzip *.gz
        perl -nle 'if (m/^>/) { m/^>(\\S+)\\s(\\S+)\\s(\\S+)\\s/ ; print ">\$2_\$3\\@\$1" } else { print }' fungi.ITS.fna > fungi.ITS-renamed.fna
        perl -i.bak -ple 's/\\]//g ; s/\\[//g; s/\\://g' fungi.ITS-renamed.fna
        grep '>' fungi.ITS-renamed.fna | cut -f1 -d'@' | cut -c2- | sort | uniq -c | sort -rn | head -n20 | perl -ple 's/\\s+\\d+\\s//; s/_/ /' > queries.txt

        # Bank
        for REFORG in `echo "${FastaFiles}" | tr ' ' '\n'`; do makeblastdb -in \$REFORG -dbtype nucl -out `basename \$REFORG .fna` -parse_seqids; done
        echo "${FastaFiles}" | tr ' ' '\n' | perl -nle 'm/(GC[AF]_\\d{9}\\.\\d+)/ ; print join "\t", \$1, \$_' > gca_file.idm
        fetch-tax.pl gca_file.idm --taxdir=${taxdir} --col=1 --sep='\t' --item-type=taxid --org-mapper
        paste gca_file.org-idm gca_file.idm | cut -f1,4 | sed 's/.fna//' > bank.idm

        # YAML generator
        yaml-generator-42.pl --run_mode=phylogenomic --SSUrRNA --out_suffix=-42 --queries queries.txt \
        --evalue=1e-05 --max_target_seqs=1000 --bank_dir \$PWD --bank_suffix=.fna --bank_mapper bank.idm --code=1 \
        --ref_brh=off --trim_homologues=on --trim_max_shift=20000 --trim_extra_margin=15 \
        --merge_orthologues=off --aligner_mode=blast --ali_skip_self=off --ali_cover_mul=1.1 --ali_keep_old_new_tags=off \
        --ali_keep_lengthened_seqs=keep --tax_reports=on --tax_min_score=0 --tax_score_mul=0 --tax_min_ident=0 --tax_min_len=0 \
        --tol_check=off

        # Run 42
        forty-two.pl fungi.ITS-renamed.fna --config=config-42.yaml --verbosity=1

        # Split results on individual FASTA
        fasta2ali.pl fungi.ITS-renamed-42.fna --degap
        grep NEW fungi.ITS-renamed-42.ali | cut -c2- > fungi.ITS-renamed-42.idl
        prune-ali.pl fungi.ITS-renamed-42.idl --out=-new
        grep '>' ./fungi.ITS-renamed-42-new.ali | perl -nle 'm/(GC[AF]_\\d{9}\\.\\d+)/ ; print \$1' | sort -u > org_gca.lis
        for f in `cat org_gca.lis`; do grep -A1 \$f fungi.ITS-renamed-42-new.ali > \$f-ITS-split.ali; done
        ali2fasta.pl *-ITS-split.ali
        perl -i.bak -ple 's/#NEW#//' *-ITS-split.fasta
        """
}

