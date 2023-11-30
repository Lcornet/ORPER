process Barrnap {
    input:
        file all_genomes

    output:
        //path 'all_barrnap.fna', emit: result
        path '*barrnap.fna', emit: result

    script:
        kingdom = params.kingdom == 'Bacteria' ? 'bac' : 'euk'
        """
        rm -f *empty-abbr.fna
        find *.fna | sed 's/-abbr.fna//' > fna.list
        for f in `cat fna.list`; do barrnap \$f-abbr.fna --outseq \$f-barrnap.fna --threads 1 --kingdom $kingdom; done
        # count the number of '*-barrnap.fna' file
        file_count=\$(find ./ -name "*-barrnap.fna" | wc -l)
        if [ \$file_count -eq 1 ]; then
            if [ ! -s *-barrnap.fna ]; then
                echo -e "[ORPER-WARNING] The only fasta file have no predicted rRNA."
                exit 1
            fi
        fi
        #cat *barrnap.fna > all_barrnap.fna
        """
}

process FilterBarrnap {
    input:
        file barrnap_files
        val rrna

    output:
        path '*-barrnap-*.fna', emit: filter_files
        path 'gca.list', emit: gca_list

    script:
        """
        find *barrnap.fna | perl -nle 'm/(GC[AF]_\\d{9}\\.\\d+)/ ; print \$1' > files_basename.list
        for f in `cat files_basename.list`; do grep -A1 --no-group-separator ${rrna} \$f*-barrnap.fna > \$f-barrnap-${rrna}.fna || true; done
        find ./ -size 0 -print -delete
        perl -i.bak -ple 's/${rrna}\\_rRNA::// ; s/\\(.*\\)//' *-barrnap-${rrna}.fna
        find *-barrnap-${rrna}.fna | perl -nle 'm/(GC[AF]_\\d{9}\\.\\d+)/ ; print \$1' > gca.list
        """
}
