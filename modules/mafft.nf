process Mafft {
    input:
        file FastaFiles

    output:
        path '*-align.fasta', emit: AlignFastaFiles

    script:
        """
        perl -i.bak -ple 's/#NEW#//' *.fasta
        for f in `ls *.fasta | sed 's/.fasta//'`; do mafft --auto --reorder \$f.fasta > \$f-align.fasta; done
        """
}
