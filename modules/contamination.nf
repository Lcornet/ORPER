process CheckM {
    input:
        file all_genomes

    output:
        path 'checkm_result.csv', emit: result

    script:
        println "[ORPER-INFO] Using CheckM as contamination detection tool"
        """
        rm -f *empty-abbr.fna
        mkdir genomes
        mv *.fna genomes/
        checkm lineage_wf -t $params.cpu -x fna genomes runc > checkm_result
        echo "#genome,completeness,contamination" > part1
        tr -s " " < checkm_result | grep "GC" | cut -f2,14,15 -d" " > part2
        sed -i -e 's/ /,/g' part2
        cat part1 part2 > checkm_result.csv
        """
}

process EukCC {
    input:
        file all_genomes
        path eukccdir

    output:
        path 'eukcc_result.csv', emit: result

    script:
        println "[ORPER-INFO] Using EukCC as contamination detection tool"
        """
        rm -f *empty-abbr.fna
        mkdir bins
        cd bins
        for f in ../*.fna; do ln -s \$f; done
        for FILE in *.fna; do mv \$FILE `basename \$FILE .fna`.fa; done
        cd ../
        # count the number of *.fa files
        fa_file_count=\$(find ./bins -name "*.fa" | wc -l)
        if [ \$fa_file_count -eq 1 ]; then
            eukcc single *.fna --out eukcc_result --threads $params.cpu --db $eukccdir
        else
            eukcc folder bins --out eukcc_result --threads $params.cpu --db $eukccdir
        fi
        cut -f1,2,3 eukcc_result/eukcc.csv | sed 's/bin/#genome/; s/.fa//; s/\t/,/g' > eukcc_result.csv
        """
}
