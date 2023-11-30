process RaxmlReferenceTree {
    input:
        file ConcatFastaFile

    output:
        path 'reference.tre', emit: Tree
        path 'GC.list', emit: IdList

    script:
        """
        #transform matrix into phylip file  map-ids
        ali2phylip.pl data-ass.ali --map-ids
        #compute tree
        # -x -f a -N 100 will do a 100x rapid bootstrap analysis
        raxmlHPC-PTHREADS-AVX -T ${params.cpu} -s data-ass.phy -n data-ass-RAXML-PROTGAMMALGF-100xRAPIDBP -m PROTGAMMALGF -N 100 -f a -x 1975021703574 -p 1975021703574
        cp data-ass.idm RAxML_bipartitions.idm
        sed -i -e 's/-abbr//g' RAxML_bipartitions.idm
        cut -f1 -d"@" RAxML_bipartitions.idm > f1
        cut -f2 RAxML_bipartitions.idm > f2
        paste f1 f2 > RAxML_bipartitions.idm
        format-tree.pl RAxML_bipartitions.data-ass-RAXML-PROTGAMMALGF-100xRAPIDBP --map-ids
        mv RAxML_bipartitions.tre reference.tre
        tree2list.pl reference.tre
        grep -v "#" reference.idl > GC.list
        """
}

process RaxmlConstrainTree {
    input:
        file RefTree
        file RefTreeId
        file BarrnaprRNA
        file rRnaFile
        file OutGcaList
        val taxdir

    output:
        path 'constrained-tree.tre', emit: Tree
        path 'constrained-tree.nex', emit: NexusTree

        //file 'Constained-tree.nex' into constTreeNexus_ch
        //file 'Constained-tree.tre' into constTreeTre_ch
        //file 'SSU-combined-ali-a2p.fasta' into alignementBmge_ch
        //file 'SSU-combined-ali.fasta' into alignement_ch
        //file 'GC.list' into gcfList_ch

    script:
        """
        # Process input rRNA fasta file, corrected unix fasta file
        rrnafile=\$(basename "${rRnaFile}")
        rrnafile="\${rrnafile%.*}"
        ali2fasta.pl \${rrnafile}.fasta
        # Shorten sequence ID to first blank
        inst-abbr-ids.pl \${rrnafile}.fasta --id-regex=:DEF
        sed -i -e 's/|//g' \${rrnafile}-abbr.fasta

        # Ensure that rejected sequence in scafos are not present in predicted rRNA file
        cat `echo "${BarrnaprRNA}"` > all_rRNA_barrnap-uniq.fasta
        fasta2ali.pl all_rRNA_barrnap-uniq.fasta
        perl -nle 'm/(GC[AF]_\\d{9}\\.\\d+)/ ; print \$1' GC.list | sort -u > RefTree.gca
        for f in `cat RefTree.gca`; do grep -A1 \$f all_rRNA_barrnap-uniq.ali >> filtered_rRNA_barrnap-uniq.ali; done
        #perl -nle 'm/(GC[AF]_\\d{9}\\.\\d+)/ ; print \$1' filtered_rRNA_barrnap-uniq.ali | sort -u > filtered_rRNA_barrnap-uniq.gca
        #fetch-tax.pl filtered_rRNA_barrnap-uniq.gca --taxdir=${taxdir} --org-mapper --item-type=taxid
        #change-ids-ali.pl filtered_rRNA_barrnap-uniq.ali --org-mapper=filtered_rRNA_barrnap-uniq.org-idm --mode=abbr2long --out=-long
        ali2fasta.pl filtered_rRNA_barrnap-uniq.ali
        perl -i.bak -ple 's/@.*//' filtered_rRNA_barrnap-uniq.fasta
    
        # Combine predicted and input rRNA files
        cat filtered_rRNA_barrnap-uniq.fasta \${rrnafile}-abbr.fasta > combined-rRNA.fasta

        # Align comibined file and alignment filtering using BMGE
        mafft --adjustdirection --anysymbol --auto --reorder combined-rRNA.fasta 2>log > combined-rRNA_align.fasta
        ali2phylip.pl combined-rRNA_align.fasta --bmge-mask=medium --max=0.6 --p80
        raxmlHPC-PTHREADS-AVX -T $params.cpu -r "${RefTree}" -s combined-rRNA_align.p80 \
        -n combined-rRNA_align-RAXML-GTRGAMMA-100xRAPIDBP -m GTRGAMMA -N 100 -f a -x 1975021703574 -p 197502170357 

        # Delete long branch
        cp RAxML_bipartitions.combined-rRNA_align-RAXML-GTRGAMMA-100xRAPIDBP raxml.tre
        run_treeshrink.py -t raxml.tre
        perl -F"\\s" -anle 'next if m/Signature/; print \$F[2] if \$F[3] >= 0.1' raxml_treeshrink/output_summary.txt > long-branch-ids.list
        grep -vf `echo "${OutGcaList}"` long-branch-ids.list || true > constrained-tree.idl
        cp RAxML_bipartitions.combined-rRNA_align-RAXML-GTRGAMMA-100xRAPIDBP constrained-tree.tre
        prune-tree.pl --negate-list constrained-tree.idl
        format-tree.pl constrained-tree.tre --figtree
        """
}
