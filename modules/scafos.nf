process Scafos {
    input:
        file AlignFastaFiles

    output:
        path 'data-ass.ali', emit: ConcatFastaFile
    
    script:
        """
        #Convert aligned FASTA into ALI format
        mkdir aligned
        mkdir a2p
        perl -i.bak -ple 's/\\|/@/' *.fasta
        mv *.fasta aligned/
        cd aligned/
        fasta2ali.pl *.fasta
        #Alignment filtering with BMGE mask
        # TODO: handle ali files with no sequence
        ali2phylip.pl *.ali --bmge-mask=medium --min=0.5 --ali
        mv *a2p* ../a2p/
        cd ../

        #Produce concat with scafos
        scafos in=a2p out=otu
        scafos in=a2p out=data otu=otu/otu-freq.otu o=ov
        # TODO: manage this error
        scafos in=data out=data-ass otu=otu/otu-freq.otu gamma=yes o=gclv g=30 format=fpm 2>log || true
        cp data-ass/data-ass.ali ./
        """
}
