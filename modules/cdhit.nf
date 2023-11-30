process CdHit {
    input:
        path Fasta
        val Threshold

    output:
        path "*-cdhit.fasta", emit: ClstFasta

    script:
        """
        for FILE in `echo "${Fasta}" | tr ' ' '\n'`; do
            filename=\$(basename "\${FILE}")
            filename="\${filename%.*}"
            cd-hit -i \${FILE} -o \${filename}-cdhit.fasta -c ${Threshold} -T ${params.cpu}
        done
        """
}
