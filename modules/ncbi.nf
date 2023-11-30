// default directory
directory = file("$workflow.projectDir/taxdump")

//TODO: add tax and level input (allowing ref taxon and out taxon)
process GetURL {
    input:
	path taxdump
        file ftp
        val taxon
        val level

    output:
	path 'ftp_*.sh',   emit: ftp_links
        path 'gc*.list', emit: gca_list

    script:
        assembly = ftp =~ '.*refseq.txt' ? 'gcf' : 'gca'

        """
	grep -v "#" $ftp | cut -f1 > all_${assembly}.list
        fetch-tax.pl all_${assembly}.list --taxdir=$taxdump --item-type=taxid --levels=$level
        grep $taxon all_${assembly}.tax | cut -f1 > ${assembly}.list
        grep -v "#" $ftp | cut -f20 | sort -u | perl -F"/" -anle 'print "wget" . " " . \$_ . "/" . "\$F[9]" . "_genomic.fna.gz"' >> ftp_${assembly}.sh
        """
}

process GenomeDownloader {
    input:
	file ftp_links
        file gca_list

    output:
        path '*-abbr.fna', emit: abbr_files

    script:
        db = gca_list.name == 'gcf.list' ? 'refseq' : 'genbank'
        
	"""
        if [ -s $gca_list ]; then
            grep -f $gca_list $ftp_links > download.sh
            bash download.sh
            gunzip *.gz
            ls *.fna | perl -nle 'm/(GC[AF]_\\d{9}\\.\\d+)/ ; print join "\t", \$_, \$1' > file.idm
            inst-abbr-ids.pl --id-prefix-mapper=file.idm --id-regex=:DEF *.fna
        else
            touch ${db}_empty-abbr.fna
        fi
	"""
}

process ProteomeDownloader {
    input:
        file ftp_links
        file gca_list
        val org_category

    output:
        path '*-abbr.faa', emit: abbr_files

    script:
        """
        cat ftp_*.sh > combined_ftp.sh
        if [ "$org_category" = 'out' ]; then
            for f in `cat $gca_list | head -n10`; do grep \$f combined_ftp.sh ; done > reliable_combined_ftp.sh
        else
            for f in `cat $gca_list`; do grep \$f combined_ftp.sh ; done > reliable_combined_ftp.sh
        fi
        perl -i.bak -ple 's/_genomic.fna.gz/_protein.faa.gz/g' reliable_combined_ftp.sh
        bash reliable_combined_ftp.sh
        find ./ -size 0 -print -delete -name '*.gz'
        gunzip *.gz
        ls *.faa | perl -nle 'm/(GC[AF]_\\d{9}\\.\\d+)/ ; print join "\t", \$_, \$1' > file.idm
        inst-abbr-ids.pl --id-prefix-mapper=file.idm --id-regex=:DEF *.faa
        """
}

process SetupTaxonomy {
    input:
        val taxdump
    
    output:
        val taxdir

    script:

        // set up taxdir variable
        taxdir = 'na'

        // if ribodb parameter is not specify..
        if (params.taxdump == 'local'){
            println "[ORPER-INFO] Taxdump not specified: searching in current directory"

            // check if 'taxdump' exist
            // ... if not, it create it
            if( !directory.exists() ) {
                println "[ORPER-INFO] Taxdump not found in project: creating taxdump in current directory"

                // exit process if cannot create directory
                if( !directory.mkdirs() )    {
                    exit 1, "[ORPER-WARN] Cannot create taxdump directory"
                }

                taxdir = directory

                """
                setup-taxdir.pl --taxdir=$directory
                echo $directory > taxdump_path.txt
                """
            }

            else {
                println "[ORPER-INFO] Taxdump available in current directory: using taxdump"

                taxdir = directory 

                """
                echo $directory > taxdump_path.txt
                """           

            }
        }

        else {
            println "[ORPER-INFO] Taxdump specified: using $params.taxdump"

            taxdir = taxdump

            """
            echo $taxdump > taxdump_path.txt
            """
        }
}
