// parameters
params.eukccdb = 'local'

// default directory
directory = file("$workflow.projectDir/ORPER-eukccdb")

process SetUpEukccDB {
    input:
        val eukccdb

    output:
        val eukccdir

    script:

        // set up eukccdir variable
        eukccdir = 'na'

        // if eukccdb parameter is not specify..
        if (params.eukccdb == 'local'){
            println "[ORPER-INFO] EukccDB not specified: searching in current directory"

            // check if 'ORPER-eukccdb' exist
            // ... if not, it create it
            if( !directory.exists() ) {
                println "[ORPER-INFO] EukccDB not found in project: creating ORPER-eukccdb in current directory"

                // exit process if cannot create directory
                if( !directory.mkdirs() )    {
                    exit 1, "Cannot create working directory"
                }

                eukccdir = directory

                """
                wget http://ftp.ebi.ac.uk/pub/databases/metagenomics/eukcc/eukcc2_db_ver_1.2.tar.gz
                tar -xzf eukcc2_db_ver_1.2.tar.gz
                #mkdir $directory/eukcc
                mv eukcc2_db_ver_1.2/* $directory/
                #mv eukcc2_db_ver_1.2 $directory/eukcc/
                """
            }

            else {
                println "[ORPER-INFO] EukccDB available in current directory: using ORPER-eukccdb"

                eukccdir = directory

                """
                echo $directory > eukccdb_path.txt
                """
            }
        }

        else {
            println "[ORPER-INFO] EukccDB specified: using $params.eukccdb"

            eukccdir = eukccdb

            """
            echo $eukccdb > eukccdb_path.txt
            """
        }
}
