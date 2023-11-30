// parameters
params.ribodb  = 'local'

// default directory
directory = file("$workflow.projectDir/ORPER-ribodb")

process SetupRiboDB {
    input:
        val ribodb

    output:
        val ribodir
    
    script:

        // set up ribodir variable
        ribodir = 'na'

        // if ribodb parameter is not specify..
        if (params.ribodb == 'local'){
            println "[ORPER-INFO] RiboDB not specified: searching in current directory"

            // check if 'ORPER-ribodb' exist
            // ... if not, it create it
            if ( !directory.exists() ) {
                println "[ORPER-INFO] RiboDB not found in project: creating ORPER-ribodb in current directory"

                // exit process if cannot create directory
                if ( !directory.mkdirs() ) {
                    exit 1, "[ORPER-WARN] Cannot create ORPER-ribodb directory"
                }

                ribodir = directory

                dir = (params.kingdom == 'Bacteria') ? 'prokaryotes/' : ((params.kingdom == 'Eukaryota') ? 'eukaryotes/' : 'na')

                """
                cd $directory
                git clone https://bitbucket.org/phylogeno/42-ribo-msas/
                mv 42-ribo-msas/MSAs/$dir/*.ali .
                ali2fasta.pl *.ali --degap 2> log
                find *.fasta | cut -f1 -d'.' > ribo.list
                rm -f *.ali
                echo $directory > ribodb_path.txt
                """
            }

            else {
                println "[ORPER-INFO] RiboDB available in current directory: using ORPER-ribodb"

                ribodir = directory

                """
                echo $directory > ribodb_path.txt
                """           
            }
        }

        else {
            println "[ORPER-INFO] RiboDB specified: using $params.ribodb"

            ribodir = ribodb

            """
            echo $ribodb > ribodb_path.txt
            """		
        }
}
