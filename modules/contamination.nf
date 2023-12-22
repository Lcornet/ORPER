process CheckM {
    input:
        path FastaFiles

    output:
        path 'checkm_result.csv', emit: result
        path 'checkm_result-*', emit: test
        path 'checkm-resume-*', emit: test2

    script:
        """
        #!/usr/bin/env perl
        use Modern::Perl '2011';
        use autodie;
        use Smart::Comments '###';
        use Template;
        use Data::UUID;
        use Parallel::Batch;
        use POSIX qw(ceil);
        use Bio::MUST::Core;
        use aliased 'Bio::MUST::Core::Ali';

        my @files = (split ' ', "${FastaFiles}");
        my @subsets = make_subset(@files);

        # determining the maximum process to run in parallel
        my \$cpu = 15;
        \$cpu = ${params.cpu} if scalar @files <= 150;
        my \$maxprocs = ceil( ${params.cpu}/15 );

        my \$batch = Parallel::Batch->new( {
            maxprocs => \$maxprocs,
            jobs     => \\@subsets,
            code     => sub {                       # closure (providing \$self)
                            my \$files = shift;
                            my \$ug = Data::UUID->new;
                            my \$uuid = \$ug->create_hex();
                            process(\$files, \$uuid);
                        },
        } );

        # launch jobs
        \$batch->run();

        # build the final results file
        system(q{echo "#genome,completeness,contamination" > title});
        system(q{cat title checkm-resume-* > checkm_result.csv});

        sub process {
            my \$files = shift;
            my \$uuid = shift;

            my \$infiles = join(' ', @\$files);

            # create bash command-line template
my \$template = <<'EOT';
mkdir genomes-[% uuid %]
mkdir results-[% uuid %]
mv [% infiles %] genomes-[% uuid %]/
rm -f genomes-[% uuid %]/*empty-abbr.fna
checkm lineage_wf -t [% cpu %] -x fna genomes-[% uuid %] results-[% uuid %] > checkm_result-[% uuid %]
tr -s " " < checkm_result-[% uuid %] | grep "GC" | cut -f2,14,15 -d" " > checkm-resume-[% uuid %]
sed -i -e 's/ /,/g' checkm-resume-[% uuid %]
EOT

            # build command
            my %vars = (
                infiles => \$infiles,
                cpu     => \$cpu,
                uuid    => \$uuid,
            );
            my \$command;
            my \$tt = Template->new;
            \$tt->process(\\\$template, \\%vars, \\\$command);
            system(\$command);

        }

        # split the list of files into subset of 150 files
        sub make_subset {
            my @files = @_;

            my @subsets;
            while (@files) {
                 push @subsets, [ splice @files, 0, 150 ];
            }

            return @subsets;
        }
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
