#!/usr/bin/perl
use strict;
use warnings;
use Parallel::ForkManager;
use Cwd;

main: {

    my $current_path = cwd();
    my @blastp_files = glob "$mblastp_dir/*blast.out";
    my $abc_file = $ARGV[0];
    foreach my $blastp_file (@blastp_files) {
        &generating_abc_file($blastp_file,$abc_file);
    }

    my  $mcxload_cmd = "mcxload -abc $abc_file -write-tab $abc_file.mci.dict -o $abc_file.mci --stream-mirror  --stream-n
eg-log10 -stream-tf 'ceil(200)'";
    system ($mcxload_cmd);
    my  @inflation_values = (1.4,1.6,1.8,2.0,,2.2,2.4,2.6);     
    &parallel_mcl_cluster( \@inflation_values, "$abc_file.mci");
}

#####subprograms#####

sub generating_abc_file {
    my ($blastp_file,$abc_file) = @_;
    if ( $blastp_file =~ /\/(\S+)_vs_(\S+)_blast.out/ ) {
        my $name1 = $1;
        my $name2 = $2;
        open IN,  "<", "$blastp_file";
        open OUT, ">>", "$abc_file";
        while ( my $line = <IN> ) {
            if ( $line =~ /^Query Id/ ) { next }
            if ( $line =~ /^(\S+.*?)\t(\S+.*?)\t(\S+)\s+/ ) {
                my $id1     = $1;
                my $id2     = $2;
                my $e_value = $3;
                $id1 =~ s/\s.*$//;
                $id2 =~ s/\s.*$//;
                if ( $e_value - 1e-10 < 0 ) {
                    print OUT $name1, "__", $id1, "\t", $name2, "__", $id2,
                      "\t", $e_value, "\n";
                }
            }
        }
        close IN;
        close OUT;
    }
}

sub parallel_mcl_cluster {
    my ( $inflation_vaules_array_ref, $abc_file ) = @_;
    my $pm = new Parallel::ForkManager(2);
    foreach my $inflation_value (@$inflation_vaules_array_ref) {
        $pm->start and next;
        &mcl_cluster( $inflation_value, $abc_file );
        $pm->finish;
    }
    $pm->wait_all_children;
}

sub mcl_cluster {
    my ( $inflation_value, $abc_file ) = @_;
    my $mcl_cmd =
"mcl  $abc_file  -I  $inflation_value  -o  $abc_file\_$inflation_value  -tf 'gq(20)' ";
    system($mcl_cmd);
    my $mcxdump_cmd =
"mcxdump -icl   $abc_file\_$inflation_value   -o   $abc_file\_$inflation_value.dump   -tabr $abc_file.dict";
    system($mcxdump_cmd);
}

__END__
