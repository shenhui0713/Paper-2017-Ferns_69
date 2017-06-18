#!/usr/bin/perl
use strict;
use warnings;
use Cwd;

main:{
	my $current_path = cwd();
        my $quality_outdir  = "$ARGV[0]"; 
	my $assembly_outdir = $ARGV[1];

	mkdir "$assembly_outdir";
	mkdir "$quality_outdir";

        #Quality control of raw reads sequence 
        my @reads_files = glob "$current_path/raw_files/*R1.fastq";
        foreach my $reads_file (@reads_files){
            &quality_control ($reads_file, $quality_outdir);                            
        }

	#De novo transcriptome assembly
	my @trimmed_files = glob "$quality_outdir/*R1.fastq.trimmed.paired1";
        foreach my $trimmed_file (@trimmed_files){
            &de_novo_assembly ($trimmed_file, $assembly_outdir);
        }

}

#####subprograms#####
sub quality_control{
    my ( $reads_file , $quality_outdir) = @_;
    my $prefix_name = $1, my $file_name = $2, if ($reads_file =~ /^(.*?\/(.*?))_R1.fastq/);
    my  $sickle_cmd = "sickle  pe  -t sanger   -f   $reads_file  -r  $prefix_name\_R2.fastq -t illumina -o $quality_outdir/$file_name\_R1.trimmed.fastq  -p $quality_outdir/$file_name\_R2.trimm.fastq -s $quality_outdir/$file_name\_s.trimmed.fastq";
    system ($sickle_cmd);

    `perl /home/Ddong/leiming/src/SolexaQA_v.2.5/LengthSort.pl      $quality_outdir/$file_name\_R1.fastq.trimmed  $quality_outdir/$file_name\_R2.fastq.trimmed -l 25 -d  $quality_outdir `;
}

sub de_novo_assembly{
    my  ( $trimmed_file, $assembly_outdir)=@_;
    my  $prefix_name = $1, my $file_name = $2, if ($trimmed_file=~/^(.*?\/(.*?))\_R1.fastq.trimmed.paired1/);
    my  $trinity_cmd="ulimit -s unlimited; perl /home/Ddong/leiming/src/trinityrnaseq_r20140413p1/Trinity --bflyHeapSpaceMax 2G  --seqType	fq   --JM  60G --left  $trimmed_file   --right $prefix_name\_R1.fastq.trimmed.paired2   --CPU  7  -o $assembly_outdir/$file_name\_assembly  >	$assembly_outdir/$file_name\_assembly.nohup.out";
    system($trinity_cmd);
}

__END__
