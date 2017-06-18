#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw/sum/;
use Parallel::ForkManager;

my @pro_files = glob "$ARGV[0]/*translation_cds.fasta";

my  $pm=new Parallel::ForkManager(4);
foreach my $pro_file (@pro_files){
	$pm->start and next; 
	&phylogenetic_analysis($pro_file);
	$pm->finish;
}
$pm->wait_all_children;

sub phylogenetic_analysis{
	my ($pro_file) = @_;
    my $dir =$1, my $prefix = $2, if($pro_file =~ /^(.*)\/(.*?)\_translation_cds.fasta/);
   	my $num = `grep ">" $pro_file -c`;
    chomp $num;
    
	open IN, "<" ,$pro_file;
    open OUT, ">", "$pro_file.v1";
	print "$pro_file.v1\n";
	while (my $line = <IN>){
	    if($line =~ />/){
		   my $seq = <IN>;
		   unless ($seq =~/^-*$/){
			   print OUT  $line,$seq;
		   }
	    } 
	}
    close IN;
	close OUT;
	
	my $raxml_cmd = "~/leiming/src/standard-RAxML-master/raxmlHPC-PTHREADS-SSE3  -T  5  -p  12345  -s $pro_file.v1  -n  $prefix  -m  PROTCATWAG  -w /mnt/lustre/users/wyunjin/leiming/Plant_phylogenomes/results/2015_1_25_orthology_assignment/$dir";
    print $raxml_cmd,"\n";
    system ($raxml_cmd);
} 

__END__
