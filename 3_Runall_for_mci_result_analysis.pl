#!/usr/bin/perl
use strict;
use warnings;
use List::Util qw/sum/;
use Parallel::ForkManager;
use List::Compare;

my $groups_file = $ARGV[0];
my @species_id = ("Amborella_trichopoda","Physcomitrella_patens","Picea_abies","Sample_RS101L","Sample_RS103S","Sample_RS107L","Sample_RS108S","Sample_RS10S","Sample_RS111S","Sample_RS112L","Sample_RS114S","Sample_RS115S","Sample_RS116S","Sample_RS11S","Sample_RS122L","Sample_RS123L","Sample_RS124L","Sample_RS127L","Sample_RS128S","Sample_RS137","Sample_RS14S","Sample_RS16","Sample_RS17","Sample_RS18","Sample_RS1","Sample_RS24L","Sample_RS27L","Sample_RS28","Sample_RS31","Sample_RS37","Sample_RS38S","Sample_RS39","Sample_RS43L","Sample_RS45S","Sample_RS46","Sample_RS47L","Sample_RS48L","Sample_RS52S","Sample_RS5S","Sample_RS70","Sample_RS77L","Sample_RS7S","Sample_RS84","Sample_RS85L","Sample_RS88S","Sample_RS89S","Sample_RS8L","Sample_RS90S","Sample_RS92","Sample_RS97","Sample_RS98S","Sample_SH_10","Sample_SH_11","Sample_SH_12","Sample_SH_13","Sample_SH_14","Sample_SH_15","Sample_SH_16","Sample_SH_17","Sample_SH_18","Sample_SH_19","Sample_SH_1","Sample_SH_20","Sample_SH_21","Sample_SH_2","Sample_SH_3","Sample_SH_4","Sample_SH_5","Sample_SH_6","Sample_SH_7","Sample_SH_8","Sample_SH_9","Selaginella_moellendorffii");
my @new_sequence = ("Sample_RS101L","Sample_RS103S","Sample_RS107L","Sample_RS108S","Sample_RS10S","Sample_RS111S","Sample_RS112L","Sample_RS114S","Sample_RS115S","Sample_RS116S","Sample_RS11S","Sample_RS122L","Sample_RS123L","Sample_RS124L","Sample_RS127L","Sample_RS128S","Sample_RS137","Sample_RS14S","Sample_RS16","Sample_RS17","Sample_RS18","Sample_RS1","Sample_RS24L","Sample_RS27L","Sample_RS28","Sample_RS31","Sample_RS37","Sample_RS38S","Sample_RS39","Sample_RS43L","Sample_RS45S","Sample_RS46","Sample_RS47L","Sample_RS48L","Sample_RS52S","Sample_RS5S","Sample_RS70","Sample_RS77L","Sample_RS7S","Sample_RS84","Sample_RS85L","Sample_RS88S","Sample_RS89S","Sample_RS8L","Sample_RS90S","Sample_RS92","Sample_RS97","Sample_RS98S","Sample_SH_10","Sample_SH_11","Sample_SH_12","Sample_SH_13","Sample_SH_14","Sample_SH_15","Sample_SH_16","Sample_SH_17","Sample_SH_18","Sample_SH_19","Sample_SH_1","Sample_SH_20","Sample_SH_21","Sample_SH_2","Sample_SH_3","Sample_SH_4","Sample_SH_5","Sample_SH_6","Sample_SH_7","Sample_SH_8","Sample_SH_9");

print scalar @species_id,"\n";
print scalar @new_sequence,"\n";

my $co_orthology_dir = $ARGV[1];
my $onetoone_orthology_dir = $ARGV[2];
my $seq_files_dir = $ARGV[3];

mkdir "$co_orthology_dir";
mkdir "$onetoone_orthology_dir";

open IN, "<", $groups_file;
my @lines = <IN>; 
close IN;

my $pm=new Parallel::ForkManager(18);
foreach my $line (@lines){
   $pm->start and next;
   &extract_homologous_genes($line,$co_orthology_dir,$onetoone_orthology_dir,$seq_files_dir);
   $pm->finish;
}
$pm->wait_all_children;

sub extract_homologous_genes{
	my ($line,$co_orthology_dir,$onetoone_orthology_dir,$seq_files_dir) = @_;
    my %hash;	
    my %names_hash;
	my %species_patterns;
	my $group_name = $1, my $genes = $2, if($line =~ /^(\S+)\s+(.*?)$/);

    while ($genes =~ /((\S+)__(\S+))/msg){
	    my $id = $1;
	    my $pattern = $3;	
	    my $species = $2;
		push @{$species_patterns{$species}}, $pattern;
	    $hash{$species}++;
    }
 
	foreach my $key (keys %hash){
        if ($hash{$key} > 10){
 	        delete $hash{$key};	
 		}	
	}

	my @orth_species = keys %hash;
    my @inter;
	&compare(\@orth_species, \@new_sequence, \@inter);

    my @values = values %hash;
	if(scalar @values >=1 ){
        my  $species_num = scalar @values;
        my $sum = sum @values; 
    	if(($sum/$species_num) == 1){
            if(scalar @inter >= 35){
		        print scalar @inter,"\n";	  
                &extract_homologous_genes_sequence($group_name,\@orth_species,\%species_patterns,$onetoone_orthology_dir,$seq_files_dir); 
			}
	   }
	   else{
        if(scalar @inter >= 35){ 
			print  scalar @inter,"\n";
             &extract_homologous_genes_sequence($group_name,\@orth_species,\%species_patterns,$co_orthology_dir,$seq_files_dir); 
	    } 
	   }
    }
}

sub compare{
	my ($select_species,$min_species,$inter)=@_;
	my $lc = List::Compare -> new($select_species,$min_species);
	@$inter= $lc -> get_intersection;
}

sub extract_homologous_genes_sequence{
    my ($group_name,$orth_species_array_ref,$species_patterns_hash_ref,$dir,$seq_files_dir)=@_; 
	#foreach my $key (keys %$species_patterns_hash_ref){
    foreach my $key (@$orth_species_array_ref){
        foreach my $element (@{$$species_patterns_hash_ref{$key}}){
            $element =~ s/\|/\\|/g;
            my $cmd = "awk -v o=$key -v c=$element '{if(\$1==c)print \"\>\"o\"__\"c\"\\n\"\$2}' $seq_files_dir/$key\_cds.fa.v1  >>    $dir/$group_name\_cds.fasta";
            print $cmd,"\n";
            system  $cmd;   
            $cmd="awk -v o=$key -v c=$element  '{if(\$1==c)print \"\>\"o\"__\"c\"\\n\"\$2}' $seq_files_dir/$key\_pep.fa.v1  >>  $dir/$group_name\_pep.fasta";
            print $cmd,"\n";
            system  $cmd;    
       }
    }
}

__END__
