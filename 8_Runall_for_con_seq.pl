#!/usr/bin/perl

use warnings;
use strict;
use List::Compare;     
use Parallel::ForkManager;


my @species_id = ("Amborella_trichopoda","Physcomitrella_patens","Picea_abies","Sample_RS101L","Sample_RS103S","Sample_RS107L","Sample_RS108S","Sample_RS10S","Sample_RS111S","Sample_RS112L","Sample_RS114S","Sample_RS115S","Sample_RS116S","Sample_RS11S","Sample_RS122L","Sample_RS123L","Sample_RS124L","Sample_RS127L","Sample_RS128S","Sample_RS137","Sample_RS14S","Sample_RS16","Sample_RS17","Sample_RS18","Sample_RS1","Sample_RS24L","Sample_RS27L","Sample_RS28","Sample_RS31","Sample_RS37","Sample_RS38S","Sample_RS39","Sample_RS43L","Sample_RS45S","Sample_RS46","Sample_RS47L","Sample_RS48L","Sample_RS52S","Sample_RS5S","Sample_RS70","Sample_RS77L","Sample_RS7S","Sample_RS84","Sample_RS85L","Sample_RS88S","Sample_RS89S","Sample_RS8L","Sample_RS90S","Sample_RS92","Sample_RS97","Sample_RS98S","Sample_SH_10","Sample_SH_11","Sample_SH_12","Sample_SH_13","Sample_SH_14","Sample_SH_15","Sample_SH_16","Sample_SH_17","Sample_SH_18","Sample_SH_19","Sample_SH_1","Sample_SH_20","Sample_SH_21","Sample_SH_2","Sample_SH_3","Sample_SH_4","Sample_SH_5","Sample_SH_6","Sample_SH_7","Sample_SH_8","Sample_SH_9","Selaginella_moellendorffii");


my @minspecies = ("Sample_RS101L","Sample_RS103S","Sample_RS107L","Sample_RS108S","Sample_RS10S","Sample_RS111S","Sample_RS112L","Sample_RS114S","Sample_RS115S","Sample_RS116S","Sample_RS11S","Sample_RS122L","Sample_RS123L","Sample_RS124L","Sample_RS127L","Sample_RS128S","Sample_RS137","Sample_RS14S","Sample_RS16","Sample_RS17","Sample_RS18","Sample_RS1","Sample_RS24L","Sample_RS27L","Sample_RS28","Sample_RS31","Sample_RS37","Sample_RS38S","Sample_RS39","Sample_RS43L","Sample_RS45S","Sample_RS46","Sample_RS47L","Sample_RS48L","Sample_RS52S","Sample_RS5S","Sample_RS70","Sample_RS77L","Sample_RS7S","Sample_RS84","Sample_RS85L","Sample_RS88S","Sample_RS89S","Sample_RS8L","Sample_RS90S","Sample_RS92","Sample_RS97","Sample_RS98S","Sample_SH_10","Sample_SH_11","Sample_SH_12","Sample_SH_13","Sample_SH_14","Sample_SH_15","Sample_SH_16","Sample_SH_17","Sample_SH_18","Sample_SH_19","Sample_SH_1","Sample_SH_20","Sample_SH_21","Sample_SH_2","Sample_SH_3","Sample_SH_4","Sample_SH_5","Sample_SH_6","Sample_SH_7","Sample_SH_8","Sample_SH_9");


print scalar @species_id,"\t",scalar @minspecies,"\n";
####main####
main:{

    my  $outdir=$ARGV[0];
    &seq_concatation($outdir,\@species_id,\@minspecies);

}

####subprogram####
sub seq_concatation{
    my ($outdir,$species,$minspecies)=@_;
    my $pwdname=$1,if($outdir=~/(\S+)$/);
    my @files=glob "$outdir/*translation_cds.fasta";my $type="PRO";
    &concatation($pwdname,\@files,$type,$species,$minspecies);
	@files=glob "$outdir/*cds_mafft_gblock.fasta";$type="DNA";
	&concatation($pwdname,\@files,$type,$species,$minspecies);
}

sub concatation{ 
    my($pwdname,$files,$type,$species,$minspecies)=@_;
    if($type eq "DNA"){
       open OUT,">>","$pwdname\_con_DNA.fasta";
    }else{
       open OUT,">>","$pwdname\_con_pro.fasta";
    }
    my %con_seq;
    my $start=1;
    foreach my $file(@$files){
        my $filename=$1,if($file=~ /^(.*?)\ .fasta/);
        my $stop_num=`grep "_"  $file -c`;chomp $stop_num;
        open FILEONE,"<",$file;
        my $num=`grep ">"  $file -c`;chomp $num;
        if($num > 0){
        my @temp_species;my @inter;
        while (my $line=<FILEONE>){
            if($line=~/>(.*?)__/){  
                push @temp_species,$1;
            }
        }
        
		&compare_2(\@temp_species,$minspecies,\@inter);
        if (scalar @inter > 0){
            my %temp_hash;
            my $length;
            open IN,"<",$file;
            my $name=$1,if($file=~/\/(.*?)\.fasta/);
            while(my $line=<IN>){
                if($line=~/>(.*?)__/){
                    my $key=$1;
                    my $seq=<IN>;$seq=~s/\s//g;
                    $length=length $seq;
                    $temp_hash{$key}=$seq;
                }   
            }
            close IN;

            my @unique;my @spe=keys (%temp_hash);
            &compare($species,\@spe,\@unique);
            if($type eq "DNA"){
                foreach my $element(@unique){
                    $con_seq{$element}.=("N"x$length);
                }
            }else{
                foreach my $element(@unique){
                    $con_seq{$element}.=("X"x$length);  
                }
            }

            foreach my $spe(@spe){
                $con_seq{$spe}.=$temp_hash{$spe};
            }        


            }
     }
     }
  
     foreach (keys %con_seq){
              print OUT ">",$_,"\n",$con_seq{$_},"\n";
     }
     close OUT;
}

sub compare{
    my ($species,$spe,$unique)=@_;
    my $lc = List::Compare -> new($species,$spe);
    @$unique=$lc -> get_unique;
}

sub compare_2{
    my ($select_species,$min_species,$inter)=@_;
    my $lc = List::Compare -> new($select_species,$min_species);
    @$inter= $lc -> get_intersection;
}

__END__
