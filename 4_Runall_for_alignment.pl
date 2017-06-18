#!/usr/bin/perl
use strict;
use warnings;
use Parallel::ForkManager;
use List::Compare;

chdir "/mnt/lustre/users/wyunjin/leiming/Plant_phylogenomes/results/2015_1_26_orthology_assignment_subtree/";

main:{
    
	my @pep_files=glob "$ARGV[0]/*pep.fasta";
	&parallel_multiply_alignment(\@pep_files);

	my @mafft_files=glob "$ARGV[0]/*cds_mafft.fasta";
	&parallel_multiply_Gblock(\@mafft_files);

}

sub parallel_multiply_alignment{
    my  ($pep_files)=@_; 
    my   $pm=new Parallel::ForkManager(7);
    foreach my $pep_file (@$pep_files){
        $pm->start and next;
        &multiply_alignment($pep_file);
        $pm->finish;
    }
    $pm->wait_all_children;
}

sub parallel_multiply_Gblock{
    my ($cds_mafft_files)=@_; 
    my  $pm=new Parallel::ForkManager(10);
    foreach my $cds_mafft_file(@$cds_mafft_files){
	    print $cds_mafft_file, "\n";
        $pm->start and next;
        &Gblock($cds_mafft_file);
        $pm->finish;
    }
    $pm->wait_all_children; 
}

sub multiply_alignment{
    my ($pep_file)=@_; 
    my $pwdname=$1,my $num = $2, if($pep_file=~/^(.*?(\d+))\_pep.fasta/);
    system("einsi --thread 3  --quiet $pep_file  > $pwdname\_pep_mafft.fasta");
    system(" /home/Ddong/leiming/bin/pal2nal.pl   $pwdname\_pep_mafft.fasta   $pwdname\_cds.fasta  -output  fasta  >  $pwdname\_cds_mafft.fasta");
}

sub Gblock{
    my ($file)=@_;  
	print $file,"\n";
	my $pwdname=$1,if($file=~/^(.*?(\d+))\_cds_mafft.fasta/);
    chomp(my $num=(`grep ">" $file -c`)); 
    system("Gblocks $file -t=c  -b5=h  -e=.gb");
    my $sequence=`grep -v -P ">" $file.gb`;
    $sequence=~s/\s//g; 
	my $length=length($sequence);
    my $seq_len=$length/$num;
	open OUT ,">", "$pwdname\_cds_mafft_gblock.fasta";
    if ($seq_len >= 99){
       open IN  ,"<","$file.gb";
       while(my $line=<IN>){
	       if($line =~ /\S/){ 
              $line=~s/>|\s//g ;
              $/=">";
              my $seq=<IN>;
              $/="\n";
              $seq=~s/>//g;
              $seq=~s/\n/ /g;
              chop $seq;
              print OUT ">",$line,"\n",$seq,"\n";
		   }
       }
    }
    
	close OUT;
    my $gblock_file="$pwdname\_cds_mafft_gblock.fasta";
    &translating_dna($pwdname, $gblock_file);    
}

sub translating_dna{
    my ($pwdname,$file)=@_;
    open IN, "<", $file;
    print  $file,"\n" ;
    open OUT,">","$pwdname\_translation_cds.fasta";
    my $seqs = q{};
    my $header_seqs = q{};
    while (my $line=<IN>){
           my $sequence=<IN>;
           $sequence=~s/\s//g;
           my $protein=&dnapeptide($sequence);
           $seqs .= $protein; 
            $header_seqs.=$line.$protein."\n";
    }
    print OUT $header_seqs,"\n";
	close OUT; 
}

sub dnapeptide { 
    my($dna) = @_ ;
    use strict;
    use warnings;
    my $protein = '';
    for(my $i=0; $i < (length($dna) - 2) ; $i += 3) {
        $protein .= &codonaa(substr($dna,$i,3));
    }
     return $protein;
} 

sub codonaa {
    my($codon) = @_;
       if ( $codon =~ /TCA/i )    { return 'S' }    # Serine 
    elsif ( $codon =~ /TCC/i )    { return 'S' }    # Serine 
    elsif ( $codon =~ /TCG/i )    { return 'S' }    # Serine 
    elsif ( $codon =~ /TCT/i )    { return 'S' }    # Serine 
    elsif ( $codon =~ /TCN/i )    { return 'S' }    # Serine 
    elsif ( $codon =~ /TTC/i )    { return 'F' }    # Phenylalanine 
    elsif ( $codon =~ /TTT/i )    { return 'F' }    # Phenylalanine 
    elsif ( $codon =~ /TTA/i )    { return 'L' }    # Leucine 
    elsif ( $codon =~ /TTG/i )    { return 'L' }    # Leucine 
    elsif ( $codon =~ /TTN/i )    { return 'X' }    # Leucine 
    elsif ( $codon =~ /TAC/i )    { return 'Y' }    # Tyrosine 
    elsif ( $codon =~ /TAT/i )    { return 'Y' }    # Tyrosine 
    elsif ( $codon =~ /TAN/i )    { return 'X' }    # Tyrosine  
    elsif ( $codon =~ /TAA/i )    { return '_' }    # Stop 
    elsif ( $codon =~ /TAG/i )    { return '_' }    # Stop 
    elsif ( $codon =~ /TGC/i )    { return 'C' }    # Cysteine 
    elsif ( $codon =~ /TGT/i )    { return 'C' }    # Cysteine 
    elsif ( $codon =~ /TGN/i )    { return 'X' }    # Cysteine 
    elsif ( $codon =~ /TGA/i )    { return '_' }    # Stop 
    elsif ( $codon =~ /TGG/i )    { return 'W' }    # Tryptophan 
    elsif ( $codon =~ /CTA/i )    { return 'L' }    # Leucine  
    elsif ( $codon =~ /CTC/i )    { return 'L' }    # Leucine 
    elsif ( $codon =~ /CTG/i )    { return 'L' }    # Leucine 
    elsif ( $codon =~ /CTT/i )    { return 'L' }    # Leucine 
    elsif ( $codon =~ /CTN/i )    { return 'L' }    # Leucine 
    elsif ( $codon =~ /CCA/i )    { return 'P' }    # Proline 
    elsif ( $codon =~ /CCC/i )    { return 'P' }    # Proline 
    elsif ( $codon =~ /CCG/i )    { return 'P' }    # Proline 
    elsif ( $codon =~ /CCT/i )    { return 'P' }    # Proline 
    elsif ( $codon =~ /CCN/i )    { return 'P' }    # Proline 
    elsif ( $codon =~ /CAC/i )    { return 'H' }    # Histidine 
    elsif ( $codon =~ /CAT/i )    { return 'H' }    # Histidine 
    elsif ( $codon =~ /CAA/i )    { return 'Q' }    # Glutamine 
    elsif ( $codon =~ /CAG/i )    { return 'Q' }    # Glutamine 
    elsif ( $codon =~ /CAN/i )    { return 'X' }    # Glutamine 
    elsif ( $codon =~ /CGA/i )    { return 'R' }    # Arginine 
    elsif ( $codon =~ /CGC/i )    { return 'R' }    # Arginine 
    elsif ( $codon =~ /CGG/i )    { return 'R' }    # Arginine 
    elsif ( $codon =~ /CGT/i )    { return 'R' }    # Arginine 
    elsif ( $codon =~ /CGN/i )    { return 'R' }    # Arginine 
    elsif ( $codon =~ /ATA/i )    { return 'I' }    # Isoleucine 
    elsif ( $codon =~ /ATC/i )    { return 'I' }    # Isoleucine  
    elsif ( $codon =~ /ATT/i )    { return 'I' }    # Isoleucine 
    elsif ( $codon =~ /ATN/i )    { return 'X' }    # Isoleucine 
    elsif ( $codon =~ /ATG/i )    { return 'M' }    # Methionine 
    elsif ( $codon =~ /ACA/i )    { return 'T' }    # Threonine 
    elsif ( $codon =~ /ACC/i )    { return 'T' }    # Threonine 
    elsif ( $codon =~ /ACG/i )    { return 'T' }    # Threonine 
    elsif ( $codon =~ /ACT/i )    { return 'T' }    # Threonine 
    elsif ( $codon =~ /ACN/i )    { return 'T' }    # Threonine 
    elsif ( $codon =~ /AAC/i )    { return 'N' }    # Asparagine 
    elsif ( $codon =~ /AAT/i )    { return 'N' }    # Asparagine 
    elsif ( $codon =~ /AAA/i )    { return 'K' }    # Lysine 
    elsif ( $codon =~ /AAG/i )    { return 'K' }    # Lysine 
    elsif ( $codon =~ /AAN/i )    { return 'X' }
    elsif ( $codon =~ /AGC/i )    { return 'S' }    # Serine 
    elsif ( $codon =~ /AGT/i )    { return 'S' }    # Serine 
    elsif ( $codon =~ /AGA/i )    { return 'R' }    # Arginine 
    elsif ( $codon =~ /AGG/i )    { return 'R' }    # Arginine 
    elsif ( $codon =~ /AGN/i )    { return 'X' }    # Arginine 
    elsif ( $codon =~ /GTA/i )    { return 'V' }    # Valine 
    elsif ( $codon =~ /GTC/i )    { return 'V' }    # Valine 
    elsif ( $codon =~ /GTG/i )    { return 'V' }    # Valine 
    elsif ( $codon =~ /GTT/i )    { return 'V' }    # Valine 
    elsif ( $codon =~ /GTN/i )    { return 'V' }    # Valine 
    elsif ( $codon =~ /GCA/i )    { return 'A' }    # Alanine 
    elsif ( $codon =~ /GCC/i )    { return 'A' }    # Alanine 
    elsif ( $codon =~ /GCG/i )    { return 'A' }    # Alanine 
    elsif ( $codon =~ /GCT/i )    { return 'A' }    # Alanine
    elsif ( $codon =~ /GCN/i )    { return 'A' }    # Alanine
    elsif ( $codon =~ /GAC/i )    { return 'D' }    # Aspartic Acid 
    elsif ( $codon =~ /GAT/i )    { return 'D' }    # Aspartic Acid 
    elsif ( $codon =~ /GAA/i )    { return 'E' }    # Glutamic Acid 
    elsif ( $codon =~ /GAG/i )    { return 'E' }    # Glutamic Acid 
    elsif ( $codon =~ /GAN/i )    { return 'X' }
    elsif ( $codon =~ /GGA/i )    { return 'G' }    # Glycine 
    elsif ( $codon =~ /GGC/i )    { return 'G' }    # Glycine 
    elsif ( $codon =~ /GGG/i )    { return 'G' }    # Glycine 
    elsif ( $codon =~ /GGT/i )    { return 'G' }    # Glycine 
    elsif ( $codon =~ /GGN/i )    { return 'G' }    # Glycine 
    elsif ( $codon =~ /.N./i )    { return 'X' }
    elsif ( $codon =~ /N../i)     { return 'X' }
    elsif ( $codon =~ /.NN/i )    { return 'X' }
    elsif ( $codon =~ /NN./i )    { return 'X' }
    elsif ( $codon =~ /N.N/i)     { return 'X' }
    elsif ( $codon =~ /NNN/i)     { return 'X' }
    elsif ( $codon =~ /---/i)     { return '-' }
    else  {
                    print STDERR "Bad codon \"$codon\"!!\n";
            exit;
    }
}

#__END__
