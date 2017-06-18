#!/usr/bin/perl
use strict;
use warnings;

my  @files = glob "$ARGV[0]/RAxML*";
mkdir "$ARGV[1]";

foreach my $file(@files){
    my $suffix = $1, if($file =~ /\.(\d+)$/);
    print $suffix,"\n";
    my $tree =`grep "\\w" $file`;
    $tree =~ s/__/@/g;
    $tree =~ s/\)\d+\:/\):/g;
    open OUT, ">", "$ARGV[1]/RAxML_bestTree\.$suffix";
    print OUT $tree;
    close OUT;		
}

`python ~/Universal_softwore_src/agalma-0.5.0/agalma/treeprune.py   $ARGV[1]   --id  $ARGV[1]    --outdir  $ARGV[1]  >  $ARGV[1].nohup.out`;

__END__
