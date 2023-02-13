#! /bin/perl
#

open(FO, "<$ARGV[0]");
open(FP, ">$ARGV[1]/tmp.txt");
my $h = <FO>;
while(<FO>)
{
        chomp;
        $_ =~ s/\r//g;
	$_ =~ s/\"//g;
        my $x = $_;
        my $dir = $x."_Red.idat";
        my $dir2 = $x."_Grn.idat";
        print FP "$dir\n";
        print FP "$dir2\n";
}
close(FP);
close(FO);

open(FO, "<$ARGV[1]/tmp.txt");
while(<FO>)
{
        chomp;
        system ("ln -s $_ $ARGV[1]/.");
}
close(FO);

system("rm $ARGV[1]/tmp.txt")
