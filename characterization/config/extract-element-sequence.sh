#!/usr/bin/perl -w

$name = "linker$ARGV[0].lmpdat";

$n = `grep "atom types" $name | awk '{print \$1}'`; chomp($n);
$ln = `grep -n "Masses" $name`; chomp($ln); ($ln,$junk) = split(":",$ln);

#print "n = $n";
#print "ln = $ln";

$nh = $ln+1+$n;
@seq = `head -n $nh $name | tail -n $n | awk '{print \$4}'`; chomp(@seq);
#print "seq = @seq"

foreach $ln (@seq){($el,$junk) = split("_",$ln); print "$el ";}
print "\n";

$filename = "in.ESequence$ARGV[0]";
open(my $fh, '>', $filename) or die "Could not open file '$filename' $!";
print $fh "variable\tESequence string \"";
foreach $ln (@seq){($el,$junk) = split("_",$ln); print $fh "$el ";}
print $fh "\"\n";
close $fh;

