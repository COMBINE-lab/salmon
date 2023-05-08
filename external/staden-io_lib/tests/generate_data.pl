#!/usr/bin/perl -w
use strict;

exit 0 if (-e "$ENV{srcdir}/.done");

# Loads a fasta file into a hash
sub load_fasta {
    my ($fn) = @_;
    my %seqs;

    open(my $fd, "<", $fn) || die "$fn: $!\n";
    my $name = undef;
    my $s = "";
    while (<$fd>) {
	chomp($_);

	if (/^>(\S+)/) {
	    if (defined($name)) {
		$seqs{$name} = $s;
	    }
	    $name = $1;
	    $s = "";
	} else {
	    $s .= $_;
	}
    }
    $seqs{$name} = $s;

    close($fd);

    return \%seqs;
}

#---- Load seq
print "Loading ce.fa\n";

my $seqs = load_fasta("$ENV{srcdir}/data/ce.fa");
my @names = keys(%$seqs);
my %len;
my $n;
my @base = qw/A C G T/;


if (! -w "$ENV{srcdir}/data") {
    chmod(0755, "$ENV{srcdir}/data") || die "$ENV{srcdir}/data: $!";
    chmod(0755, "$ENV{srcdir}")      || die "$ENV{srcdir}: $!";
}

#---- Generate sorted data
print "Generating ce#sorted.sam\n";
srand 15551;
open(my $out, ">", "$ENV{srcdir}/data/ce#sorted.sam") ||
    die "$ENV{srcdir}/data/ce#sorted.sam: $!";

# Create @SQ headers
foreach (sort @names) {
    my $len = length($seqs->{$_});
    print $out "\@SQ\tSN:$_\tLN:$len\n";
    $len{$_} = $len;
}

# Sequence lines
$n = 1;
foreach my $chr (sort @names) {
    my $len = $len{$chr};
    my $s = $seqs->{$chr};
    for (my $i=0; $i<$len-100; $i++) {
	if (rand() < 0.5) { #50x coverage
	    my $dna = substr($s,$i,100);
	    for (my $j=0; $j<5; $j++) {
		substr($dna, 100*rand(), 1) = $base[4*rand()];
	    }

           print $out "Seq",$n++,"\t0\t$chr\t",$i+1,"\t2\t100M\t*\t0\t0\t$dna\t*\n";
	}
    }
}
close($out) || die;


#---- Generate unsorted data
print "Generating ce#unsorted.sam\n";
srand 15551;
open($out, ">", "$ENV{srcdir}/data/ce#unsorted.sam") ||
    die "$ENV{srcdir}/data/ce#unsorted.sam: $!";

# Create @SQ headers
foreach (sort @names) {
    my $len = length($seqs->{$_});
    print $out "\@SQ\tSN:$_\tLN:$len\n";
    $len{$_} = $len;
}

# Sequence lines
$n = 1;
for (my $i = 0; $i < 500000; $i++) {
    my $chr = $names[$#names*rand()];
    my $pos = int(rand() * ($len{$chr}-100));
    my $dna = substr($seqs->{$chr},$pos,100);
    for (my $j=0; $j<5; $j++) {
	substr($dna, 100*rand(), 1) = $base[4*rand()];
    }

    print $out "Seq",$n++,"\t0\t$chr\t",$pos+1,"\t2\t100M\t*\t0\t0\t$dna\t*\n";
}
close($out) || die;


#---- Generate unsorted data

open(DONE, ">$ENV{srcdir}/.done")||die;
close(DONE);

