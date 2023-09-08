#!/usr/bin/perl

# PLEASE SEE LICENSE FILE FOR COPYRIGHT (all rights reserved) AND LICENSE (UoI/NCSA)

our (@name,@fa,@qu);

#TODO support unpaired fasta / fastq
my $state;
my $isFasta = undef;
my $lastName;
while (<>) {
    if (not defined $isFasta) {
        if (!/^([>\@])/) {
            die "Invalid fasta or fastq!";
        }
        if ($1 eq '>') {
            $isFasta = 1;
        } else {
            $isFasta = 0;
        }
    }
    if ($isFasta) {
        if (!/^>(\S+)\/(\d)$/) {
            die "Invalid fasta line: $_!";
        }
        if ($2 eq '1' || (defined $lastName && $lastName ne $1)) {
             printCrossbow();
        }
        $state = $2;
        $name[$state] = $lastName = $1;
        chomp($fa[$state] = <>);
    } else {
        if (!/^@(\S+)\/(\d)$/) {
            die "Invalid fastq line: $_!";
        }
        if ($2 eq '1') {
             printCrossbow();
        }
        $state = $2;
        $name[$state] = $lastName = $1;
        chomp($fa[$state] = <>);
        my $dummy = <>;
        chomp($qu[$state] = <>);
    }        
}
printCrossbow();

sub printCrossbow {
    return unless defined $name[1] && defined $fa[1] && defined $name[2] && defined $fa[2];
    print join("\t", $name[1], $fa[1], defined $qu[1] ? $qu[1] :'a' x length($fa[1]), $fa[2], defined $qu[2] ? $qu[2] : 'a' x length($fa[2])) . "\n";
    @name = ();
    @fa = ();
    @qu = ();
}
