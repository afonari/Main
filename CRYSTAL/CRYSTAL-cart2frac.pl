#!/usr/bin/env perl
#;-*- Perl -*-

# Copyright (c) "2012, by Georgia Institute of Technology
#                Contributors: Alexandr Fonari
#                Affiliation: Dr. Bredas group
#                URL: https://github.com/alexandr-fonari/Main/tree/master/CRYSTAL
#                License: MIT License
# Version 0.99
#
# it is heavily based on (some SUBs are copy/pasted):
# http://theory.cm.utexas.edu/vasp/scripts/src/Vasp.pm
#
#### I hope you will like it! ####

use strict; 
use warnings;
use Data::Dumper;
use Math::Trig qw(acos_real);

if ( scalar( @ARGV ) < 1 )
{
    die( "\nUse: $0 <freq.out>\n" );
}
open( my $outcar_fh, "<", $ARGV[0]) || die "Can't open $ARGV[0]: $!\n";

my (%f_labels, $f_atoms);
open( my $f_fh, "<", "atoms.list") || die "Can't open atoms.list: $!\n";

my $iter = 0;
while(my $line = <$f_fh>)
{
    $line = trim($line);
    my ($label, @rest) = split /\s+/, $line;
    $f_labels{$label} = $iter;

    $f_atoms->[$iter][0] = $rest[1];
    $f_atoms->[$iter][1] = $rest[2];
    $f_atoms->[$iter][2] = $rest[3];

    $iter++;
}
close($f_fh);

my ($atoms, $total_atoms, $basis, $coordinates);
while(my $line = <$outcar_fh>)
{
    if($line =~ m/ATOMS IN THE UNIT CELL:\s+(\d+)/)
    {
        $total_atoms = $1;
        print "Number of atoms = $total_atoms;\n";
    }

    if($line =~ /DIRECT LATTICE VECTORS CARTESIAN COMPONENTS/)
    {
        <$outcar_fh>; # next line

        for (my $i=0; $i<3; $i++)
        {
            my $line = <$outcar_fh>; $line = trim($line);
            my @line = split(/\s+/,$line);
            # This is how Vasp reads in the basis
            for (my $j=0; $j<3; $j++)
            {
                $basis->[$j][$i] = $line[$j];
            }
        }

        <$outcar_fh>;<$outcar_fh>;<$outcar_fh>;
        <$outcar_fh>;<$outcar_fh>;<$outcar_fh>;

        for(my $i=0; $i < $total_atoms; $i++)
        {
            my $line = <$outcar_fh>; $line = trim($line);

            #  1     6 C     3.261901092826E+00 -9.881441284132E-01  9.509079807780E-01
            my @temp = split('\s+', $line);
            $atoms->[$i]{"indx"} = $temp[0];
            $atoms->[$i]{"number"} = $temp[1];
            $atoms->[$i]{"label"} = $temp[2];
            $atoms->[$i]{"mass"} = $m{$temp[2]};

            $coordinates->[$i][0] = $temp[3]; # x
            $coordinates->[$i][1] = $temp[4]; # y
            $coordinates->[$i][2] = $temp[5]; # z
        }
    }
}
close($outcar_fh);

# my @x = ($a_cart_pos_x[2], $a_cart_pos_y[2], $a_cart_pos_z[2]);my @v = cart2frac(\@x, \@g);printf("%10.6f %10.6f %10.6f || %10.6f %10.6f %10.6f\n", $a_frac_pos_x[2], $a_frac_pos_y[2], $a_frac_pos_z[2], $v[0], $v[1], $v[2]);
# print Dumper(@a_labels, @a_masses, @a_cart_pos_x, @a_cart_pos_y, @a_cart_pos_z);
# print Dumper(@e_values, @e_vectors);

$f1_atoms = dirkar($f1_atoms,$basis,scalar(keys(%f1_labels)));
$f2_atoms = dirkar($f2_atoms,$basis,scalar(keys(%f2_labels)));

open( my $poscar_cart_fh, ">", "DISPCAR_cart" ) || die "Can't open DISPCAR_cart file: $!";

my @displacements = (-2, -1, 1, 2); # hard-coded so far for 5-point stencil finite difference 1st derivative.
for( my $i = 0; $i < scalar(@e_values); $i++)
{
    print "processing ".($i+1)." out of ".scalar(@e_values)." eigenvalues\n";
    my $ev = $e_values[$i];
    if($ev < 0.0){next;} # skip imaginary frequency

    my $qi0 = sqrt((HBAR*CL)**2/(AM*$ev*CM2EV)); # characteristic length of a normal mode.
                                                 # Will divide by square root of atomic mass later in the code.
    my $amp = sqrt(2*HBAR*CL/sqrt(AM));          # classical amplitude

    my @disps = split('\s+', trim($e_vectors[$i]));
    foreach my $d (@displacements)
    {
        my ($frag1_out, $frag2_out);
        my $header = sprintf("POSCAR: disp= %d, w= %8.5f meV, qi0= %5e, amp= %5e, f= %6.4f\n", $d, $ev*CM2EV*1000, $qi0, $amp, $fact);
        print $poscar_cart_fh $header;

        for( my $j = 0; $j < $total_atoms; $j++)
        {
            my $sqrtm = sqrt($atoms->[$j]{"mass"});
            # in CRYSTAL normal modes are normlaized (divided by) to classical amplitudes: A=Sqrt(2E/k)
            # thus, displacement will be: dx = dx_orig*A*qi0
            # also, dividing by square roots of the mass as follows from the characteristic length and amplitude formulas
            my $qi0_cry = $fact*$qi0*$amp*$d/($sqrtm*sqrt($sqrtm));
            my @dv = ($disps[3*$j], $disps[3*$j+1], $disps[3*$j+2]);

            my $atom_name = $atoms->[$j]{"label"}.$atoms->[$j]{"indx"};
            if(exists($f1_labels{$atom_name}))
            {
                my $i = $f1_labels{$atom_name};
                my @pos = ($f1_atoms->[$i][0]+$dv[0], $f1_atoms->[$i][1]+$dv[1], $f1_atoms->[$i][2]+$dv[2]);
                $frag1_out .= sprintf("%s %9.5f %9.5f %9.5f\n", $atom_name, @pos);
            }
            elsif(exists($f2_labels{$atom_name}))
            {
                my $i = $f2_labels{$atom_name};
                my @pos = ($f2_atoms->[$i][0]+$dv[0], $f2_atoms->[$i][1]+$dv[1], $f2_atoms->[$i][2]+$dv[2]);
                $frag2_out .= sprintf("%s %9.5f %9.5f %9.5f\n", $atom_name, @pos);
            }
        }
        print $poscar_cart_fh $frag1_out;
        print $poscar_cart_fh "NEXT_FRAG\n";
        print $poscar_cart_fh $frag2_out;
    }
}

print "DISPCAR files have been created.\n\n";
close($poscar_cart_fh);
# close($poscar_recp_fh);

# SUBs
sub trim{ my $s=shift; $s =~ s/^\s+|\s+$//g; return $s;}
sub inverse {
    my $basis = shift;

    my $inverse;

    my $omega = 0;
    my ($i,$ii,$iii,$j,$jj,$jjj);

    for ($i=0; $i<3; $i++) {
        $ii = $i+1;
        if ($ii>2) { $ii-=3; }
        $iii = $ii+1;
        if ($iii>2) { $iii-=3; }
        for ($j=0;$j<3;$j++) {
            $jj = $j+1;
            if ($jj>2) { $jj-=3; }
            $jjj = $jj + 1;
            if ($jjj>2) { $jjj-=3; }
            $inverse->[$j][$i] = $basis->[$jj][$ii]*$basis->[$jjj][$iii] - $basis->[$jjj][$ii]*$basis->[$jj][$iii];
            # print "$i $ii $iii $j $jj $jjj: ".$inverse->[$j][$i]."\n";
        }
    }

    $omega = $inverse->[0][0]*$basis->[0][0] + $inverse->[1][0]*$basis->[1][0] + $inverse->[2][0]*$basis->[2][0];

    for ($i=0; $i<3; $i++) {
        for ($j=0; $j<3; $j++) {
            $inverse->[$i][$j] /= $omega;
        }
    }

    return($inverse);
}

sub kardir {
    my $vector = shift;
    my $basis = shift;
    my $lattice = shift;
    my $total_atoms = shift;
    my $recip_basis;
    my ($v1,$v2,$v3,$i,$j);
  
    $recip_basis = inverse($basis);
  
    for ($i=0; $i<$total_atoms; $i++) {
        $v1 = $vector->[$i][0]*$recip_basis->[0][0] + $vector->[$i][1]*$recip_basis->[1][0] + $vector->[$i][2]*$recip_basis->[2][0];
        $v2 = $vector->[$i][0]*$recip_basis->[0][1] + $vector->[$i][1]*$recip_basis->[1][1] + $vector->[$i][2]*$recip_basis->[2][1];
        $v3 = $vector->[$i][0]*$recip_basis->[0][2] + $vector->[$i][1]*$recip_basis->[1][2] + $vector->[$i][2]*$recip_basis->[2][2];

        # move atoms to primative cell
        $vector->[$i][0] = $v1+60-int($v1+60);
        $vector->[$i][1] = $v2+60-int($v2+60);
        $vector->[$i][2] = $v3+60-int($v3+60);
    }
    for ($i=0;$i<3;$i++) {
        for ($j=0;$j<3;$j++) {
            # print $basis->[$i][$j]." ";
        }
        # print " .... ";
        for ($j=0;$j<3;$j++) {
            #  print $recip_basis->[$i][$j]." ";
        }
        # print "\n";
    }
    return ($vector);
}

