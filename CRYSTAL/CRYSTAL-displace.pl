#!/usr/bin/env perl

# Copyright (c) "2012, by Georgia Institute of Technology
#                Contributors: Alexandr Fonari
#                Affiliation: Dr. Bredas group
#                URL: https://github.com/alexandr-fonari/Main/tree/master/CRYSTAL
#                License: MIT License
# Version 0.99
####### ===== ## =====
# This is script is a fork of:
#
#   Program: xdatgen
#   Sung Sakong, PhD
#   sung.sakong _at_ uni-due.de
#   ver 0.5  4. June 2007
#   url: http://www.uni-due.de/~hp0058/vmdplugins/utilities/xdatgen.c
#
# it is heavily based on:
# http://theory.cm.utexas.edu/vasp/scripts/src/Vasp.pm
#
#### I hope you will like it! ####

use strict; 
use warnings;
use Data::Dumper;
use Math::Trig qw(acos_real);

# physical constants in eV, Ang and s
use constant PI => 4*atan2(1, 1);
use constant PlanckConstant => 4.135667516e-15; # [eV s]
use constant HBAR => PlanckConstant/(2*PI); # [eV s]
use constant CL => 2.99792458e18;    # [A/s]
use constant AM => 931.494043e6;     # [eV/c^2]
use constant CM2EV => 0.0001239842573148; # [eV]

if ( scalar( @ARGV ) < 1 )
{
    die( "\nUse: $0 <freq.out> <float factor>\n" );
}
open( my $outcar_fh, "<", $ARGV[0]) || die "Can't open $ARGV[0]: $!\n";

my %frag1_atoms;
open( my $frag1_fh, "<", frag1.list) || die "Can't open frag1.list: $!\n";
while(my $line = <$frag1_fh>)
{
    $line = trim($line);
    my ($label, $rest) = split /\s+/, $line, 2;
    my @coords = split '\s+', $rest;
    $frag1_atoms{$label} = [ @coords ];
}
close($frag1_fh);

my %frag2_atoms;
open( my $frag2_fh, "<", frag2.list) || die "Can't open frag2.list: $!\n";
while(my $line = <$frag2_fh>)
{
    $line = trim($line);
    my ($label, $rest) = split /\s+/, $line, 2;
    my @coords = split '\s+', $rest;
    $frag2_atoms{$label} = [ @coords ];
}
close($frag2_fh);

my $fact = 1.0;
if ( scalar( @ARGV ) == 2 )
{
    $fact = $ARGV[1];
}

# from http://cpansearch.perl.org/src/MKHRAPOV/Chemistry-MolecularMass-0.1/MolecularMass/MolecularMass.pm
my %m = ("H" => 1.00794, "D" => 2.014101, "T" => 3.016049, "He" => 4.002602, "Li" => 6.941, "Be" => 9.012182, "B" => 10.811, "C" => 12.0107, "N" => 14.00674, "O" => 15.9994, "F" => 18.9984032, "Ne" => 20.1797, "Na" => 22.989770, "Mg" => 24.3050, "Al" => 26.981538, "Si" => 28.0855, "P" => 30.973761, "S" => 32.066, "Cl" => 35.4527, "Ar" => 39.948, "K" => 39.0983, "Ca" => 40.078, "Sc" => 44.955910, "Ti" => 47.867, "V" => 50.9415, "Cr" => 51.9961, "Mn" => 54.938049, "Fe" => 55.845, "Co" => 58.933200, "Ni" => 58.6934, "Cu" => 63.546, "Zn" => 65.39, "Ga" => 69.723, "Ge" => 72.61, "As" => 74.92160, "Se" => 78.96, "Br" => 79.904, "Kr" => 83.80, "Rb" => 85.4678, "Sr" => 87.62, "Y" => 88.90585, "Zr" => 91.224, "Nb" => 92.90638, "Mo" => 95.94, "Tc" => 98, "Ru" => 101.07, "Rh" => 102.90550, "Pd" => 106.42, "Ag" => 107.8682, "Cd" => 112.411, "In" => 114.818, "Sn" => 118.710, "Sb" => 121.760, "Te" => 127.60, "I" => 126.90447, "Xe" => 131.29, "Cs" => 132.90545, "Ba" => 137.327, "La" => 138.9055, "Ce" => 140.116, "Pr" => 140.90765, "Nd" => 144.24, "Pm" => 145, "Sm" => 150.36, "Eu" => 151.964, "Gd" => 157.25, "Tb" => 158.92534, "Dy" => 162.50, "Ho" => 164.93032, "Er" => 167.26, "Tm" => 168.93421, "Yb" => 173.04, "Lu" => 174.967, "Hf" => 178.49, "Ta" => 180.9479, "W" => 183.84, "Re" => 186.207, "Os" => 190.23, "Ir" => 192.217, "Pt" => 195.078, "Au" => 196.96655, "Hg" => 200.59, "Tl" => 204.3833, "Pb" => 207.2, "Bi" => 208.98038, "Po" => 209, "At" => 210, "Rn" => 222, "Fr" => 223, "Ra" => 226, "Ac" => 227, "Th" => 232.038, "Pa" => 231.03588, "U" => 238.0289, "Np" => 237, "Pu" => 244, "Am" => 243, "Cm" => 247, "Bk" => 247, "Cf" => 251, "Es" => 252, "Fm" => 257, "Md" => 258, "No" => 259, "Lr" => 262, "Rf" => 261, "Db" => 262, "Sg" => 266, "Bh" => 264, "Hs" => 269, "Mt" => 268, "Uun" => 271, "Uuu" => 272,);

my (@a_masses, @a_labels, $total_atoms, $total_freqs, @a_indx, $basis, $coordinates, @e_values, @e_vectors);

while(my $line = <$outcar_fh>)
{
    if($line =~ m/ATOMS IN THE UNIT CELL:\s+(\d+)/)
    {
        $total_atoms = $1; $total_freqs = 3*$total_atoms;
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
            push(@a_indx, $temp[1]);
            push(@a_labels, ucfirst($temp[2]));
            push(@a_masses, ucfirst($m{$temp[2]}));

            $coordinates->[$i][0] = $temp[3]; # x
            $coordinates->[$i][1] = $temp[4]; # y
            $coordinates->[$i][2] = $temp[5]; # z
        }
    }
    if($line =~ /NORMAL MODES NORMALIZED TO CLASSICAL AMPLITUDES/)
    {
        # suppose to work with arbitrary number of frequencies
        for(my $k = 0; $k < $total_freqs;)
        {
            <$outcar_fh>;
            my $line = <$outcar_fh>;
            $line = trim($line);
            my ($trash, @freqs) = split('\s+', $line);
            push(@e_values, @freqs);
            # print Dumper(@freqs);

            <$outcar_fh>;
            for(my $j = 0; $j < $total_freqs; $j++)
            {
                my @t = (<$outcar_fh> =~ m/(-?\d+\.\d+)/g);
                # print Dumper(@t);
                for(my $i = 0; $i < scalar(@freqs); $i++)
                {
                    $e_vectors[$k+$i] .= ' '.$t[$i];
                }
            }
            $k += scalar(@freqs);
        }
    }
}
close($outcar_fh);

# my @x = ($a_cart_pos_x[2], $a_cart_pos_y[2], $a_cart_pos_z[2]);my @v = cart2frac(\@x, \@g);printf("%10.6f %10.6f %10.6f || %10.6f %10.6f %10.6f\n", $a_frac_pos_x[2], $a_frac_pos_y[2], $a_frac_pos_z[2], $v[0], $v[1], $v[2]);
# print Dumper(@a_labels, @a_masses, @a_cart_pos_x, @a_cart_pos_y, @a_cart_pos_z);
# print Dumper(@e_values, @e_vectors);

open( my $poscar_cart_fh, ">", "DISPCAR_cart" ) || die "Can't open DISPCAR_cart file: $!";
# open( my $poscar_recp_fh, ">", "DISPCAR_recp" ) || die "Can't open DISPCAR_resp file: $!";

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
        my $header = sprintf("POSCAR: disp= %d, w= %8.5f meV, qi0= %5e, amp= %5e, f= %6.4f\n", $d, $ev*CM2EV*1000, $qi0, $amp, $fact);
        print $poscar_cart_fh $header;

        for( my $j = 0; $j < $total_atoms; $j++)
        {
            my $sqrtm = sqrt($a_masses[$j]);
            my @pos = ($coordinates->[$j][0], $coordinates->[$j][1], $coordinates->[$j][2]);
            my @dv = ($disps[3*$j], $disps[3*$j+1], $disps[3*$j+2]);

            # in CRYSTAL normal modes are normlaized (divided by) to classical amplitudes: A=Sqrt(2E/k)
            # thus, displacement will be: dx = dx_orig*A*qi0
            my $qi0_cry = $fact*$qi0*sqrt(2*HBAR*CL/sqrt(AM));

            # dividing by square roots of the mass as follows from the characteristic length and amplitude formulas
            @dv = map { $_*$qi0_cry*$d/($sqrtm*sqrt($sqrtm)) } @dv;
            @pos = ($pos[0]+$dv[0], $pos[1]+$dv[1], $pos[2]+$dv[2]);
            print $poscar_cart_fh sprintf("%s %3d %15.12f %15.12f %15.12f\n", $a_labels[$j], $j+1, @pos);

        }
        print $poscar_cart_fh "\n";
    }
}

print "DISPCAR files have been created.\n\n";
close($poscar_cart_fh);
# close($poscar_recp_fh);

# SUBs
sub trim{ my $s=shift; $s =~ s/^\s+|\s+$//g; return $s;}
# http://stackoverflow.com/questions/10193750/perl-count-unique-elements-in-array
# counts consistently uniq and count_uniq, important for VASP
sub uniq       {my %c; my @u = grep !$c{$_}++, @_; keys(%c);}
sub count_uniq {my %c; my @u = grep !$c{$_}++, @_; values(%c);}

