#!/usr/bin/env perl

# Copyright (c) "2012, by Georgia Institute of Technology
#                Contributors: Alexandr Fonari
#                Affiliation: Dr. Bredas group
#                URL: https://github.com/alexandr-fonari/Main/tree/master/CRYSTAL
#                License: MIT License

use strict; 
use warnings;
use Data::Dumper;

# physical constants in eV, Ang and s
use constant PI => 4 * atan2(1, 1);
use constant PlanckConstant => 4.135667516e-15; # [eV s]
use constant HBAR => PlanckConstant/(2*PI); # [eV s]
use constant CL => 2.99792458e18;    # [A/s]
use constant AM => 931.494043e6;
use constant CM2EV => 0.0001239842573148; # [eV]

if ( scalar( @ARGV ) < 1 )
{
    die( "\nUse: $0 <freq.out>\n" );
}

# from http://cpansearch.perl.org/src/MKHRAPOV/Chemistry-MolecularMass-0.1/MolecularMass/MolecularMass.pm
my %m = ("H" => 1.00794, "D" => 2.014101, "T" => 3.016049, "He" => 4.002602, "Li" => 6.941, "Be" => 9.012182, "B" => 10.811, "C" => 12.0107, "N" => 14.00674, "O" => 15.9994, "F" => 18.9984032, "Ne" => 20.1797, "Na" => 22.989770, "Mg" => 24.3050, "Al" => 26.981538, "Si" => 28.0855, "P" => 30.973761, "S" => 32.066, "Cl" => 35.4527, "Ar" => 39.948, "K" => 39.0983, "Ca" => 40.078, "Sc" => 44.955910, "Ti" => 47.867, "V" => 50.9415, "Cr" => 51.9961, "Mn" => 54.938049, "Fe" => 55.845, "Co" => 58.933200, "Ni" => 58.6934, "Cu" => 63.546, "Zn" => 65.39, "Ga" => 69.723, "Ge" => 72.61, "As" => 74.92160, "Se" => 78.96, "Br" => 79.904, "Kr" => 83.80, "Rb" => 85.4678, "Sr" => 87.62, "Y" => 88.90585, "Zr" => 91.224, "Nb" => 92.90638, "Mo" => 95.94, "Tc" => 98, "Ru" => 101.07, "Rh" => 102.90550, "Pd" => 106.42, "Ag" => 107.8682, "Cd" => 112.411, "In" => 114.818, "Sn" => 118.710, "Sb" => 121.760, "Te" => 127.60, "I" => 126.90447, "Xe" => 131.29, "Cs" => 132.90545, "Ba" => 137.327, "La" => 138.9055, "Ce" => 140.116, "Pr" => 140.90765, "Nd" => 144.24, "Pm" => 145, "Sm" => 150.36, "Eu" => 151.964, "Gd" => 157.25, "Tb" => 158.92534, "Dy" => 162.50, "Ho" => 164.93032, "Er" => 167.26, "Tm" => 168.93421, "Yb" => 173.04, "Lu" => 174.967, "Hf" => 178.49, "Ta" => 180.9479, "W" => 183.84, "Re" => 186.207, "Os" => 190.23, "Ir" => 192.217, "Pt" => 195.078, "Au" => 196.96655, "Hg" => 200.59, "Tl" => 204.3833, "Pb" => 207.2, "Bi" => 208.98038, "Po" => 209, "At" => 210, "Rn" => 222, "Fr" => 223, "Ra" => 226, "Ac" => 227, "Th" => 232.038, "Pa" => 231.03588, "U" => 238.0289, "Np" => 237, "Pu" => 244, "Am" => 243, "Cm" => 247, "Bk" => 247, "Cf" => 251, "Es" => 252, "Fm" => 257, "Md" => 258, "No" => 259, "Lr" => 262, "Rf" => 261, "Db" => 262, "Sg" => 266, "Bh" => 264, "Hs" => 269, "Mt" => 268, "Uun" => 271, "Uuu" => 272,);

open( my $outcar_fh, "<", $ARGV[0]) || die "$!\n";

# atomic masses and labels
my (@a_masses, @a_labels, @a_count, $natoms, $nfreqs);

# basis vectors
my @f;

# fractional atomic positions
my (@a_frac_pos_x, @a_frac_pos_y, @a_frac_pos_z);

# obtain Cartesian coordinates: Transpose(basis)*v:
my (@a_cart_pos_x, @a_cart_pos_y, @a_cart_pos_z);

# hessian eigenvalues and eigenvectors
my (@e_values, @e_vectors);

while(my $line = <$outcar_fh>)
{
    if($line =~ m/ATOMS IN THE UNIT CELL:\s+(\d+)/)
    {
        $natoms = $1;
        print "Number of atoms = $natoms;\n";

        $nfreqs = 3*$natoms;
    }

    if($line =~ /DIRECT LATTICE VECTORS CARTESIAN COMPONENTS/)
    {
        $_ = <$outcar_fh>; # next line

        ($f[1], $f[2], $f[3]) = (<$outcar_fh> =~ m/([\w\d\.\-\+]+)/g);
        ($f[4], $f[5], $f[6]) = (<$outcar_fh> =~ m/([\w\d\.\-\+]+)/g);
        ($f[7], $f[8], $f[9]) = (<$outcar_fh> =~ m/([\w\d\.\-\+]+)/g);

        $_ = <$outcar_fh>;$_ = <$outcar_fh>;$_ = <$outcar_fh>;
        $_ = <$outcar_fh>;$_ = <$outcar_fh>;$_ = <$outcar_fh>;

        for(my $i = 1; $i <= $natoms; $i++)
        {
            my $line = <$outcar_fh>;
            $line = trim($line);

            #  1     6 C     3.261901092826E+00 -9.881441284132E-01  9.509079807780E-01
            my @t = split('\s+', $line);
            push(@a_masses, ucfirst($m{$t[2]}));
            push(@a_labels, ucfirst($t[2]));
            push(@a_cart_pos_x, $t[3]);
            push(@a_cart_pos_y, $t[4]);
            push(@a_cart_pos_z, $t[5]);
        }
    }

    if($line =~ /NORMAL MODES NORMALIZED TO CLASSICAL AMPLITUDES/)
    {
        for(my $k = 0; $k < $nfreqs;)
        {
            <$outcar_fh>;
            my $line = <$outcar_fh>;
            $line = trim($line);
            my ($trash, @freqs) = split('\s+', $line);
            push(@e_values, @freqs);
            # print Dumper(@freqs);

            <$outcar_fh>;            
            for(my $j = 1; $j <=  $nfreqs; $j++)
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

# print Dumper(@a_labels, @a_masses, @a_cart_pos_x, @a_cart_pos_y, @a_cart_pos_z);
# print Dumper(@e_values, @e_vectors);
close($outcar_fh);

open( my $poscar_fh, ">", "DISPCAR" ) || die "Can't open DISPCAR file: $!";

my @displacements = (-2, -1, 1, 2); # hard-coded so far for 5-point stencil finite difference 1st deriv.
for( my $i = 0; $i < scalar(@e_values); $i++)
{
    print "processing ".($i+1)." out ".scalar(@e_values)." eigenvalues\n";
    my $ev = $e_values[$i];
    if($ev < 0.0){next;} # skip imaginary frequency

    my $qi0 = sqrt((HBAR*CL)**2/(AM*$ev*CM2EV)); # mode quanta
    # printf("%15.12f %15.12f\n", sqrt($ev)*VaspToEv, $qi0);

    my @disps = split('\s+', trim($e_vectors[$i]));

    foreach (@displacements)
    {
        print $poscar_fh sprintf("POSCAR: disp=%d, w=%8.5f meV, qi0=%5e\n", $_, $ev*CM2EV, $qi0);
        print $poscar_fh "1.00000\n";
        print $poscar_fh sprintf("%15.12f %15.12f %15.12f\n", $f[1], $f[2], $f[3]);
        print $poscar_fh sprintf("%15.12f %15.12f %15.12f\n", $f[4], $f[5], $f[6]);
        print $poscar_fh sprintf("%15.12f %15.12f %15.12f\n", $f[7], $f[8], $f[9]);
        print $poscar_fh join(" ", uniq(@a_labels))."\n";
        print $poscar_fh join(" ", count_uniq(@a_labels))."\n";
        print $poscar_fh "Cartesian\n";

        for( my $j = 0; $j < scalar(@a_cart_pos_x); $j++)
        {
            my $sqrtm = sqrt($a_masses[$j]);
            my($dx, $dy, $dz) = ($disps[3*$j]*$qi0*$_/$sqrtm, $disps[3*$j+1]*$qi0*$_/$sqrtm, $disps[3*$j+2]*$qi0*$_/$sqrtm);
            print $poscar_fh sprintf("%s %15.12f %15.12f %15.12f\n", $a_labels[$j], $a_cart_pos_x[$j]+$dx, $a_cart_pos_y[$j]+$dy, $a_cart_pos_z[$j]+$dz);
        }
        print $poscar_fh "\n";
    }
    
}

print "DISPCAR created\n";
close($poscar_fh);

sub trim{ my $s=shift; $s =~ s/^\s+|\s+$//g; return $s;}

# http://stackoverflow.com/questions/10193750/perl-count-unique-elements-in-array
# consistent counting
sub uniq {my %c; my @u = grep !$c{$_}++, @_; keys(%c);}
sub count_uniq {my %c; my @u = grep !$c{$_}++, @_; values(%c);}

