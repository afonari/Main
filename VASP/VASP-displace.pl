#!/usr/bin/env perl
#;-*- Perl -*-

# Copyright (c) "2012, Alexandr Fonari
#                URL: https://github.com/alexandr-fonari/Main/tree/master/CRYSTAL
#                License: MIT License
# Version 0.99
#
# ===== ## =====
# This script is a fork of (copyrighted accordingly):
#
#   Program: xdatgen
#   Sung Sakong, PhD
#   sung.sakong _at_ uni-due.de
#   ver 0.5  4. June 2007
#   url: http://www.uni-due.de/~hp0058/vmdplugins/utilities/xdatgen.c
#
# Crystallographic conversions are heavily based on and copyrighted accordingly (some SUBs are copy/pasted):
# http://theory.cm.utexas.edu/vasp/scripts/src/Vasp.pm
# ===== ## =====
#
# I hope you will like it,
# Sasha.

use strict;
use warnings;
use XML::Simple;
use Data::Dumper;

# physical constants
use constant PI => 4 * atan2(1, 1);
use constant PlanckConstant => 4.135667516e-15; # [eV s]
use constant HBAR => PlanckConstant/(2*PI); # [eV s]
use constant CL => 2.99792458e18;    # [A/s]
use constant AM => 931.494043e6;
use constant Angstrom => 1.0e-10;    # [m]
use constant EV => 1.60217733e-19;   # [J]
use constant AMU => 1.6605402e-27;   # [kg]
use constant VaspToEv => sqrt(EV/AMU)/Angstrom/(2*PI)*PlanckConstant; # [eV] 6.46541380e-2
use constant VaspToCm =>  VaspToEv/(1.2398419e-4); # [cm^-1] 521.47083

my $xml = new XML::Simple();
my $data = $xml->XMLin("vasprun.xml");

my (@a_labels, @a_indx, @a_masses, @a_count, %label_mass);

my @t = @{$data->{"atominfo"}->{"array"}->{"atomtypes"}->{"set"}->{"rc"}};
foreach (@t)
{
    my ($label, $mass) = ( trim($_->{"c"}[1]), trim($_->{"c"}[2]) );
    $label_mass{$label} = $mass;
}
# print Dumper( \% label_mass);

@t = @{$data->{"atominfo"}->{"array"}->{"atoms"}->{"set"}->{"rc"}};
foreach (@t)
{
    push(@a_labels, trim($_->{"c"}[0]));
    push(@a_indx,   trim($_->{"c"}[1]));
    push(@a_masses, $label_mass{trim($_->{"c"}[0])});
}
# print Dumper( @a_labels, @a_indx, @a_masses );

# basis vectors
my $basis;
@t = @{$data->{"structure"}->{"initialpos"}->{"crystal"}->{"varray"}->{"basis"}->{"v"}};
for (my $i=0; $i<3; $i++)
{
    my @line = split(/\s+/, trim($t[$i]));
    # This is how Vasp reads in the basis
    for (my $j=0; $j<3; $j++)
    {
        $basis->[$j][$i] = $line[$j];
    }
}
# print Dumper( \%basis );

# fractional atomic positions and # of atoms
my ($coord_frac, $coord_cart, $natoms);

@t = @{$data->{"structure"}->{"initialpos"}->{"varray"}->{"v"}};
$natoms = scalar(@t);
for (my $i=0; $i<$natoms; $i++)
{
    my @t1 = split(/\s+/, trim($t[$i]));
    $coord_frac->[$i][0] = $t1[0]; # x
    $coord_frac->[$i][1] = $t1[1]; # y
    $coord_frac->[$i][2] = $t1[2]; # z
}
# print Dumper( \%coord_frac );

# obtain Cartesian coordinates: Transpose(basis)*v:
$coord_cart = dirkar($coord_frac, $basis, $natoms);
# print Dumper( \%$coord_cart );

# hessian eigenvalues
my @e_values = split(/\s+/, trim($data->{"calculation"}->{"dynmat"}->{"v"}->{"content"}));

# hessian eigenvectors
my @e_vectors = @{$data->{"calculation"}->{"dynmat"}->{"varray"}->{"eigenvectors"}->{"v"}};
# print Dumper(@e_vectors);

open( my $poscar_fh, ">", "DISPCAR" ) || die "Can't open DISPCAR file: $!";

for( my $i = 0; $i < scalar(@e_values); $i++)
{
    print "processing ".($i+1)." out ".scalar(@e_values)." eigenvalues\n";
    my $ev = $e_values[$i]*(-1.0);
    if($ev < 0.0){next;} # skip imaginary frequency

    my $qi0 = sqrt((HBAR*CL)**2/(AM*sqrt($ev)*VaspToEv)); # mode quanta
    # printf("%15.12f %15.12f\n", sqrt($ev)*VaspToEv, $qi0);

    my @disps = split('\s+', trim($e_vectors[$i]));

    my @displacements = (-1, -0.5, -0.1, -0.01, 0.5, 1); # hard-coded so far for the five-point stencil 1st deriv.
    foreach (@displacements)
    {
        print $poscar_fh sprintf("POSCAR: disp=%f, w=%8.5f meV\n", $_, sqrt($ev)*VaspToEv*1000);
        print $poscar_fh "1.00000\n";
        print $poscar_fh sprintf("%15.12f %15.12f %15.12f\n", $basis->[0][0], $basis->[1][0], $basis->[2][0]);
        print $poscar_fh sprintf("%15.12f %15.12f %15.12f\n", $basis->[0][1], $basis->[1][1], $basis->[2][1]);
        print $poscar_fh sprintf("%15.12f %15.12f %15.12f\n", $basis->[0][2], $basis->[1][2], $basis->[2][2]);
        print $poscar_fh "_VASP SPECIFIC ATOM LABELS_\n"; # following KISS, will not generate those lines
        print $poscar_fh "_VASP SPECIFIC ATOM COUNTS_\n"; # ...
        print $poscar_fh "Cartesian\n";

        for( my $j = 0; $j < $natoms; $j++)
        {
            my $sqrtm = sqrt($a_masses[$j]);
            my($dx, $dy, $dz) = ($disps[3*$j]*$qi0*$_/$sqrtm, $disps[3*$j+1]*$qi0*$_/$sqrtm, $disps[3*$j+2]*$qi0*$_/$sqrtm);
            print $poscar_fh sprintf("%15.12f %15.12f %15.12f\n", $coord_cart->[$j][0]+$dx, $coord_cart->[$j][1]+$dy, $coord_cart->[$j][2]+$dz);
        }
        print $poscar_fh "\n";
    }
}

print "DISPCAR created\n";
close($poscar_fh);

sub trim{ my $s=shift; $s =~ s/^\s+|\s+$//g; return $s;}

# http://stackoverflow.com/questions/439647/how-do-i-print-unique-elements-in-perl-array
# sub uniq {local %_; grep {!$_{$_}++} @_}

sub dirkar {
    my $vector = shift;
    my $basis = shift;
    my $total_atoms = shift;
    my ($i,$v1,$v2,$v3);

    for ($i=0; $i<$total_atoms; $i++) {
        $v1 = $vector->[$i][0]*$basis->[0][0] + $vector->[$i][1]*$basis->[0][1] + $vector->[$i][2]*$basis->[0][2];
        $v2 = $vector->[$i][0]*$basis->[1][0] + $vector->[$i][1]*$basis->[1][1] + $vector->[$i][2]*$basis->[1][2];
        $v3 = $vector->[$i][0]*$basis->[2][0] + $vector->[$i][1]*$basis->[2][1] + $vector->[$i][2]*$basis->[2][2];
        $vector->[$i][0] = $v1;
        $vector->[$i][1] = $v2;
        $vector->[$i][2] = $v3;
    }

    return ($vector);
}
