#!/usr/bin/env perl

# Copyright (c) "2012, Alexandr Fonari
#                URL: https://github.com/alexandr-fonari/Main/tree/master/VASP
#                License: MIT License
#
# ===== ## =====
# This script is a fork of (copyrighted accordingly):
#
#   Program: xdatgen
#   Sung Sakong, PhD
#   sung.sakong _at_ uni-due.de
#   ver 0.5  4. June 2007
#   url: http://www.uni-due.de/~hp0058/vmdplugins/utilities/xdatgen.c
###
#
#   Crystallographic conversions are heavily based on and copyrighted accordingly (some SUBs are copy/pasted):
#   http://theory.cm.utexas.edu/vasp/scripts/src/Vasp.pm
# ===== ## =====
##
#
#
my @displacements = (-1.0); # hard-coded so far for the five-point stencil 1st deriv.
#
##  ============== EDIT BELOW WITH CAUTION!! ============

use strict;
use warnings;
use Data::Dumper;

# physical constants in eV, Ang and s
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
use constant A2B => 1.88972613289; # ANGSTROM_TO_BOHR
use constant CmToEv => VaspToEv/VaspToCm;

die( "\nUse: $0 <input.out> <dyn.out>\n" ) if ( scalar( @ARGV ) < 2 );

open( my $outcar_fh, "<", $ARGV[0]) || die "Can't open $ARGV[0]: $!\n";

my ($alat, $basis, $coord_frac, $natoms, @a_labels, @a_masses, %label_mass, $coord_cart, @r);
while(my $line = <$outcar_fh>)
{
    if($line =~ /lattice parameter \(alat\)\s+=\s+([\d\.]+)/)
    {
        $alat = $1;
        printf("Found alat = %6.5f au\n", $alat)
    }

    if($line =~ /crystal axes: \(cart\. coord\. in units of alat\)/)
    {
        for (my $i=0; $i<3; $i++)
        {
                                    # a(1) = (   1.000000   0.000000   0.000000 )
            my @t = (<$outcar_fh> =~ m/.+?(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)/);
            # This is how VASP reads in the basis
            for (my $j=0; $j<3; $j++)
            {
                $basis->[$j][$i] = $t[$j]*$alat/A2B;
            }
            $r[$i] = sqrt($basis->[0][$i]**2 + $basis->[1][$i]**2 + $basis->[2][$i]**2)
        }
    }

    if($line =~ /atomic species\s+valence\s+mass\s+pseudopotential/)
    {
        while(my $next = <$outcar_fh>)
        {
            # C              4.00    12.01100      C( 1.00)
            if($next  =~ m/^\s+(\w+)\s+(\d+\.\d+)\s+(\d+\.\d+)/)
            {
                $label_mass{$1} = $3;
                print "Found atom type $1 with mass $3\n";
            }else{last;}
        }
    }

    if($line =~ /Crystallographic axes/)
    {
        my $next = <$outcar_fh>; # empty line
        $next = <$outcar_fh>; # site n.     atom                  positions (alat units)

        while($next = <$outcar_fh>)
        {
            # 1           C   tau(   1) = (  -0.1047125   0.2990902   0.0108898  )
            if($next  =~ m/^\s*\d+\s+(\w+).+?(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)/)
            {
                push(@a_labels, $1);
                push(@a_masses, $label_mass{$1});
                $coord_frac->[$#a_labels][0] = $2; # x
                $coord_frac->[$#a_labels][1] = $3; # y
                $coord_frac->[$#a_labels][2] = $4; # z
            }else{last;}
        }
        $natoms = scalar(@a_labels);
        print "Found $natoms atoms\n";
    }
}

$coord_cart = dirkar($coord_frac, $basis, $natoms);

open( my $dyncar_fh, "<", $ARGV[1]) || die "Can't open $ARGV[1]: $!\n";

my (@q, @e_values, @e_vectors);
while(my $line = <$dyncar_fh>)
{
    # q =       0.0000      0.0000     -0.6802
    if($line =~ /\s*q\s+=\s+([-\d\.]+)\s+([-\d\.]+)\s+([-\d\.]+)/) # its Cartesian (1/length)
    {
        @q = ($1, $2, $3);
        #print "Found q point: $qx, $qy, $qz 2pi/alat\n";
        @q = map { $_ *  2*PI/($alat/A2B) } @q;
        #print "Found q point: $qx, $qy, $qz 1/A\n";

    }

    # omega( 1) =       1.939432 [THz] =      64.692486 [cm-1]
    if($line =~ /\s*omega.+?([-\d\.]+)\s*\[cm-1\]/)
    {
        push(@e_values, $1);

        my $vec;
        for (my $i=0; $i<$natoms; $i++)
        {
            my ($xr, $xi, $yr, $yi, $zr, $zi) = (<$dyncar_fh> =~ m/(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)/);
            my @ph = (atan2($xi, $xr), atan2($yi, $yr), atan2($zi, $zr)); # phases
            my @f = (d0z(2*PI,$q[0]*$r[0]), d0z(2*PI,$q[1]*$r[1]), d0z(2*PI,$q[2]*$r[2]));

            #my $tz_i =(sin($phz+$qz*$r[2]*$fz)-sin($phz))/2.0; this is THE CHECK!
            #print $s[2]."\n";

            my @s;
            for (my $j=0; $j<3; $j++){
                $s[$j] = (cos($ph[$j]+$q[$j]*$r[$j]*$f[$j])+cos($ph[$j]))/2.0;
            }

            $vec .= sprintf("%10.6f %10.6f %10.6f", abs($xr)*$s[0], abs($yr)*$s[1], abs($zr)*$s[2])." ";
            #print sprintf("%10.6f %10.6f %10.6f", abs($xr)*$s[0], abs($yr)*$s[1], abs($zr)*$s[2])."\n";
            # print get_phase($t[0])." ".get_phase($t[2])." ".get_phase($t[4])."\n";
            #die;
        }
        push(@e_vectors, $vec);
    }
}

open( my $poscar_fh, ">", "DISPCAR" ) || die "Can't open DISPCAR file: $!";

## the differences between QE and VASP:
#  1. frequencies are positive and given in cm-1
#  2. eigenvectors are already normalized by sqrt(mass)
##
for( my $i = 0; $i < scalar(@e_values); $i++)
{
    print "processing ".($i+1)." out ".scalar(@e_values)." eigenvalues\n";
    my $ev = $e_values[$i];
    if($ev < 0.0){next;} # skip imaginary frequency

    my $qi0 = sqrt((HBAR*CL)**2/(AM*$ev*CmToEv)); # mode quanta
    # printf("%15.12f %15.12f\n", sqrt($ev)*VaspToEv, $qi0);

    my @disps = split('\s+', trim($e_vectors[$i]));

    foreach (@displacements)
    {
        print $poscar_fh sprintf("QE flavored POSCAR: disp=%f, hw=%8.5f meV\n", $_, $ev*CmToEv*1000);
        print $poscar_fh "1.00000\n";
        print $poscar_fh sprintf("%10.6f %10.6f %10.6f\n", $basis->[0][0], $basis->[1][0], $basis->[2][0]);
        print $poscar_fh sprintf("%10.6f %10.6f %10.6f\n", $basis->[0][1], $basis->[1][1], $basis->[2][1]);
        print $poscar_fh sprintf("%10.6f %10.6f %10.6f\n", $basis->[0][2], $basis->[1][2], $basis->[2][2]);
        print $poscar_fh "Cartesian\n";

        for( my $j = 0; $j < $natoms; $j++)
        {
            my $sqrtm = sqrt($a_masses[$j]);
            my($dx, $dy, $dz) = ($disps[3*$j]*$qi0*$_/$sqrtm, $disps[3*$j+1]*$qi0*$_/$sqrtm, $disps[3*$j+2]*$qi0*$_/$sqrtm);
            print $poscar_fh sprintf("%s %10.7f %10.7f %10.7f\n", $a_labels[$j], $coord_cart->[$j][0]+$dx, $coord_cart->[$j][1]+$dy, $coord_cart->[$j][2]+$dz);
        }
        print $poscar_fh "\n";
    }
}

print "DISPCAR created\n";
close($poscar_fh);

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

sub dot_product {
    my $coordinates1 = shift;
    my $coordinates2 = shift;
    my $total_atoms = shift;

    my ($i,$j);
    my $mag = 0;

    for ($i=0; $i<$total_atoms; $i++) {
        for ($j=0; $j<3; $j++) {
            $mag += $coordinates1->[$i][$j]*$coordinates2->[$i][$j];
        }
    }
    return ($mag);
}

sub trim{my $s=shift; $s =~ s/^\s+|\s+$//g; return $s;}
#sub rad2deg{my $rad=shift; return ($rad/PI)*180.0;}
sub d0z{my ($y, $x)=@_; if(!$x){return 0}else{return $y/$x};} # divide or return 0
# http://stackoverflow.com/questions/439647/how-do-i-print-unique-elements-in-perl-array
sub uniq {local %_; grep {!$_{$_}++} @_}
