#!/usr/bin/env perl

# Copyright (c) "2012, Alexandr Fonari
#                URL: https://github.com/alexandr-fonari/Main/tree/master/CRYSTAL
#                License: MIT License
# Version 1.0
#
# Program converts CRYSTAL output into CIF file.
# Output from ionic relaxation is supported (variable unit cell relaxations have not been tested),
# and will result in CIF with many structures (as many as ionic steps)
# CIF file can be opened with CSD Mercury (c).
#


use strict; 
use warnings;
use Data::Dumper;
use Math::Trig; # for deg2rad
use Math::Complex; # for Re function

use constant Bohr => 0.52917721092;

if ( scalar( @ARGV ) < 1 )
{
    die( "\nUse: $0 <scf.out>\n" );
}

open( my $outcar_fh, "<", $ARGV[0]) || die "Can't open $ARGV[0]: $!\n";
open( my $cif_fh, ">", "$ARGV[0].cif" ) || die "Can't open $ARGV[0].cif: $!\n";

my $count = 0;
my @atoms_lines;

my ($alat);
while(my $line = <$outcar_fh>)
{
    if($line =~ /lattice parameter \(alat\)\s+=\s+([\d\.]+)/)
    {
        $alat = $1;
        printf("Found alat = %6.5f au\n", $alat)
    }

    if($line =~ /crystal axes: \(cart\. coord\. in units of alat\)/)
    {
        my @f;
        ($f[1], $f[2], $f[3]) = (<$outcar_fh> =~ m/.+?(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)/);
        ($f[4], $f[5], $f[6]) = (<$outcar_fh> =~ m/.+?(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)/);
        ($f[7], $f[8], $f[9]) = (<$outcar_fh> =~ m/.+?(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)/);
        foreach my $x (@f[1..9]) { $x = $x * $alat * Bohr; }

        # unit cell dimensions
        my $a = sqrt($f[1]**2 + $f[2]**2 + $f[3]**2);
        my $b = sqrt($f[4]**2 + $f[5]**2 + $f[6]**2);
        my $c = sqrt($f[7]**2 + $f[8]**2 + $f[9]**2);
        my $alpha = Re(rad2deg(acos(($f[4]*$f[7] + $f[5]*$f[8] + $f[6]*$f[9])/($b*$c))));
        my $beta = Re(rad2deg(acos(($f[1]*$f[7] + $f[2]*$f[8] + $f[3]*$f[9])/($a*$c))));
        my $gamma = Re(rad2deg(acos(($f[1]*$f[4] + $f[2]*$f[5] + $f[3]*$f[6])/($a*$b))));

        print_cif_header($cif_fh, $a, $b, $c, $alpha, $beta, $gamma);
    }

    if($line =~ /Crystallographic axes/)
    {
        my $next = <$outcar_fh>; # empty line
        $next = <$outcar_fh>; # site n.     atom                  positions (alat units)

        while($next = <$outcar_fh>)
        {
            # 1           C   tau(   1) = (  -0.1047125   0.2990902   0.0108898  )
            if($next  =~ m/^\s*(\d+)\s+(\w+).+?(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)/)
            {
                my ($index, $label, $x, $y, $z) = ($1, $2, $3, $4, $5);
                print $cif_fh sprintf("%s%d %s %9.5f %9.5f %9.5f\n", $label, $index, $label, $x, $y, $z);
            }else{last;}
        }
    }
}

close($outcar_fh);
close($cif_fh);

print "\nqe2cif.pl V. 1.0.\nFile $ARGV[0].cif created.\nJob done.\n\n";

# SUBs
sub print_cif_header
{
    my ($cif_fh, $a, $b, $c, $alpha, $beta, $gamma) = @_;
    print $cif_fh "\n\ndata_a\n";
    print $cif_fh "\n";
    print $cif_fh "_cell_length_a     $a\n";
    print $cif_fh "_cell_length_b     $b\n";
    print $cif_fh "_cell_length_c     $c\n";
    print $cif_fh "_cell_angle_alpha  $alpha\n";
    print $cif_fh "_cell_angle_beta   $beta\n";
    print $cif_fh "_cell_angle_gamma  $gamma\n";
    print $cif_fh "\n";
    print $cif_fh "loop_\n";
    print $cif_fh "_space_group_symop_id\n";
    print $cif_fh "_space_group_symop_operation_xyz\n";
    print $cif_fh "1 'x, y, z'\n";
    print $cif_fh "\n";
    print $cif_fh "loop_\n";
    print $cif_fh "_atom_site_label\n";
    print $cif_fh "_atom_site_type_symbol\n";
    print $cif_fh "_atom_site_fract_x\n";
    print $cif_fh "_atom_site_fract_y\n";
    print $cif_fh "_atom_site_fract_z\n";
}
sub trim{ my $s=shift; $s =~ s/^\s+|\s+$//g; return $s;}
