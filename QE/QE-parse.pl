#!/usr/bin/env perl

# Copyright (c) "2012, Alexandr Fonari
#                URL: https://github.com/alexandr-fonari/Main/tree/master/VASP
#                License: MIT License

use strict;
use warnings;
use XML::Simple;
use Data::Dumper;
use Math::Complex;

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

die( "\nUse: $0 <input.out> <dyn.out>\n" ) if ( scalar( @ARGV ) < 2 );

open( my $outcar_fh, "<", $ARGV[0]) || die "Can't open $ARGV[0]: $!\n";

my ($alat, $basis, $coord_frac, $natoms, @a_labels, @a_masses, %label_mass);
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
            # This is how QE reads in the basis
            for (my $j=0; $j<3; $j++)
            {
                $basis->[$i][$j] = $t[$j]*$alat/A2B;
            }
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

        my $c = 0;
        while($next = <$outcar_fh>)
        {
            # 1           C   tau(   1) = (  -0.1047125   0.2990902   0.0108898  )
            if($next  =~ m/^\s*\d+\s+(\w+).+?(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)/)
            {
                push(@a_labels, $1);
                push(@a_masses, $label_mass{$1});
                $coord_frac->[$c][0] = $2; # x
                $coord_frac->[$c][1] = $3; # y
                $coord_frac->[$c][2] = $4; # z
                $c += 1;
            }else{last;}
        }
        $natoms = scalar(@a_labels);
        print "Found $natoms atoms\n";
    }
}

open( my $dyncar_fh, "<", $ARGV[1]) || die "Can't open $ARGV[1]: $!\n";

my (@eigen_vals, @eigen_vecs);
while(my $line = <$dyncar_fh>)
{
    # omega( 1) =       1.939432 [THz] =      64.692486 [cm-1]
    if($line =~ /\s*omega.+?([\d\.]+)\s*\[cm-1\]/)
    {
        push(@eigen_vals, $1);

        my $vec;
        for (my $i=0; $i<$natoms; $i++)
        {
            my @t = (<$dyncar_fh> =~ m/(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)\s+(-*\d+\.\d+)/);
            my $xi = Math::Complex->make($t[0], $t[3]);
            my $yi = Math::Complex->make($t[1], $t[4]);
            my $zi = Math::Complex->make($t[2], $t[5]);

            $vec .= sprintf("%10.8e %10.8e %10.8e", $xi->abs(), $yi->abs(), $zi->abs())." ";
        }
        push(@eigen_vecs, $vec);
    }
}

## creating XML file
my $xml_labels_indx;
for (my $i=0; $i<$natoms; $i++)
{
    $xml_labels_indx .= "<rc><c>".$a_labels[$i]." </c><c>   1</c></rc>\n";
}
my $xml_labels_masses;
foreach my $key (keys %label_mass)
{
    $xml_labels_masses .= "<rc><c>1</c><c> ".$key." </c><c> ".$label_mass{$key}." </c><c>4.0</c><c> GT rules </c></rc>\n";
}
my $xml_basis;
for (my $i=0; $i<3; $i++)
{
    $xml_basis .= "<v> ";
    for (my $j=0; $j<3; $j++)
    {
        $xml_basis .= sprintf("%15.12f", $basis->[$j][$i])." ";
    }
    $xml_basis .= "</v>\n";
}

my $xml_positions;
for (my $i=0; $i<$natoms; $i++)
{
    $xml_positions .= "<v> ".sprintf("%15.12f %15.12f %15.12f", $coord_frac->[$i][0], $coord_frac->[$i][1], $coord_frac->[$i][2])." </v>\n";
}

my $xml_eigen_vals = join(' ', @eigen_vals);

my $xml_eigen_vecs;
for (my $i=0; $i<scalar(@eigen_vals); $i++)
{
    $xml_eigen_vecs .= "<v> ".$eigen_vecs[$i]." </v>\n";
}
my $xml_out = <<END;
<?xml version="1.0" encoding="ISO-8859-1"?>
<modeling>
 <atominfo>
  <array name="atoms" >
   <dimension dim="1">ion</dimension>
   <field type="string">element</field>
   <field type="int">atomtype</field>
   <set>
$xml_labels_indx</set>
  </array>
  <array name="atomtypes" >
   <dimension dim="1">type</dimension>
   <field type="int">atomspertype</field>
   <field type="string">element</field>
   <field>mass</field>
   <field>valence</field>
   <field type="string">pseudopotential</field>
   <set>
$xml_labels_masses</set>
  </array>
 </atominfo>
 <structure name="initialpos" >
  <crystal>
   <varray name="basis" >
$xml_basis</varray>
  </crystal>
  <varray name="positions" >
$xml_positions</varray>
 </structure>
 <calculation>
  <dynmat>
   <v name="eigenvalues">$xml_eigen_vals</v>
   <varray name="eigenvectors" >
$xml_eigen_vecs </varray>
  </dynmat>
 </calculation>
</modeling>
END

print $xml_out;

sub trim{ my $s=shift; $s =~ s/^\s+|\s+$//g; return $s;}

# http://stackoverflow.com/questions/439647/how-do-i-print-unique-elements-in-perl-array
sub uniq {local %_; grep {!$_{$_}++} @_}
