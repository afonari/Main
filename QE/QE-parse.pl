#!/usr/bin/env perl

# Copyright (c) "2012, Alexandr Fonari
#                URL: https://github.com/alexandr-fonari/Main/tree/master/VASP
#                License: MIT License

use strict;
use warnings;
use XML::Simple;
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
my $xml = new XML::Simple();
my $data = $xml->XMLin("vasprun.xml");

# atomic masses and labels
my (@a_labels, @a_indx, @a_masses, @a_count);

my @t = @{$data->{"atominfo"}->{"array"}->{"atoms"}->{"set"}->{"rc"}};
foreach (@t)
{
    #push(@a_masses, (trim($_->{"c"}[2])) x trim($_->{"c"}[0]));
    push(@a_labels, (trim($_->{"c"}[0]))); # x trim($_->{"c"}[0]));
    push(@a_indx,   (trim($_->{"c"}[1]))); # x trim($_->{"c"}[0]));
    #push(@a_count, trim($_->{"c"}[0]));
}
# print Dumper( @a_labels, @a_indx );

@t = @{$data->{"atominfo"}->{"array"}->{"atomtypes"}->{"set"}->{"rc"}};
my %label_mass;
foreach (@t)
{
    my ($label, $mass) = ( trim($_->{"c"}[1]), trim($_->{"c"}[2]) );
    $label_mass{$label} = $mass;
}
# print Dumper( \% label_mass);

foreach my $label (@a_labels)
{
    push( @a_masses, $label_mass{$label} );
}
print Dumper( @a_masses );
die;
# basis vectors
my @f;
@t = @{$data->{"structure"}->{"initialpos"}->{"crystal"}->{"varray"}->{"basis"}->{"v"}};
($f[1], $f[2], $f[3]) = split('\s+', trim($t[0]));
($f[4], $f[5], $f[6]) = split('\s+', trim($t[1]));
($f[7], $f[8], $f[9]) = split('\s+', trim($t[2]));
# print Dumper( @f );

# fractional atomic positions
my (@a_frac_pos_x, @a_frac_pos_y, @a_frac_pos_z);
@t = @{$data->{"structure"}->{"initialpos"}->{"varray"}->{"v"}};
foreach (@t)
{
    my @t = split('\s+', trim($_));
    push(@a_frac_pos_x, $t[0]);
    push(@a_frac_pos_y, $t[1]);
    push(@a_frac_pos_z, $t[2]);
}
# print Dumper(@a_frac_pos_x, @a_frac_pos_y, @a_frac_pos_z);

# obtain Cartesian coordinates: Transpose(basis)*v:
my (@a_cart_pos_x, @a_cart_pos_y, @a_cart_pos_z);
for(my $i = 0; $i < scalar(@a_frac_pos_x); $i++ )
{
    my @v = ($a_frac_pos_x[$i], $a_frac_pos_y[$i], $a_frac_pos_z[$i]);
    push(@a_cart_pos_x, $f[1]*$v[0] + $f[4]*$v[1] + $f[7]*$v[2]);
    push(@a_cart_pos_y, $f[2]*$v[0] + $f[5]*$v[1] + $f[8]*$v[2]);
    push(@a_cart_pos_z, $f[3]*$v[0] + $f[6]*$v[1] + $f[9]*$v[2]);
}
# print Dumper(@a_cart_pos_x, @a_cart_pos_y, @a_cart_pos_z);

# hessian eigenvalues
my @e_values = split('\s+', trim($data->{"calculation"}->{"dynmat"}->{"v"}->{"content"}));

# hessian eigenvectors
my @e_vectors = @{$data->{"calculation"}->{"dynmat"}->{"varray"}->{"eigenvectors"}->{"v"}};
# print Dumper(@e_vectors);

open( my $poscar_fh, ">", "DISPCAR" ) || die "Can't open DISPCAR file: $!";

my @displacements = (-2, -1, 0, 1, 2); # hard-coded so far for 5-point stencil finite difference 1st deriv.
for( my $i = 0; $i < scalar(@e_values); $i++)
{
    print "processing ".($i+1)." out ".scalar(@e_values)." eigenvalues\n";
    my $ev = $e_values[$i]*(-1.0);
    if($ev < 0.0){next;} # skip imaginary frequency

    my $qi0 = sqrt((HBAR*CL)**2/(AM*sqrt($ev)*VaspToEv)); # mode quanta
    # printf("%15.12f %15.12f\n", sqrt($ev)*VaspToEv, $qi0);

    my @disps = split('\s+', trim($e_vectors[$i]));

    foreach (@displacements)
    {
        print $poscar_fh sprintf("POSCAR: disp=%d, w=%8.5f meV, qi0=%5e\n", $_, sqrt($ev)*VaspToEv*1000, $qi0);
        print $poscar_fh "1.00000\n";
        print $poscar_fh sprintf("%15.12f %15.12f %15.12f\n", $f[1], $f[2], $f[3]);
        print $poscar_fh sprintf("%15.12f %15.12f %15.12f\n", $f[4], $f[5], $f[6]);
        print $poscar_fh sprintf("%15.12f %15.12f %15.12f\n", $f[7], $f[8], $f[9]);
        print $poscar_fh join(" ", uniq(@a_labels))."\n";
        print $poscar_fh join(" ", @a_count)."\n";
        print $poscar_fh "Cartesian\n";

        for( my $j = 0; $j < scalar(@a_cart_pos_x); $j++)
        {
            my $sqrtm = sqrt($a_masses[$j]);
            my($dx, $dy, $dz) = ($disps[3*$j]*$qi0*$_/$sqrtm, $disps[3*$j+1]*$qi0*$_/$sqrtm, $disps[3*$j+2]*$qi0*$_/$sqrtm);
            print $poscar_fh sprintf("%15.12f %15.12f %15.12f\n", $a_cart_pos_x[$j]+$dx, $a_cart_pos_y[$j]+$dy, $a_cart_pos_z[$j]+$dz);
        }
        print $poscar_fh "\n";
    }
    
}

print "DISPCAR created\n";
close($poscar_fh);

sub trim{ my $s=shift; $s =~ s/^\s+|\s+$//g; return $s;}

# http://stackoverflow.com/questions/439647/how-do-i-print-unique-elements-in-perl-array
sub uniq {local %_; grep {!$_{$_}++} @_}
