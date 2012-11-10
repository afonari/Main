#!/usr/bin/env perl

use strict;
use warnings;
use XML::Simple;
use Data::Dumper;


# physical constants in eV and Ang
use constant HBAR => 6.58211915e-16;
use constant CL => 2.99792458e18;
use constant AM => 931.494043e6;

use constant Angstrom => 1.0e-10;  # [m]
use constant EV => 1.60217733e-19; # [J]
use constant AMU => 1.6605402e-27; # [kg]
use constant PI    => 4 * atan2(1, 1);
use constant PlanckConstant => 4.13566733e-15; # [eV s]
use constant VaspToEv => sqrt(EV/AMU)/Angstrom/(2*PI)*PlanckConstant; # [eV] 6.46541380e-2

my $xml = new XML::Simple();
my $data = $xml->XMLin("vasprun.xml");

# getting atomic masses
my %a_masses;
my @t = @{$data->{"atominfo"}->{"array"}->{"atomtypes"}->{"set"}->{"rc"}};
foreach (@t)
{
    $a_masses{trim($_->{"c"}[1])} = trim($_->{"c"}[2]);
}
# print Dumper( \%a_masses );

my @a;
@t = @{$data->{"atominfo"}->{"array"}->{"atoms"}->{"set"}->{"rc"}};
foreach (@t)
{
    push(@a, $a_masses{trim($_->{"c"}[1])});
}
# print Dumper( \%a_masses );

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
#print Dumper(@a_frac_pos_x);
#print Dumper(@a_frac_pos_y);
#print Dumper(@a_frac_pos_z);

# obtain Cartesian coordinates: Transpose(basis)*v:
my (@a_cart_pos_x, @a_cart_pos_y, @a_cart_pos_z);
for(my $i = 0; $i < scalar(@a_frac_pos_x); $i++ )
{
    my @v = ($a_frac_pos_x[$i], $a_frac_pos_y[$i], $a_frac_pos_z[$i]);
    push(@a_cart_pos_x, $f[1]*$v[0] + $f[4]*$v[1] + $f[7]*$v[2]);
    push(@a_cart_pos_y, $f[2]*$v[0] + $f[5]*$v[1] + $f[8]*$v[2]);
    push(@a_cart_pos_z, $f[3]*$v[0] + $f[6]*$v[1] + $f[9]*$v[2]);
}
#print Dumper(@a_cart_pos_x);
#print Dumper(@a_cart_pos_y);
#print Dumper(@a_cart_pos_z);

# hessian eigenvalues
my @e_values = split('\s+', trim($data->{"calculation"}->{"dynmat"}->{"v"}->{"content"}));

# hessian eigenvectors
my @e_vectors = @{$data->{"calculation"}->{"dynmat"}->{"varray"}->{"eigenvectors"}->{"v"}};
#print Dumper(@e_vectors);

my @disp = (-1, 1);
for( my $i = 0; $i < scalar(@e_values); $i++)
{
    my $ev = $e_values[$i]*(-1.0);
    if($ev < 0.0){next;} # skip imaginary frequency

    my $qi0 = sqrt((HBAR*CL)**2/(AM*$ev*VaspToEv)); # a quanta
    my @disps = split('\s+', trim($e_vectors[$i]));

    for( my $j = 0; $j < scalar(@a_cart_pos_x); $j++)
    {
        #printf("%10.6f %10.6f %10.6f %10.6f\n", $disps[3*$j], $disps[3*$j+1], $disps[3*$j+2], $qi0);
    }
    #print "\n";
    
}
#foreach my $attributes (keys %{$data}){
#  print"$attributes : ${$data}{$attributes}\n";
#}

sub trim
{
   my $string = shift;
   $string =~ s/^\s+|\s+$//g;
   return $string;
}

