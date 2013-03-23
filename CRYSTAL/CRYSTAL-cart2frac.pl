#!/usr/bin/env perl
#;-*- Perl -*-

# Copyright (c) "2012, Alexandr Fonari
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

# from http://cpansearch.perl.org/src/BDFOY/Chemistry-Elements-1.07/Elements.pm
my %e=('H'=>'1','He'=>'2','Li'=>'3','Be'=>'4','B'=>'5','C'=>'6','N'=>'7','O'=>'8','F'=>'9','Ne'=>'10','Na'=>'11','Mg'=>'12','Al'=>'13','Si'=>'14','P'=>'15','S'=>'16','Cl'=>'17','Ar'=>'18','K'=>'19','Ca'=>'20','Sc'=>'21','Ti'=>'22','V'=>'23','Cr'=>'24','Mn'=>'25','Fe'=>'26','Co'=>'27','Ni'=>'28','Cu'=>'29','Zn'=>'30','Ga'=>'31','Ge'=>'32','As'=>'33','Se'=>'34','Br'=>'35','Kr'=>'36','Rb'=>'37','Sr'=>'38','Y'=>'39','Zr'=>'40','Nb'=>'41','Mo'=>'42','Tc'=>'43','Ru'=>'44','Rh'=>'45','Pd'=>'46','Ag'=>'47','Cd'=>'48','In'=>'49','Sn'=>'50','Sb'=>'51','Te'=>'52','I'=>'53','Xe'=>'54','Cs'=>'55','Ba'=>'56','La'=>'57','Ce'=>'58','Pr'=>'59','Nd'=>'60','Pm'=>'61','Sm'=>'62','Eu'=>'63','Gd'=>'64','Tb'=>'65','Dy'=>'66','Ho'=>'67','Er'=>'68','Tm'=>'69','Yb'=>'70','Lu'=>'71','Hf'=>'72','Ta'=>'73','W'=>'74','Re'=>'75','Os'=>'76','Ir'=>'77','Pt'=>'78','Au'=>'79','Hg'=>'80','Tl'=>'81','Pb'=>'82','Bi'=>'83','Po'=>'84','At'=>'85','Rn'=>'86','Fr'=>'87','Ra'=>'88','Ac'=>'89','Th'=>'90','Pa'=>'91','U'=>'92','Np'=>'93','Pu'=>'94','Am'=>'95','Cm'=>'96','Bk'=>'97','Cf'=>'98','Es'=>'99','Fm'=>'100','Md'=>'101','No'=>'102','Lr'=>'103','Rf'=>'104','Ha'=>'105','Sg'=>'106','Bh'=>'107','Hs'=>'108','Mt'=>'109');

if ( scalar( @ARGV ) < 2 )
{
    die( "\nUse: $0 <output.out> <atoms.list>\n" );
}
open( my $outcar_fh, "<", $ARGV[0]) || die "Can't open $ARGV[0]: $!\n";
open( my $f_fh, "<", $ARGV[1]) || die "Can't open $ARGV[1]: $!\n";

my (@f_labels, $f_atoms);
my $iter = 0;
while(my $line = <$f_fh>)
{
    next if $line =~ /^\s*$/;

    $line = trim($line);
    my ($label, @rest) = split /\s+/, $line;
    push(@f_labels, $label);

    $f_atoms->[$iter][0] = $rest[0];
    $f_atoms->[$iter][1] = $rest[1];
    $f_atoms->[$iter][2] = $rest[2];

    $iter++;
}
close($f_fh);

my $basis;
while(my $line = <$outcar_fh>)
{
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
    }
}
close($outcar_fh);

$f_atoms = kardir($f_atoms,$basis,scalar(@f_labels));

open( my $poscar_frac_fh, ">", $ARGV[1].".frac" ) || die "Can't open $ARGV[1].frac file: $!";
for( my $i = 0; $i < scalar(@f_labels); $i++ )
{
    my $label = $f_labels[$i]; $label =~ s/[0-9]//g;
    my $number = $e{$label};
    my @coords = ($f_atoms->[$i][0], $f_atoms->[$i][1], $f_atoms->[$i][2]);

    print $poscar_frac_fh sprintf("%s %9.5f %9.5f %9.5f\n", $number, @coords);
}

close($poscar_frac_fh);
print "$ARGV[1].frac file has been created.\n\n";

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
    my $total_atoms = shift;
    my $recip_basis;
    my ($v1,$v2,$v3,$i,$j);
  
    $recip_basis = inverse($basis);
  
    for ($i=0; $i<$total_atoms; $i++) {
        $v1 = $vector->[$i][0]*$recip_basis->[0][0] + $vector->[$i][1]*$recip_basis->[1][0] + $vector->[$i][2]*$recip_basis->[2][0];
        $v2 = $vector->[$i][0]*$recip_basis->[0][1] + $vector->[$i][1]*$recip_basis->[1][1] + $vector->[$i][2]*$recip_basis->[2][1];
        $v3 = $vector->[$i][0]*$recip_basis->[0][2] + $vector->[$i][1]*$recip_basis->[1][2] + $vector->[$i][2]*$recip_basis->[2][2];

        # move atoms to primative cell
        $vector->[$i][0] = $v1;#+60-int($v1+60);
        $vector->[$i][1] = $v2;#+60-int($v2+60);
        $vector->[$i][2] = $v3;#+60-int($v3+60);
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

