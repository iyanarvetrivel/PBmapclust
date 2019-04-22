##########################    In script documentation	##########################
=head

Description : 
This script reads a "MD.pb" file with multiple PB sequences in fasta format and generates a corresponding MD.xpm file.
This script can be run with options:
perl pb2xpm.pl MD.pb start_frame end_frame start_residue end_residue color_scheme >MD.xpm
or as default:
perl pb2xpm.pl MD.pb color_scheme >MD.xpm

Usage : example- 
perl pb2xpm.pl 2000.pb 1500 2000 1 225 propensity >1500-2000.xpm
or
perl pb2xpm.pl 2000.pb propensity >2000.xpm

Subsequently, the .xpm file can be converted to an .eps file using the xpm2ps program from GROMACS.
gmx xpm2ps -f 4000_3.xpm -o 4000_3.eps -di PBmapclust.m2p

=cut
##################################################################################

#!/usr/bin/perl


%PBvalue = () ;
open(FIN,"<","$ARGV[0]")||die "Cannot open the file $ARGV[0]\n";
@pb=<FIN>;
if(scalar (@ARGV)>2)
{
	$start_frame=$ARGV[1]*2-1;
	$end_frame=$ARGV[2]*2-1;
	$start_residue=$ARGV[3];
	$end_residue=$ARGV[4];
	$nb_frames=$ARGV[2]-$ARGV[1]+1;
	$nb_residues=$ARGV[4]-$ARGV[3]+1;
	chomp($ARGV[5]);
	$color_scheme=$ARGV[5];	
}
else
{
	chomp($pb[1]);
	$start_frame=1;
	$end_frame=scalar @pb-1;
	$start_residue=1;
	$end_residue=length($pb[1]);
	$nb_frames=scalar @pb/2;
	$nb_residues=length($pb[1]);
	chomp($ARGV[1]);	
	$color_scheme=$ARGV[1];
}

@x_axis=($start_frame..$nb_frames);
@y_axis=($start_residue..$end_residue);
print  "\/\* XPM \*\/\n" ;
print  "\/\* title:   \"", $ARGV[0], " PBmapclust plot\" \*\/\n";
print  "\/\* legend:  \"PB\" \*\/\n";
print  "\/\* x-label: \"Frames\" \*\/\n";
print  "\/\* y-label: \"Residue\" \*\/\n";
print  "\/\* type:    \"Discrete\" \*\/\n" ;
print  "static char \* PB_matrix\[\] = \{\n";
print  "\"$nb_frames ", $nb_residues, " 18 1\",\n";
if($color_scheme eq 'random')
{
print  "\"A c \#000000\"  \/\* \"A\" \*\/,\n";
print  "\"B c \#FA58F4\"  \/\* \"B\" \*\/,\n";
print  "\"C c \#00561B\"  \/\* \"C\" \*\/,\n";
print  "\"D c \#FFFF00\"  \/\* \"D\" \*\/,\n";
print  "\"E c \#CECECE\"  \/\* \"E\" \*\/,\n";
print  "\"F c \#0000FF\"  \/\* \"F\" \*\/,\n";
print  "\"G c \#9D3E0C\"  \/\* \"G\" \*\/,\n";
print  "\"H c \#9E0E40\"  \/\* \"H\" \*\/,\n";
print  "\"I c \#FE9A2E\"  \/\* \"I\" \*\/,\n";
print  "\"J c \#00FFFF\"  \/\* \"J\" \*\/,\n";
print  "\"K c \#00FF00\"  \/\* \"K\" \*\/,\n";
print  "\"L c \#0B0B61\"  \/\* \"L\" \*\/,\n";
print  "\"M c \#FF0000\"  \/\* \"M\" \*\/,\n";
print  "\"N c \#C8AD7F\"  \/\* \"N\" \*\/,\n";
print  "\"O c \#B0F2B6\"  \/\* \"O\" \*\/,\n";
print  "\"P c \#606060\"  \/\* \"P\" \*\/,\n";
print  "\"Z c \#000000\"  \/\* \"Z\" \*\/,\n";
print  "\"X c \#FFFFFF\"  \/\* \"X\" \*\/,\n";
}
elsif($color_scheme eq 'propensity')
{
print  "\"A c \#00C43B\"  \/\* \"A\" \*\/,\n";
print  "\"B c \#00DE21\"  \/\* \"B\" \*\/,\n";
print  "\"C c \#00946B\"  \/\* \"C\" \*\/,\n";
print  "\"D c \#0047B8\"  \/\* \"D\" \*\/,\n";
print  "\"E c \#007F7F\"  \/\* \"E\" \*\/,\n";
print  "\"F c \#00B847\"  \/\* \"F\" \*\/,\n";
print  "\"G c \#12D617\"  \/\* \"G\" \*\/,\n";
print  "\"H c \#00CF30\"  \/\* \"H\" \*\/,\n";
print  "\"I c \#00F00F\"  \/\* \"I\" \*\/,\n";
print  "\"J c \#0AE014\"  \/\* \"J\" \*\/,\n";
print  "\"K c \#59A303\"  \/\* \"K\" \*\/,\n";
print  "\"L c \#708C03\"  \/\* \"L\" \*\/,\n";
print  "\"M c \#DE2100\"  \/\* \"M\" \*\/,\n";
print  "\"N c \#AD4F00\"  \/\* \"N\" \*\/,\n";
print  "\"O c \#6E9100\"  \/\* \"O\" \*\/,\n";
print  "\"P c \#1CE003\"  \/\* \"P\" \*\/,\n";
print  "\"Z c \#000000\"  \/\* \"Z\" \*\/,\n";
print  "\"X c \#FFFFFF\"  \/\* \"X\" \*\/,\n";
}
$nb_iterations=int(scalar @x_axis/50);
$remainder=scalar @x_axis%50;
#print $remainder, "\n";
for ($iteration=0;$iteration<$nb_iterations;$iteration++)
{
	$start_x_axis=$iteration*50;
	$end_x_axis=($iteration+1)*50-1;
	print "\/* x-axis: @x_axis[$start_x_axis..$end_x_axis] \*\/,\n";
}

if ($remainder!=0)
{
	$start_x_axis_last=(scalar @x_axis*50)+1;
	$end_x_axis_last=(scalar @x_axis*50)+$remainder;
	print "\/* x-axis: @x_axis[$start_x_axis_last..$end_x_axis_last] \*\/,\n";
}

#print  "\/* x-axis: @x_axis \*\/,\n";
print  "\/* y-axis: @y_axis \*\/,\n";

for ($i=1;$i<=scalar @pb-1;$i+=2)
{
	chomp($pb[$i]) ;
	@PBseq = split (//,$pb[$i]) ;
	for ($j=0 ; $j <=scalar(@PBseq) ; $j++) 
	{
		$pos = $j + 1 ;
		$PBvalue{$i."frame,residue".$pos} = $PBseq[$j] ;
	}
}
for ($residue = $end_residue; $residue>=$start_residue; $residue--)
{
	print  "\"";
	for ($frame=$start_frame ; $frame <= $end_frame ; $frame+=2) {
		print  $PBvalue{$frame."frame,residue".$residue};
	}
	if ($residue!=$start_residue) {
		print  "\",\n";
	} 
	else {
		print  "\"";
	}
}
print  "\};";
