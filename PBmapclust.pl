##########################    In script documentation	##########################
=head

Description :
This script is the main program for PBmapclust, it takes GROMACS run input file and trajectory files (.tpr & .trr/.xtc) and generates a MD.pb file with PB sequences corresponding to the frames from MD simulation.
It is a wrapper program for running the other scripts- dssp_separator.pl, dssp_to_pb.pl and pb2xpm.pl.
It requires DSSP, GROMACS and Weblogo programs to be installed.

Switch comments (#) for choosing options in the script.

Usage : example- perl PBmapclust.pl MD.tpr MD.trr

=cut
##################################################################################


#!/usr/bin/perl

@steps=`gmx dump -s $ARGV[0] -quiet| grep -E "(nstxout|dt|nstep)"`;
@line1=split /=\s/,$steps[0];
@line2=split /=\s/,$steps[1];
@line3=split /=\s/,$steps[2];
@line4=split /=\s/,$steps[3];
$dt= $line1[1];
$nsteps= $line2[1];
$nstxout= $line3[1];
$nstxout_compressed= $line4[1];

$MD_duration_ps=$nsteps*$dt;
$num_frames=$nsteps/$nstxout;
$one_frame_every_ps=$MD_duration_ps/$num_frames;
$skip=10;#Number of frames to skip from the MD trajectory between generating PB sequences.
print"Analysing $num_frames frames from a $MD_duration_ps (ps) simulation\nPlease wait...\n";
for($frame=0;$frame<=$MD_duration_ps;$frame+=($one_frame_every_ps*$skip))
{
	if($frame==0)
	{
		`echo 1|gmx trjconv -s $ARGV[0] -f $ARGV[1] -dump $frame -o $frame.pdb -quiet >/dev/null 2>&1`;
		`dssp $frame.pdb $frame.pdb.dssp`;
		open(DSSP,"<","0.pdb.dssp")||print"Cannot open the file 0.pdb.dssp\n";
		@dssp=<DSSP>;
		close(DSSP);	
		@check_chain_present=split(//,$dssp[$#dssp]);		
		if($check_chain_present[11] eq " ")
		{
			push(@chains,'A');
			`rm *.pdb`;
			`mv $frame.pdb.dssp $frame.pdb_A.dssp`;
			$check_chain_flag=1;
		}
		else	
		{
			`perl dssp_separator.pl $frame.pdb.dssp ./`;
			`rm *.pdb`;
			`rm *.pdb.dssp`;
			@chains=`ls 0.pdb_*.dssp`;		
			foreach (@chains)
			{
				chomp $_;
				$_=~ s/0.pdb_//;
				$_=~ s/.dssp//;
			}
		}
		opendir DIR, "./" or die "cannot open dir $dir: $!";
		@dssp_files_seperated = grep { -f && /\.dssp$/ } readdir DIR;
		for $dssp_chain (@dssp_files_seperated)
		{
			`perl dssp_to_pb.pl $dssp_chain ./`;
		}
		`rm *.dssp`;
		foreach $chain (@chains)
		{
			open(FOUT,">","$num_frames\_$chain.pb")||die "Cannot open the file $num_frames\_$chain.pb\n";
			open(FIN,"<","$frame.pdb_$chain.dssp.pb")||die "Cannot open the file $frame.pdb_$chain.dssp.pb\n";
			$pb_seq=<FIN>;
			close FIN;
			print FOUT ">","$frame.pdb_$chain.dssp.pb","\n","$pb_seq\n";
		}
		`rm *.dssp.pb`;
	}
	else
	{
		$frame_minus1=$frame-1;
		`echo 1|gmx trjconv -s $ARGV[0] -f $ARGV[1] -b $frame_minus1 -dump $frame -o $frame.pdb -skip $skip -quiet >/dev/null 2>&1`;
		`dssp $frame.pdb $frame.pdb.dssp`;
		if($check_chain_flag==1)
		{
			`rm *.pdb`;
			`mv $frame.pdb.dssp $frame.pdb_A.dssp`;
		}
		else	
		{
		`perl dssp_separator.pl $frame.pdb.dssp ./`;
		`rm *.pdb`;
		`rm *.pdb.dssp`;
		}		
		opendir DIR, "./" or die "cannot open dir $dir: $!";
		@dssp_files_seperated = grep { -f && /\.dssp$/ } readdir DIR;
		for $dssp_chain (@dssp_files_seperated)
		{
			`perl dssp_to_pb.pl $dssp_chain ./`;
		}
		`rm *.dssp`;
		foreach $chain (@chains)
		{
			open(FOUT,">>","$num_frames\_$chain.pb")||die "Cannot open the file $num_frames\_$chain.pb\n";
			open(FIN,"<","$frame.pdb_$chain.dssp.pb")||die "Cannot open the file $frame.pdb_$chain.dssp.pb\n";
			$pb_seq=<FIN>;
			close FIN;
			print FOUT ">","$frame.pdb_$chain.dssp.pb","\n","$pb_seq\n";
		}
		`rm *.dssp.pb`;
	}
close FOUT;	
}

#Switch comment to choose color_scheme
$color_scheme='propensity';
#$color_scheme='random';

foreach $chain (@chains)
{
	#Switch comment to run pb2xpm.pl with options 
	#`perl pb2xpm.pl $num_frames\_$chain.pb 1500 2000 1 225 $color_scheme >$num_frames\_$chain.xpm`;
	`perl pb2xpm.pl $num_frames\_$chain.pb $color_scheme >$num_frames\_$chain.xpm`;
	
        `gmx xpm2ps -f $num_frames\_$chain.xpm -o $num_frames\_$chain.eps -di PBmapclust.m2p >/dev/null 2>&1`;
	if($color_scheme eq 'random')
	{
		`weblogo -f $num_frames\_$chain.pb -D fasta -o $num_frames\_$chain.png -F png_print --units probability -a ABCDEFGHIJKLMNOPZX --color '#FFFFFF' X 'unassigned' --color '#000000' Z 'first-last' --color '#000000' A 'A' --color '#FA58F4' B 'B' --color '#00561B' C 'C' --color '#FFFF00' D 'D' --color '#CECECE' E 'E' --color '#0000FF' F 'F' --color '#9D3E0C' G 'G' --color '#9E0E40' H 'H' --color '#FE9A2E' I 'I' --color '#00FFFF' J 'J' --color '#00FF00' K 'K' --color '#0B0B61' L 'L' --color '#FF0000' M 'M' --color '#C8AD7F' N 'N' --color '#B0F2B6' O 'O' --color '#606060' P 'P' --composition none -t "Position wise PB distribution during MD" -y 'PB Distribution' -n 100`;
	}
	elsif($color_scheme eq 'propensity')
	{
		`weblogo -f $num_frames\_$chain.pb -D fasta -o $num_frames\_$chain.png -F png_print --units probability -a ABCDEFGHIJKLMNOPZX --color '#FFFFFF' X 'unassigned' --color '#000000' Z 'first-last' --color '#00C43B' A 'A' --color '#00DE21' B 'B' --color '#00946B' C 'C' --color '#0047B8' D 'D' --color '#007F7F' E 'E' --color '#00B847' F 'F' --color '#12D617' G 'G' --color '#00CF30' H 'H' --color '#00F00F' I 'I' --color '#0AE014' J 'J' --color '#59A303' K 'K' --color '#708C03' L 'L' --color '#DE2100' M 'M' --color '#AD4F00' N 'N' --color '#6E9100' O 'O' --color '#1CE003' P 'P' --composition none -t "Position wise PB distribution during MD" -y 'PB Distribution' -n 100`;
	}
}

print"Done !\n";



