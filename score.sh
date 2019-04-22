#!/bin/bash

#THIS SCRIPT CALCULATES KIMURA DISTANCE BETWEEN ALL PAIRWISE DISTANCES BETWEEN PB SEQUENCES EXTRACTED FROM MD SIMULATION
# CONTACT DETAILS: bernard.offmann@univ-nantes.fr

#TO EXECUTE THIS SCRIPT, A TEXT FILE CONTAINTING PB SEQUENCES IN FASTA FORMAT IS REQUIRED
# EXAMPLE OF CONTENT OF PB SEQUENCE FILE :
#>0.pdb_A.dssp.pb
#ZZACDDDDFBDFKBCCDDDDDEDDEHIAFBDCFKBCCDFBDGHKLPCCFBDCDDDDDDFBACDFKLCCDDDFBGNHIACDDDDDFBFKLCCFBFBFBDCFBFKBCCDFKLGMGHPLMMMMZZ
#>1000.pdb_A.dssp.pb
#ZZADDDDDFBDFKLCCDDDDDEHIADDDFKLCFKLMBACFBEGHJACCDCDDDDFKPAFBDCDFKLGCCDDDEEHHIACDDDDDFBFKLCCFBFBFBACFBFKBCCDFKLOCGHJOPMLMZZ
#>2000.pdb_A.dssp.pb
#ZZDCDDDDFBDFKLCCDDDDDEHIADDDFBDCFKBMBACDDEEHJBMCDCDDDDFKPAFBACDFKBGCDDDDEDEHIACDDDDDFBFKBCCFBFBFBDCFBFKBCCDFKLPCGHJOJMLMZZ
#>3000.pdb_A.dssp.pb
#ZZDCDDDDFBDFKLCCDDDDDEHIADDDFBDCFKBMCDCDDEEEHIACDCDDDDFKPAFBDCDFKBCCDDDDEHEHIACDDDDDFBFKLCCFBFBLBACFBFKBCCDFKLOCGHJOJKLMZZ
#>4000.pdb_A.dssp.pb
#ZZACDDDDFBDFKLCCDDDDDDDDEHIAFBDCFKBFBDCDDEEHHIABDCDDDDDDCCFBACDFKBGCDDDDDDEHIACDDDDDFBFKBCCFBFBFBEHKAFKBCCDFKLGCGHJAJKLMZZ
#>5000.pdb_A.dssp.pb
#ZZDCDDDDDDDFKBCCDDDDDEHIADDDFBDCFKBFBGCIAEEEHIABEHIADDDDDDFBACDFKBCCDDDDEDEHIACDDDDDFBFKLCCFBFBLBGCKBFKBCCDFKLGCGHJLJMLMZZ
#>6000.pdb_A.dssp.pb
#ZZACDDDDDDDFKBCCDDDDDEHIADDDFBDCFKBFBGHIAEEJHIACEHIADDFBDCFBDCDFKBCCDDDDFBGHIACDDDDDFBFKLCCFBFBFBGCFBFKBCCDFKLGCGHJLJMLMZZ

#IT ALSO REQUIRES A CSV FILE CONTAINING THE SCORES FOR PB TRANSITIONS DERIVED FROM PB SUBSTITUTION MATRIX
#THIS CSV FILE IS PROVIDED: matrix.csv

#THE COMMAND LINE TO EXECUTE THE SCRIPT IS AS FOLLOWS:
# nohup score.sh myfile.pb &
# EXECUTION TAKES TIME (several hours for 200 snapshots on 1 CPU)

# IT GENERATES ONE OUPUT FILE CALLED kimura.txt
# THIS FILE IS CSV FORMATTED:
#  0;  0;122;118;  4;0.967;0.033;238.09;0.034
#  0;  1;122; 86; 36;0.705;0.295;123.33;0.375
#  0;  2;122; 91; 31;0.746;0.254;129.14;0.311
#  0;  3;122; 88; 34;0.721;0.279;130.73;0.348
#  0;  4;122; 96; 26;0.787;0.213;154.07;0.251
#  0;  5;122; 90; 32;0.738;0.262;139.21;0.323
#  0;  6;122; 93; 29;0.762;0.238;142.75;0.286
# EACH LINE CONTAINS A PAIRWISE COMPARISON
# THE 9 COLUMNS ARE AS FOLLOWS:
# ID OF PB SEQUENCE 1; 
# ID OF PB SEQUENCE 2;
# LENGTH OF PB SEQUENCE = NUMBER OF POSITIONS (npos);
# COUNT OF POSITIONS WITH POSITIVE PB TRANSITIONS;
# COUNT OF POSITIONS WITH NEGATIVE PB TRANSITIONS;
# FRACTION OF POSITIVE TRANSITIONS;
# FRACTION OF NEGATIVE TRANSITIONS;
# SUM OF SCORES;
# KIMURA DISTANCE

# IT ALSO GENERATE TEXT FILE mymd.txt CONTAINING ONE PB SEQUENCES PER LINE WITHOUT HEADERS
# TEMPORARY FILES ARE ALSO GENERATED, ONE PER SNAPSHOT, WITH .tmp EXTENSION. THESE ARE REMOVED WHEN THE SCRIPT TERMINATES

N=2
grep -v "^>" $1 > mymd.txt #extract from fasta sequences, only PB sequences
mycount=`cat mymd.txt | wc -l` #count the number of PB sequences in MD file
#echo "$mycount"
mymd=`cat mymd.txt` #stores in a variable all the PB sequences
myarr=($mymd) #stores in an array all PB sequences from MD

rm kimura.txt titi.txt scores.txt

#-------------------------
calculate_pair() #function to score a pair of PB sequences using the PB substitution matrix (requires input file matrix.csv)
#it calculates the Kimura protein dissimilarity distance
{
	#LOCAL VARIABLES DEFINITION
	local mypairs=("${@}") #get all PB pairs between two PB sequences
	local mylength=${#mypairs[@]} #gets number of positions - should be length of PB sequence
	local mylength=$((mylength-2)) #To account for the last two ZZ
	local P=0 # this will store the count of positions with positive PB substitution scores
	local Q=0 # this will store the count of positions with negative PB substitution scores
	local myscore=0 # this will store the raw score i.e the sum of PB substitution scores between two PB sequences
	#local mystart=$2
	#local myend=$3

	#EXTRACTING PB SUBSTITUTIONS SCORES
	for k in `seq 2 $((mylength-1))`; #for each position up to length of PB sequence -1, starting from 2 to exclude the first ZZ
	#for k in `seq $mystart $((myend-1))`
	do
		local score=$(grep "${mypairs[$k]}" matrix.csv | awk -F ";" '{print $2}') #grep PB pairs in MD file and get the corresponding scores from PB substitution matrix
		#printf "%5.2f;" $score >> toto.txt
		#echo " Pair: score ${mypairs[$k]} $score" 
		if (( $(echo "$score" | awk '{print ($1 > 0)}') )); then
			P=$((P+1)) #counts the number of positive PB substitutions
		elif (( $(echo "$score" | awk '{print ($1 < 0)}') )); then
				Q=$((Q+1)) #counts the number of negative PB substitutions
		fi
		myscore=`echo $myscore $score | awk '{print $1 + $2}'` #sum the pb substitution scores
	done

	#CALCULATION OF KIMURA PROTEIN DISSIMILARITY DISTANCE 
	local myP=`echo $P $((mylength-2)) | awk '{print $1 / $2}'` #calculates the fraction of positions with positive scores
	local myQ=`echo $Q $((mylength-2)) | awk '{print $1 / $2}'` #calculates the fraction of positions with negative scores
	local mykimura=`awk -v D=$myQ 'BEGIN{print (-1*log(1 - D - 0.2 * D * D))}'` #calculates Kimura protein distance

	#PRINTS THE RESULT TO A TEXT FILE
	printf "%3d;%3d;%3d;%5.3f;%5.3f;%5.2f;%5.3f" $((mylength-2)) $P $Q $myP $myQ $myscore $mykimura #>> kimura.txt

} # end of function
#-------------------------

#for i in `seq 0 $((mycount-1))`
for i in `seq 0 $((mycount-1))`
do
	#rm $i.txt
	#create a temporary file with .tmp extension for each snapshot where the PB sequence is written in a column
	grep -o . <<< ${myarr[$i]} > $i.tmp
done
 
#for i in `seq 0 $((mycount-2))` # this loops over all snapshots
for i in `seq 0 $((mycount-1))` # this loops over all snapshots
do 
	for j in `seq $(($i)) $((mycount-1))` #this again loops over all snapshots - to create all possible pairs of PB sequences
	do
		if [[ $i -le $j ]]; then
			mypbpairs=($(paste -d , $i.tmp $j.tmp | tr -d ',')) #get list of PBs pairs arranged in an array
			#printf "%3d;%3d;" $i $j >> kimura.txt 
			myresults=`calculate_pair "${mypbpairs[@]}"` #score the dissimilarity between a pair of PB sequence
			printf "%3d;%3d;$myresults\n" $i $j >> kimura.txt
		fi			
	done
done

rm *.tmp #removes temporary tmp .files




