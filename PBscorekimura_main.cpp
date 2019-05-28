//THIS C++ PROGRAM CALCULATES KIMURA DISTANCE BETWEEN ALL PAIRWISE DISTANCES BETWEEN PB SEQUENCES EXTRACTED FROM MD SIMULATION
// CONTACT DETAILS: Lionel.hoffmann@univ-nantes.fr

//TO EXECUTE THIS SCRIPT, A TEXT FILE CONTAINTING PB SEQUENCES IN FASTA FORMAT IS REQUIRED
// EXAMPLE OF CONTENT OF PB SEQUENCE FILE :
//>0.pdb_A.dssp.pb
//ZZACDDDDFBDFKBCCDDDDDEDDEHIAFBDCFKBCCDFBDGHKLPCCFBDCDDDDDDFBACDFKLCCDDDFBGNHIACDDDDDFBFKLCCFBFBFBDCFBFKBCCDFKLGMGHPLMMMMZZ
//>1000.pdb_A.dssp.pb
//ZZADDDDDFBDFKLCCDDDDDEHIADDDFKLCFKLMBACFBEGHJACCDCDDDDFKPAFBDCDFKLGCCDDDEEHHIACDDDDDFBFKLCCFBFBFBACFBFKBCCDFKLOCGHJOPMLMZZ
//>2000.pdb_A.dssp.pb
//ZZDCDDDDFBDFKLCCDDDDDEHIADDDFBDCFKBMBACDDEEHJBMCDCDDDDFKPAFBACDFKBGCDDDDEDEHIACDDDDDFBFKBCCFBFBFBDCFBFKBCCDFKLPCGHJOJMLMZZ
//>3000.pdb_A.dssp.pb
//ZZDCDDDDFBDFKLCCDDDDDEHIADDDFBDCFKBMCDCDDEEEHIACDCDDDDFKPAFBDCDFKBCCDDDDEHEHIACDDDDDFBFKLCCFBFBLBACFBFKBCCDFKLOCGHJOJKLMZZ
//>4000.pdb_A.dssp.pb
//ZZACDDDDFBDFKLCCDDDDDDDDEHIAFBDCFKBFBDCDDEEHHIABDCDDDDDDCCFBACDFKBGCDDDDDDEHIACDDDDDFBFKBCCFBFBFBEHKAFKBCCDFKLGCGHJAJKLMZZ
//>5000.pdb_A.dssp.pb
//ZZDCDDDDDDDFKBCCDDDDDEHIADDDFBDCFKBFBGCIAEEEHIABEHIADDDDDDFBACDFKBCCDDDDEDEHIACDDDDDFBFKLCCFBFBLBGCKBFKBCCDFKLGCGHJLJMLMZZ
//>6000.pdb_A.dssp.pb
//ZZACDDDDDDDFKBCCDDDDDEHIADDDFBDCFKBFBGHIAEEJHIACEHIADDFBDCFBDCDFKBCCDDDDFBGHIACDDDDDFBFKLCCFBFBFBGCFBFKBCCDFKLGCGHJLJMLMZZ

//IT ALSO REQUIRES A CSV FILE CONTAINING THE SCORES FOR PB TRANSITIONS DERIVED FROM PB SUBSTITUTION MATRIX
//THIS CSV FILE IS PROVIDED: matrix.csv

//THE COMMAND LINE TO EXECUTE THE SCRIPT IS AS FOLLOWS:
// nohup score.sh myfile.pb &
// EXECUTION TAKES TIME (several hours for 200 snapshots on 1 CPU)

// IT GENERATES ONE OUPUT FILE CALLED kimura.txt
// THIS FILE IS CSV FORMATTED:
//  0;  0;122;118;  4;0.967;0.033;238.09;0.034
//  0;  1;122; 86; 36;0.705;0.295;123.33;0.375
//  0;  2;122; 91; 31;0.746;0.254;129.14;0.311
//  0;  3;122; 88; 34;0.721;0.279;130.73;0.348
//  0;  4;122; 96; 26;0.787;0.213;154.07;0.251
//  0;  5;122; 90; 32;0.738;0.262;139.21;0.323
//  0;  6;122; 93; 29;0.762;0.238;142.75;0.286
// EACH LINE CONTAINS A PAIRWISE COMPARISON
// THE 9 COLUMNS ARE AS FOLLOWS:
// ID OF PB SEQUENCE 1; 
// ID OF PB SEQUENCE 2;
// LENGTH OF PB SEQUENCE = NUMBER OF POSITIONS (npos);
// COUNT OF POSITIONS WITH POSITIVE PB TRANSITIONS;
// COUNT OF POSITIONS WITH NEGATIVE PB TRANSITIONS;
// FRACTION OF POSITIVE TRANSITIONS;
// FRACTION OF NEGATIVE TRANSITIONS;
// SUM OF SCORES;
// KIMURA DISTANCE

#include "PBscoreClasses.hpp"
#include <iostream>
#include <iomanip>      // std::setprecision
#include <cstdlib>			// std::system
#include <algorithm>    // std::count
#include <fstream>
#include <vector>
#include <map>
#include <math.h>

using namespace std ;

int main (int argc, char *argv[])
{
  // reading mode opening of the file containing the reference substitution matrix
  
  string fichier = "PBsubstitutionmatrix.csv" ;
	string commande = "sed -i \"s/;/ /g\" " + fichier ;		// changing file separator : switching from ";" to " "
	system(commande.c_str()) ;														
  ifstream flux_in (fichier.c_str(), ios::in) ;
  if (!flux_in) {
    cout << "oups !!" << endl ;
    return 0 ;
  }

	string PBtrans ; // this variable contains a single PB transition name, for example, A to B is AB
	double transscore ; // this variable contains a PB transition score 
	map <string, double> matrix ; // this dynamic indexed table contains the substitution matrix 
	while (flux_in >> PBtrans >> transscore) {
		matrix[PBtrans] = transscore ;						// reading the substitution matrix from the substitution matrix file 
	}
	flux_in.close() ;
	commande = "sed -i \"s/ /;/g\" " + fichier ;		// putting back the reference substitution matrix file 
	system(commande.c_str()) ;											// in its initial state  

  // reading mode opening of the file containing PB sequences  
  
  string fichier2 = string(argv[1]) ;
	commande = "grep -v \"^>\" " + fichier2 + "> PBsequences.txt" ; // creating new sequence file where sequence ID lines 
	system(commande.c_str()) ;                                      // have been removed  
	string fichier3 = "PBsequences.txt" ;
  ifstream flux_in2 (fichier3.c_str(), ios::in) ;
  if (!flux_in2) {
    cout << "zut !!" << endl ;
		return 0 ;
  }

  // reading PB sequences 
  
	Sequence PBseq ; // this variable contains a single PB sequence 
	string pb1, pb2 ; // pb1 is the initial state and pb2 is the final state of the PB transition 
	vector<string> seqlist ; // this dynamic table contains all the PB sequences 
	Dataline datalinetemp ; // buffer variable containing all the results for one PB sequence 
	vector<Dataline> result ; // this dynamic table contains all the results for all PB sequences 
	double sumscore = 0 ; // sum of all scores for one PB sequence 
	double numberofZ = 0 ; // number of Z in one PB sequence 
	double D ; // the D from the kimura distance formula : -ln(1 - D - 0.2D²) 
	while (flux_in2 >> PBseq.WholeSequence) {
		seqlist.push_back(PBseq.WholeSequence) ; // reading the PB sequences 
	}

	// treating PB sequences 

	string seq1 = seqlist[0] ;
	string seq2 ;
	for (unsigned int cpt = 0 ; cpt < seqlist.size() ; cpt++) { // the first loop is to load 2 consecutive sequences 
		seq2 = seqlist[cpt] ;																			// into appropriate variables 
		for (unsigned int k=0 ; k < seq2.length() ; k++) {				// the second loop is to calculate the different values 
			pb1 = seq1[k] ;																					// for each sequence 
			pb2 = seq2[k] ;
			PBtrans =  pb1 + pb2 ; 
			transscore = matrix[PBtrans] ;
			sumscore += transscore ;
			datalinetemp.sumscore = sumscore ;
			datalinetemp.pbid1 = 0 ;
			datalinetemp.pbid2 = cpt ;
			datalinetemp.pblength = seq2.length() ;
			if (transscore > 0) {++datalinetemp.postrans ; }
			if (transscore < 0) {++datalinetemp.negtrans ; }
			if (transscore == 0) {++numberofZ ; }
		}
		datalinetemp.posfrac = datalinetemp.postrans/(datalinetemp.pblength - numberofZ) ;
		datalinetemp.negfrac = datalinetemp.negtrans/(datalinetemp.pblength - numberofZ) ;
		D = 1 - datalinetemp.posfrac ;
		double X = 1 - D - (0.2 * D * D) ;											// those 2 lines are to calculate kimura distance 
		datalinetemp.kimura = -1*log(X) ;												// for each PB sequence 
		result.push_back(datalinetemp) ;
		sumscore = 0 ;
		datalinetemp.postrans = 0 ;
		datalinetemp.negtrans = 0 ;
		numberofZ = 0 ;
	}
	result.pop_back() ;

  // writing mode opening of a file to print the results 

  string resultfile = "kimura.txt" ;		                  
  ofstream flux_out (resultfile.c_str(), ios::out | ios::app) ; 

	// printing the results in "kimura.txt" file for each PB sequence 

	for (unsigned int cpt =0 ; cpt <= result.size() ; cpt++) {
		flux_out << result[cpt] ;
	}
	matrix.clear () ;
	seqlist.clear () ;
	result.clear () ;
	return 0 ;
}
