//THIS C++ PROGRAM CALCULATES KIMURA DISTANCE BETWEEN ALL PAIRWISE DISTANCES BETWEEN PB SEQUENCES EXTRACTED FROM MD SIMULATION
// CONTACT DETAILS: Lionel.hoffmann@univ-nantes.fr

//TO EXECUTE THIS SCRIPT, A TEXT FILE CONTAINTING PB SEQUENCES IN FASTA FORMAT IS REQUIRED
// EXAMPLE OF CONTENT OF 5 PB SEQUENCES FILE :
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

//IT ALSO REQUIRES A CSV FILE CONTAINING THE SCORES FOR PB TRANSITIONS DERIVED FROM PB SUBSTITUTION MATRIX
//THIS CSV FILE IS PROVIDED: PB substitutionmatrix.csv
// THE COMMAND LINE TO EXECUTE THE SCRIPT IS AS FOLLOWS:
// nohup ./PBscorekimura mysequencefile.pb &
// EXECUTION IS QUICK (a few seconds for 200 snapshots on 1 CPU)

// IT GENERATES ONE OUPUT FILE CALLED kimura.txt
// THIS FILE IS CSV FORMATTED:
// 0;1;122;86;32;0.728814;0.271186;123.33;0.336725
// 0;2;122;91;27;0.771186;0.228814;129.14;0.273496
// 0;3;122;88;30;0.745763;0.254237;130.73;0.310834
// 0;4;122;96;22;0.813559;0.186441;154.07;0.214918
// 1;2;122;108;10;0.915254;0.0847458;200.47;0.090124
// 1;3;122;109;9;0.923729;0.0762712;196.84;0.0805971
// 1;4;122;91;27;0.771186;0.228814;142.31;0.273496
// 2;3;122;110;8;0.932203;0.0677966;207.87;0.0711909
// 2;4;122;96;22;0.813559;0.186441;168.96;0.214918
// 3;4;122;94;24;0.79661;0.20339;163.84;0.23783
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

	string PBtrans ; // this variable contains a single PB transition name ; for example, A to B is AB
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
	vector<string> seqlist ; // this dynamic table contains all the PB sequences

	while (flux_in2 >> PBseq.WholeSequence) {
		seqlist.push_back(PBseq.WholeSequence) ; // reading the PB sequences 
	}

  // writing mode opening of a file to print the results 

  string resultfile = "kimura.txt" ;		                  
  ofstream flux_out (resultfile.c_str(), ios::out | ios::app) ; 

	// treating PB sequences
 
	string pb1, pb2 ; // pb1 is the initial state and pb2 is the final state of the PB transition  
	Dataline datalinetemp ; // buffer variable containing all the results for one PB sequence 
	vector<Dataline> result ; // this dynamic table contains all the results for all PB sequences 
	double sumscore = 0 ; // sum of all scores for one PB sequence 
	double numberofZ = 0 ; // number of Z in one PB sequence 
	double D ; // the D from the kimura distance formula : -ln(1 - D - 0.2D²)
	string seq1 ; // first sequence to compare
	string seq2 ; // second sequence to compare 

	for (unsigned int i = 0 ; i < seqlist.size() ; i++) {					// the first loop is to select the first sequence to compare
		seq1 = seqlist[i] ;																					// and load it into an appropriate variable 
		for (unsigned int j = i + 1 ; j < seqlist.size() ; j++) { 			// the second loop is to load the second sequence to compare 
			seq2 = seqlist[j] ;																				// and load it into an appropriate variable 
			for (unsigned int k=0 ; k < seq1.length() ; k++) {				// the third loop is to calculate the different values 
				pb1 = seq1[k] ;																					// for each sequence and for each pair 
				pb2 = seq2[k] ;
				PBtrans =  pb1 + pb2 ; 
				transscore = matrix[PBtrans] ;
				sumscore += transscore ;
				datalinetemp.sumscore = sumscore ;
				datalinetemp.pbid1 = i ;
				datalinetemp.pbid2 = j ;
				datalinetemp.pblength = seq1.length() ;
				if ((pb1 == "-") || (pb2 == "-")) {++datalinetemp.negtrans ; }
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
	}
	result.pop_back() ; 																				// removing extra line at the end of the table 

	// printing the results in "kimura.txt" file for each PB sequence 

	for (unsigned int cpt =0 ; cpt <= result.size() ; cpt++) {
		flux_out << result[cpt] ;
	}
	matrix.clear() ;
	seqlist.clear() ;
	result.clear() ;
	flux_in2.close() ;
	flux_out.close() ;
	commande = "rm -f PBsequences.txt" ;												// removing the file created at the beginning 
	system(commande.c_str()) ;																	// of the program which contains only PB sequences 
	return 0 ;
}
