#include <string.h>
#include <iostream>
#include <math.h>
#include "PBscoreClasses.hpp"

using namespace std ;

// constructor function  of Sequence class 

Sequence::Sequence ()
{ 
  WholeSequence = "ABC" ;
}

// constructor function of Dataline class 

Dataline::Dataline ()
{ 
	pbid1 = 0 ;
	pbid2 = 0 ;
	pblength = 0 ;
  postrans = 0 ;
	negtrans = 0 ;
	posfrac = 0 ;
	negfrac = 0 ;
	sumscore = 0 ;
	kimura = 0 ;
}

// overload of <<, >> and = operator of Dataline class 

ostream& operator<< (ostream& flux_out, Dataline& dataline)
{
  flux_out << dataline.pbid1 << ";" << dataline.pbid2 << ";" << dataline.pblength << ";" << dataline.postrans << ";" << dataline.negtrans << ";" << dataline.posfrac << ";" << dataline.negfrac << ";" << dataline.sumscore << ";" << dataline.kimura << endl ;
  return flux_out ;
}

istream& operator>> (istream& flux_in, Dataline& dataline)
{
  flux_in >> dataline.pbid1 >> dataline.pbid2 >> dataline.pblength >> dataline.postrans >> dataline.negtrans >> dataline.posfrac >> dataline.negfrac >> dataline.sumscore >> dataline.kimura ;
  return flux_in ;
}
