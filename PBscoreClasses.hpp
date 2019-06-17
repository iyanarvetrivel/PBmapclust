#include <string>

using namespace std ;


class Sequence
{
  // Attribute
  public:

	string WholeSequence ;

  // constructor 
  Sequence () ;

};

class Dataline 
{
  // Attributes
  public:

	int pbid1, pbid2 ; // ID OF PB SEQUENCES
	double pblength, postrans, negtrans ; 
// LENGTH OF PB SEQUENCE = NUMBER OF POSITIONS (npos);
// COUNT OF POSITIONS WITH POSITIVE PB TRANSITIONS;
// COUNT OF POSITIONS WITH NEGATIVE PB TRANSITIONS;
	double posfrac, negfrac, sumscore, kimura ;
// FRACTION OF POSITIVE TRANSITIONS;
// FRACTION OF NEGATIVE TRANSITIONS;
// SUM OF SCORES;
// KIMURA DISTANCE

  // constructor 
  Dataline () ;

  // overload of <<, >> and = operator 
  friend ostream& operator<< (ostream&, Dataline&) ;
  friend istream& operator>> (istream&, Dataline&) ;
};
