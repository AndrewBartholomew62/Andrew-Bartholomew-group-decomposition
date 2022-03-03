/***********************************************************************
			   surd.h

Original coding as surd.h A. Bartholomew 9th May 2008

Split into surd.spp and surd.h Mat 2010


***********************************************************************/

#include <list>
#include <sstream>
#include <rational.h>
#include <bigreal.h>

bigreal operator * (const bigreal& a, rational<int> r);
bigreal operator += (bigreal& a, rational<int> r);




#include <surd_control.h>

struct Surd
{
	int prime;
	bigreal approx;
	rational<int> coeff;
	
	Surd(): prime(1),approx(0),coeff(0){}
	Surd (const int p);

	Surd operator *= (const rational<int> r) { coeff *= r; return *this;}
	void dump(ostream&) const;
	
    friend istream& operator >> (istream& is, Surd& S);
};

class surd
{

	list<Surd> terms;	
	
public:
	
    surd(){}
	void dump(ostream&) const;
	bool zero() const {return !terms.size();}
	bigreal estimate();
	
    friend ostream& operator << (ostream& os, const surd& s);
    friend surd operator + (const surd& a, const surd& b);
    friend surd operator + (const surd& s, const Surd& t);
    friend surd operator += (const surd& s, const Surd& t);
    friend surd operator += (surd& s, const Surd& t);
	friend surd operator - (const surd& a, const surd& b);
    friend surd operator - (const surd& s, const Surd& t);
    friend surd operator -= (const surd& s, const Surd& t);
    friend surd operator -= (surd& s, const Surd& t);
    friend bool operator == (const surd& a, const surd& b);
    friend bool operator != (const surd& a, const surd& b);
    friend bool operator > (const surd& a, const surd& b);
    friend istream& operator >> (istream& is, surd& s);

};

surd operator + (const surd& a, const surd& b);
surd operator - (const surd& a, const surd& b);
surd operator += (surd& a, const surd& b);
surd operator -= (surd& a, const surd& b);
bool operator == (const surd& a, const surd& b);
bool operator != (const surd& a, const surd& b);
bool operator > (const surd& aa, const surd& b);
bool operator < (const surd& a, const surd& b);
bool operator >= (const surd& a, const surd& b);
bool operator <= (const surd& a, const surd& b);



