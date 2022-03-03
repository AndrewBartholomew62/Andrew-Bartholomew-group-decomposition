/***********************************************************************
			   bigreal
	       A. Bartholomew 4th May 2008

    bigreal is an arbitrary precision real number class intended to support
	a class of linear combinations of surds.
	
	A bigreal is comprised of an integral part and a decimal part, each of which 
	are recorded as bigints.

    The sign is maintained separately as a boolean.


***********************************************************************/

#include <iostream>

#include <bigint.h>

class bigreal
{
    bigint integral;
    bigint decimal;
	bool negative;
	int places;      // number of decimal places in decimal
	static int div_num_places;  // number of additional digits of precision for division operator
	
	
public:
	
    bigreal():integral(0),decimal(0),negative(false),places(0){}
	bigreal(const int num): integral(abs(num)),decimal(0),negative((num<0?true:false)),places(0){}
	
	bigint& integral_part() {return integral;}
	bigint& decimal_part() {return decimal;}
	const int decimal_places() const {return places;}
	void set_decimal_places(int num_places) {places = num_places;}
	bigint get_div_num_places() {return div_num_places;}
	void set_div_num_places(int num_places) {div_num_places = num_places;}
	
//	operator int(); // conversion operator to an int
//	operator bool(); // conversion operator to a bool

    bigreal& operator ++ (); //prefix
    bigreal operator ++ (int); //postfix
    bigreal& operator -- ();
    bigreal operator -- (int);

	void dump(ostream& s) const;
//	void sanitize();

	friend void sum(bigreal& result,const bigreal& a, const bigreal& b);
	friend void diff(bigreal& result,const bigreal& a, const bigreal& b);

    friend ostream& operator << (ostream& os, const bigreal& a);
    friend istream& operator >> (istream& is, bigreal& a);


    friend bigreal operator + (const bigreal& a, const bigreal& b);
    friend bigreal operator - (const bigreal& a, const bigreal& b);
    friend bigreal operator * (const bigreal& a, const bigreal& b);
    friend bigreal operator / (const bigreal& a, const bigreal& b);
//    friend bigreal operator % (const bigreal& a, const bigreal& b);
    friend bigreal operator += (bigreal& a, const bigreal& b);
    friend bigreal operator -= (bigreal& a, const bigreal& b);
    friend bigreal operator *= (bigreal& a, const bigreal& b);
    friend bigreal operator /= (bigreal& a, const bigreal& b);
//    friend bigreal operator %= (bigreal& a, const bigreal& b);
//	friend bigreal operator >>= (bigreal& a, int n);
//	friend bigreal operator <<= (bigreal& a, int n);
    friend bool operator == (const bigreal& a, const bigreal& b);
    friend bool operator != (const bigreal& a, const bigreal& b);
    friend bool operator > (const bigreal& a, const bigreal& b);
    friend bool operator < (const bigreal& a, const bigreal& b);
    friend bool operator >= (const bigreal& a, const bigreal& b);
    friend bool operator <= (const bigreal& a, const bigreal& b);
    friend bigreal abs (const bigreal& a);
    friend int num_len(const bigreal& a);
//	friend bigreal gcd (const bigreal& u, const bigreal& v);

};

#include <bigreal_control.h>

struct bigreal_error {
	bigreal_error (string message) {cout << "\nbigreal error!" << message << endl;}
};

void carry(bigint& src, bigint& dst, int n);
