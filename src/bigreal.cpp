/*********************************************************************** 

	bigreal A. Bartholomew 4th May 2008


***********************************************************************/
using namespace std;
#include <fstream>
#include <iostream>
#include <cstdlib>
//#include <vector>
//#include <string>
//#include <iomanip>
#include <sstream>

extern ofstream debug;

//#include <util.h>
#include <bigreal.h>

/*********** Function prototypes **************/
void carry(bigint& src, bigint& dst, int n);
ostream& operator << (ostream& os, const bigreal& a);
istream& operator >> (istream& is, bigreal& a);

bigreal& bigreal::operator ++ () //prefix
{
    bigreal local = 1;
    *this += local;
    return *this;
}

bigreal bigreal::operator ++ (int) //postfix
{
    bigreal tmp = *this;
    ++*this;
    return tmp;
}

bigreal& bigreal::operator -- () //prefix
{
    bigreal local = 1;
    *this -= local;
    return *this;
}

bigreal bigreal::operator -- (int) //postfix
{
    bigreal tmp = *this;
    --*this;
    return tmp;
}

void bigreal::dump (ostream& os) const
{
	os << "bigreal(int=" << integral << ",dec=" << decimal << ", negative=" << (negative?"true":"false") 
	   << ",places=" << places << ")" << flush;
}

/* sum calculates the sum of the bigreals at aa and bb and stores it at
   result.  The signs of aa and bb must have been taken care of elsewhere, 
   sum asumes the negative flag in result is correct.
*/
void sum(bigreal& result,const bigreal& aa, const bigreal& bb)
{

if (bigreal_control::DEBUG & bigreal_control::sum)
{
	debug << "bigreal::sum : a = " << aa << 
	         "\n               b = " << bb << endl;
	debug << "bigreal::sum : a = ";
	aa.dump(debug);
	debug << endl;
	debug << "bigreal::sum : b = ";
	bb.dump(debug);
	debug << endl;
}

    /* check for a == 0 */
    if (aa == bigreal(0))
    {
		result = bb;
		return;
    }

    if (bb == bigreal(0))
	{
		result = aa;
		return;
	}

	bigreal a = aa;
	bigreal b = bb;

	result.integral = a.integral + b.integral;
	result.places = a.places;

	if (a.places > b.places)
	{
		for (int i=b.places; i < a.places; i++)
			b.decimal *= 10; 
	}
	else if (b.places > a.places)
	{
		for (int i=a.places; i < b.places; i++)
			a.decimal *= 10; 
		
		result.places = b.places;
	}

if (bigreal_control::DEBUG & bigreal_control::sum)
{
	debug << "bigreal::sum : after adjusting, a.decimal = " << a.decimal << ", b.decimal = " << b.decimal << endl;
	debug << "bigreal::sum : result = ";
	result.dump(debug);
	debug << endl;
}
	result.decimal = a.decimal + b.decimal;
	carry (result.decimal, result.integral, result.places);

if (bigreal_control::DEBUG & bigreal_control::sum)
{
	debug << "bigreal::sum : result = " << result << endl;
	debug << "bigreal::sum : result = ";
	result.dump(debug);
	debug << endl;
}
	
}

/* diff calculates the difference between the bigreals at a and b and stores it at
   result.  The signs must have been taken care of elsewhere and
   a must be greater than or equal to b.

   diff clears the negative flag in result if the difference is 0
*/
void diff(bigreal& result, const bigreal& aa, const bigreal& bb)
{

if (bigreal_control::DEBUG & bigreal_control::diff)
{
	debug << "bigreal::diff : a = " << aa << 
	         "\n                b = " << bb << endl;
	debug << "bigreal::diff : a = ";
	aa.dump(debug);
	debug << endl;
	debug << "bigreal::diff : b = ";
	bb.dump(debug);
	debug << endl;
}

	bigreal a = aa;
	bigreal b = bb;

	result.places = a.places; // has to be here in case places of a and b are equal

	if (a.places > b.places)
	{
		for (int i=b.places; i < a.places; i++)
			b.decimal *= 10; 
	}
	else if (b.places > a.places)
	{
		for (int i=a.places; i < b.places; i++)
			a.decimal *= 10; 
		
		result.places = b.places;
	}

	if (b.decimal > a.decimal)
	{
		/* need to carry from a.integral */
		bigint borrow(1);
		for (int i=0; i< result.places; i++)
			borrow*= 10;
		a.decimal += borrow;	

if (bigreal_control::DEBUG & bigreal_control::diff)
	debug << "bigreal::diff : b.decimal > a.decimal, after borrow a.decimal = " << a.decimal << endl;

		a.integral--;
		
	}

	result.decimal = a.decimal-b.decimal;	
	result.integral = a.integral-b.integral;
	
	if (result.integral < bigint(0))
		throw bigreal_error("diff attempting to subtract b from a where b>a");

    /* check for the zero result */
    bool loc_sign = result.negative;
    result.negative = false;
    if (!(result == bigreal(0)))
		result.negative = loc_sign;

if (bigreal_control::DEBUG & bigreal_control::diff)
{
	debug << "bigreal::diff : result = " << result << endl;
	debug << "bigreal::diff : result = ";
	result.dump(debug);
	debug << endl;
}

//	result.sanitize();

}

ostream& operator << (ostream& os, const bigreal& a)
{
	
	if (a.negative)
		os << "-";
		
	os << a.integral;
	if (a.decimal != bigint(0))
	{
		os << '.';
		int len = num_len(a.decimal);
		if (len < a.places)
		{
			for (int i=len; i< a.places; i++)
				os << "0";
		}
		os << a.decimal;
	}
	
	return os;
}

istream& operator >> (istream& is, bigreal& a)
{
	is >> a.integral;
	a.decimal = bigint(0);
	a.places = 0;

	char ch;
	is.get(ch);
	
	if (is.fail())
	{
		is.clear(); // no more characters to read, decimal is zero
	}
	else
	{
//cout << "\nis good = " << is.good();
//cout << "\nis eof = " << is.eof();
//cout << "\nis fail = " << is.fail();
//cout << "\nis bad = " << is.bad();
//cout << endl;

		if (ch == '.')
		{
			int leading_zeros = 0;
			do
			{
				is.get(ch);
				
				if (ch == '0')
					leading_zeros++;
					
			} while (ch == '0');
			is.putback(ch);
			is >> a.decimal;
			a.places = num_len(a.decimal)+leading_zeros;		
		}
		else
		{
			is.putback(ch);
		}
	}
	
	/* ajust the sign */
	if (a.integral < bigint(0))
	{
		a.negative = true;
		a.integral *= -1;
	}
	else
	{
		a.negative = false;
	}
	
	return is;
}

bigreal operator + (const bigreal& a, const bigreal& b)
{

if (bigreal_control::DEBUG & bigreal_control::add)
{
	debug << "bigreal::operator + : a = " << a << 
	         "\n                      b = " << b << endl;
}
    bigreal result;

    if (a.negative && b.negative)
    {
if (bigreal_control::DEBUG & bigreal_control::add)
	debug << "bigreal::operator + : both negative" << endl;
		result.negative = true;
		sum(result,a,b);
    }
    else if (a.negative)
    {
		if ( abs(a) >= abs(b) )
		{
if (bigreal_control::DEBUG & bigreal_control::add)
	debug << "bigreal::operator + : a negative b positive result negative" << endl;
	    	result.negative = true;
	    	diff (result,a,b);
		}
		else
		{
if (bigreal_control::DEBUG & bigreal_control::add)
	debug << "bigreal::operator + : a negative b positive result positive" << endl;
		    result.negative = false;
		    diff (result,b,a);
		}
    }
    else if (b.negative)
    {
		if ( abs(a) >= abs(b) )
		{
if (bigreal_control::DEBUG & bigreal_control::add)
	debug << "bigreal::operator + : a positive b negative result positive" << endl;
	    	result.negative = false;
	    	diff (result,a,b);
		}
		else
		{
if (bigreal_control::DEBUG & bigreal_control::add)
	debug << "bigreal::operator + : a positive b negative result negative" << endl;
	    	result.negative = true;
	    	diff (result,b,a);
		}
    }
    else
    {
if (bigreal_control::DEBUG & bigreal_control::add)
	debug << "bigreal::operator + : both positive" << endl;
		result.negative = false;
		sum (result,a,b);
    }
    return result;
}

bigreal operator - (const bigreal& a, const bigreal& b)
{

if (bigreal_control::DEBUG & bigreal_control::subtract)
{
	debug << "bigreal::operator - : a = " << a << 
	         "\n                      b = " << b << endl;
}
	
	bigreal minusb = b*bigreal(-1);
	return a+minusb;
}

bigreal operator * (const bigreal& a, const bigreal& b)
{
    bigreal result;
	
if (bigreal_control::DEBUG & bigreal_control::multiply)
{
	debug << "bigreal::operator * : a = " << a << 
             "\n                      b = " << b << endl;
	debug << "bigreal::operator * : a = ";
	a.dump(debug);
	debug << endl;
	debug << "bigreal::operator * : b = ";
	b.dump(debug);
	debug << endl;
}

    /* first check for zero product */
    if (a == bigreal(0) || b == bigreal(0))
    {

if (bigreal_control::DEBUG & bigreal_control::multiply)
	debug << "bigreal::operator * : zero operand, returning zero" << endl;	

		return bigreal(0);
    }

    if (a == bigreal(1) )
    {

if (bigreal_control::DEBUG & bigreal_control::multiply)
	debug << "bigreal::operator * : a = 1, returning b" << endl;	

		return b;
    }

    if (b == bigreal(1) )
    {

if (bigreal_control::DEBUG & bigreal_control::multiply)
	debug << "bigreal::operator * : b = 1, returning a" << endl;	

		return a;
    }

    if ((a.negative || b.negative) && !(a.negative && b.negative))
        result.negative = true;
    else
        result.negative = false;
        
    /* now multiply b through by a */
	result.integral = a.integral * b.integral;

if (bigreal_control::DEBUG & bigreal_control::multiply)
	debug << "bigreal::operator * : integral product = " << result.integral << endl;
	
	result.decimal = a.integral * b.decimal;

if (bigreal_control::DEBUG & bigreal_control::multiply)
	debug << "bigreal::operator * : a.integral * b.decimal = " << result.decimal << endl;

	// add zero to result.decimal to produce a.places+b.places digits
	for (int i=0; i < a.places; i++)
		result.decimal *=10;
	
if (bigreal_control::DEBUG & bigreal_control::multiply)
	debug << "bigreal::operator * : result.decimal = " << result.decimal << endl;

	bigint temp = a.decimal * b.integral;
	
if (bigreal_control::DEBUG & bigreal_control::multiply)
	debug << "bigreal::operator * : a.decimal * b.integral = " << temp << endl;

	// add zero to temp to produce places+b.places digits
	for (int i=0; i < b.places; i++)
		temp *=10;

	result.decimal += temp;

if (bigreal_control::DEBUG & bigreal_control::multiply)
	debug << "bigreal::operator * : result.decimal = " << result.decimal << endl;
	
	temp = a.decimal * b.decimal;

if (bigreal_control::DEBUG & bigreal_control::multiply)
	debug << "bigreal::operator * : a.decimal * b.decimal = " << temp << endl;

	result.decimal += temp;
	
if (bigreal_control::DEBUG & bigreal_control::multiply)
	debug << "bigreal::operator * : before carry result.integral = " << result.integral << ", result.decimal = " << result.decimal << endl;

	carry (result.decimal, result.integral, a.places+b.places);

if (bigreal_control::DEBUG & bigreal_control::multiply)
	debug << "bigreal::operator * : after carry result.integral = " << result.integral << ", result.decimal = " << result.decimal << endl;

	result.places = a.places+b.places;


if (bigreal_control::DEBUG & bigreal_control::multiply)
{
	debug << "bigreal::operator * : result = " << result << endl;
	debug << "bigreal::operator * : result = ";
	result.dump(debug);
	debug << endl;
}
   return result;
}

bigreal operator / (const bigreal& a, const bigreal& b)
{
	bigreal u = a;
	bigreal v = b;

if (bigreal_control::DEBUG & bigreal_control::divide)
{
	debug << "bigreal::operator / : a = " << a << 
             "\n                      b = " << b << endl;
	debug << "bigreal::operator / : a = ";
	a.dump(debug);
	debug << endl;
	debug << "bigreal::operator / : b = ";
	b.dump(debug);
	debug << endl;
}

   	if (v == bigreal(0))
   	{
    	cout << "bigreal::operator / : zero divide error" << endl;
       	exit(0);
   	}

   	if (v == bigreal(1))
   	{

if (bigint_control::DEBUG & bigint_control::divide)
{
	debug << "bigreal::operator / : b = 1, returning a" << endl;
}
       	return u;
   	}

   	if (u == v)
   	{
if (bigint_control::DEBUG & bigint_control::divide)
{
	debug << "bigreal::operator / : u = v, returning 1" << endl;
}
    	bigreal result = bigreal(1);
       	return result;
   	}

	/* multiply up the integral parts and add the decimals */
	for (int i=0; i< u.places; i++)
		u.integral *= 10;
	
	u.integral += u.decimal;

	for (int i=0; i< v.places; i++)
		v.integral *= 10;
	
	v.integral += v.decimal;
	
if (bigreal_control::DEBUG & bigreal_control::divide)
{
	debug << "bigreal::operator / : after multiplying up " << endl; 
	debug << "bigreal::operator / :   a = " << u << 
             "\n                        b = " << v << endl;
	debug << "bigreal::operator / :   a = ";
	u.dump(debug);
	debug << endl;
	debug << "bigreal::operator / :   b = ";
	v.dump(debug);
	debug << endl;
}

	/* multiply u to get and additional div_num_places precision in the quotient */
	for (int i=0; i< bigreal::div_num_places; i++)
		u.integral *= 10;	
	
if (bigreal_control::DEBUG & bigreal_control::divide)
{
	debug << "bigreal::operator / : after adding div_num_places " << endl; 
	debug << "bigreal::operator / :   a = " << u << endl;
	debug << "bigreal::operator / :   a = ";
	u.dump(debug);
	debug << endl;
}

	/* now divide */
	bigint quotient = u.integral / v.integral;

if (bigreal_control::DEBUG & bigreal_control::divide)
	debug << "bigreal::operator / : quotient =  " << quotient << endl; 
	
	bigreal result;
	
 	if ((u.negative || v.negative) && !(u.negative && v.negative))
   		result.negative = true;
		
	result.places = u.places - v.places + bigreal::div_num_places;

if (bigreal_control::DEBUG & bigreal_control::divide)
	debug << "bigreal::operator / : result.places =  " << result.places << endl; 
	
	if (result.places < 0)
	{
		result.integral = quotient;
		for (; result.places < 0; result.places++)
			result.integral *= 10;

if (bigreal_control::DEBUG & bigreal_control::divide)
{
	debug << "bigreal::operator / : result.places negative, result.integral set to " 
	      << result.integral << endl; 
}

	}
	else
	{
		bigint factor = 1;
		for (int i=0; i< result.places; i++)
			factor *=10;
		
if (bigreal_control::DEBUG & bigreal_control::divide)
	debug << "bigreal::operator / : result.places positive, factor = " << factor << endl; 

		result.integral = quotient / factor;
		result.decimal = quotient % factor;
		
if (bigreal_control::DEBUG & bigreal_control::divide)
{
	debug << "bigreal::operator / : result before adjustment = ";
	result.dump(debug);
	debug << endl;
}

		/* adjust the number of decimal places to remove trailing zeros */
		if (result.decimal == bigint(0))
		{
			result.places = 0;
		}
		else
		{
			while (!(result.decimal % bigint(10)))
			{
				result.decimal /= bigint(10);
				result.places--;
			}
		}
	}
	
if (bigreal_control::DEBUG & bigreal_control::divide)
{
	debug << "bigreal::operator / : result = " << result << endl;
	debug << "bigreal::operator / : result = ";
	result.dump(debug);
	debug << endl;
}

	return result;
}


bigreal operator += (bigreal& a, const bigreal& b)
{
    return a = a+b;
}

bigreal operator -= (bigreal& a, const bigreal& b)
{
    return a = a-b;
}

bigreal operator *= (bigreal& a, const bigreal& b)
{
    return a = a*b;
}

bigreal operator /= (bigreal& a, const bigreal& b)
{
    return a = a/b;
}

bool operator == (const bigreal& a, const bigreal& b)
{

if (bigreal_control::DEBUG & bigreal_control::equal)
	debug << "bigreal::operator == : a = " << a << "\n                       b = " << b << endl;

    if (a.integral == b.integral && a.decimal == b.decimal && a.negative == b.negative)
		return true;
	else
		return false;
}

bool operator != (const bigreal& a, const bigreal& b)
{
	return (!(a == b));
}

bool operator > (const bigreal& a, const bigreal& b)
{

if (bigreal_control::DEBUG & bigreal_control::greater)
{
	debug << "bigreal::operator > : a = " << a << 
	         "\n                      b = " << b << endl;
	debug << "bigreal::operator > : a = ";
	a.dump(debug);
	debug << endl;
	debug << "bigreal::operator > : b = ";
	b.dump(debug);
	debug << endl;
}
	if  (a == b)
	{
if (bigreal_control::DEBUG & bigreal_control::greater)
	debug << "bigreal::operator > : a = b, return false" << endl; 
		return false;
	}
		
	if (a.negative && !b.negative)
	{
if (bigreal_control::DEBUG & bigreal_control::greater)
	debug << "bigreal::operator > : a negative b positive, return false" << endl; 
		return false;
	}
	
	if (!a.negative && b.negative)
	{
if (bigreal_control::DEBUG & bigreal_control::greater)
	debug << "bigreal::operator > : a positive b negative, return true" << endl; 
		return true;
	}
	
	/* take local copies */ //and sanitize to remove leading zeros 
	bigreal aa = a;
	bigreal bb = b;
//	aa.sanitize();
//	bb.sanitize();

	const bigreal* A;
	const bigreal* B;
	
	if (a.negative) // b.negative also, check |b| > |a|
	{
		A = &bb;
		B = &aa;
	}
	else // !a.negative && !b.negative, check |a| > |b|
	{
		A = &aa;
		B = &bb;
	}

if (bigreal_control::DEBUG & bigreal_control::greater)
{
	debug << "bigreal::operator > : checking A > B where A = " << *A << 
	       "\n                                           B = " << *B << endl; 
}
	
	/* so here we return true iff |A| > |B| */
	
	if (A->integral > B->integral)
	{
if (bigreal_control::DEBUG & bigreal_control::greater)
	debug << "bigreal::operator > : A->integral > B->integral, return true" << endl; 
		return true;
	}
	else if (A->integral < B->integral)
	{
if (bigreal_control::DEBUG & bigreal_control::greater)
	debug << "bigreal::operator > : A->integral < B->integral, return false" << endl; 
		return false;
	}

	/* even up the number of places in the decimal parts */
	if (aa.places < bb.places)
	{
		for (int i=aa.places; i<bb.places; i++)
			aa.decimal *= 10;
	}
	else if (aa.places > bb.places)
	{
		for (int i=b.places; i<aa.places; i++)
			bb.decimal *= 10;
	}
	
	if (A->decimal > B->decimal)
	{
if (bigreal_control::DEBUG & bigreal_control::greater)
	debug << "bigreal::operator > : A->decimal > B->decimal, return true" << endl; 
		return true;
	}
	else if (A->decimal < B->decimal)
	{
if (bigreal_control::DEBUG & bigreal_control::greater)
	debug << "bigreal::operator > : A->decimal < B->decimal, return false" << endl; 
		return false;
	}
	
	return false;
}

bool operator < (const bigreal& a, const bigreal& b)
{
  if (a == b || a > b)
	  return false;
  else
  	return true;
}

bool operator >= (const bigreal& a, const bigreal& b)
{
	return (!(a<b));
}

bool operator <= (const bigreal& a, const bigreal& b)
{
	return (!(a>b));
}

bigreal abs (const bigreal& a)
{
   bigreal result = a;
   result.negative = false;
   return result;
}

int num_len(const bigreal& a)
{

if (bigreal_control::DEBUG & bigreal_control::num_len)
	debug << "bigreal::num_len : a = " << a << endl;
	
	int length = num_len(a.integral);
	
	if (a.decimal != bigint(0))
	{
		length += num_len(a.decimal) + 1;
	}
	
	return  length;
}

/* carry moves the carry for n significant places in src to dst */
void carry(bigint& src, bigint& dst, int n)
{
	ostringstream oss;
	oss << src;
	string digits = oss.str();
if (bigreal_control::DEBUG & bigreal_control::carry)
	debug << "bigreal::carry: original decimal digits = " << digits << endl;
	
	int carry_length = digits.length() - n;
	
	if (carry_length <0)
	{
		// add leading zeros to the digits to create the correct number of places	
		for (int i=0; i< abs(carry_length); i++)
			digits.insert(0,"0");
if (bigreal_control::DEBUG & bigreal_control::carry)
	debug << "bigreal::carry: adding " << abs(carry_length) << " zeros to give decimal digits = " << digits << endl;
		
		carry_length = 0;
		
	}

	string carry_digits = digits.substr(0,carry_length);

if (bigreal_control::DEBUG & bigreal_control::carry)
	debug << "bigreal::carry: carry_digits = " << carry_digits << endl;

	digits.erase(0,carry_length);

if (bigreal_control::DEBUG & bigreal_control::carry)
	debug << "bigreal::carry: remaining decimal digits = " << digits << endl;	

	istringstream iss_i(carry_digits);
	bigint int_part;
	iss_i >> int_part;

	istringstream iss_d(digits);
	bigint dec_part;
	iss_d >> dec_part;

	dst += int_part;
	src = dec_part;
}





#ifdef TAKEOUT

bigreal::bigreal(const int nn)
{

	int num = nn;
    if (num <0)
    {
		negative = true;
		num *= -1;
    }
    else
		negative = false;

	int digits = (num > BASE ? 2 : 1);

    n.resize(digits);
	n[0] = num%BASE;
	
	if (digits > 1)
		n[1] = num/BASE;
}

bigreal::bigreal(const unsigned long num)
{
	negative = false;

	int digits = (num > BASE ? 2 : 1);
    n.resize(digits);
	n[0] = num%BASE;
	
	if (digits > 1)
		n[1] = num/BASE;
}

bigreal::operator bool() // conversion operator to a bool
{

if (bigreal_control::DEBUG & bigreal_control::bool_conv)
{
	debug << "bigreal::operator bool() : ";
	(*this).dump(debug);
	debug << endl;
}

	if (*this != bigreal(0))
		return true;
	else
		return false;
}

bigreal::operator int() // conversion operator to an int
{
	if (negative)
		return -1 * n[0];
	else
		return n[0];
}



/* sanitize removes any leading zeros from the list of digits */

void bigreal::sanitize()
{

if (bigreal_control::DEBUG & bigreal_control::sanitize)
{
	debug << "bigreal::sanizize : ";
	(*this).dump(debug);
	debug << endl;
}
	size_t size = n.size();

if (bigreal_control::DEBUG & bigreal_control::sanitize)
{
	debug << "bigreal::sanizize : initial size = " << size << endl;
}
	
	for (vector<unsigned int>::reverse_iterator n_ptr=n.rbegin(); n_ptr != n.rend(); n_ptr++)
	{
		if (*n_ptr == 0)
			size--;
		else
			break;	
	}

if (bigreal_control::DEBUG & bigreal_control::sanitize)
{
	debug << "bigreal::sanizize : final size = " << size << endl;
}
	
	n.resize(size); // discards high end zeros.

if (bigreal_control::DEBUG & bigreal_control::sanitize)
{
	debug << "bigreal::sanizize : after sanitizing (*this) = ";
	(*this).dump(debug);
	debug << endl;
}

}


bigreal operator % (const bigreal& a, const bigreal& b)
{
if (bigreal_control::DEBUG & bigreal_control::remainder)
	debug << "bigreal::operator % : " << endl;
    return a - a/b * b;
}

bigreal operator %= (bigreal& a, const bigreal& b)
{
    return a = a%b;
}

bigreal operator >>= (bigreal& a, int n)
{
	vector<unsigned int>::iterator aptr;

if (bigreal_control::DEBUG & bigreal_control::r_shift)
{
	debug << "bigreal::operator >>= : a = " << a  << ", n = " << n << endl;
	debug << "bigreal::operator >>= : a = ";
	a.dump(debug);
	debug << endl;	
}
	 

//cout << "\nin shift operator" << endl;	

	for (int i=0; i < n; i++)
	{
		aptr = a.n.begin();	
		
//cout << "\n\niteration " << i+1 << " of " << n	<< endl;

		if (a.n.size() == 1)
		{
//cout << "\nsize = 1, *aptr = " << *aptr << endl;			
			*aptr >>= 1;
//cout << "\nsize = 1, after right shift *aptr = " << *aptr << endl;			
		}
		else
		{
			// r.n.size() >= 2
			for ( ;aptr != a.n.end()-1; aptr++)
			{
//cout << "\n*aptr = " << *aptr << endl;			
				*aptr >>= 1;
//cout << "\nafter right shift *aptr = " << *aptr << endl;			
				if (*(aptr+1) & 1)
					*aptr |= BIGINT_MSB;
			}
//cout << "\nafter loop *aptr = " << *aptr << endl;			
			*aptr >>= 1;
//cout << "\nafter right shift *aptr = " << *aptr << endl;			
		}
	}

	a.sanitize();
	
if (bigreal_control::DEBUG & bigreal_control::r_shift)
{
	debug << "bigreal::operator >>= : returning a = " << a << endl;
	debug << "bigreal::operator >>= : a = ";
	a.dump(debug);
	debug << endl;
}
	return a;
}

bigreal operator <<= (bigreal& a, int n)
{
	vector<unsigned int>::reverse_iterator aptr;
	
if (bigreal_control::DEBUG & bigreal_control::l_shift)
{
	debug << "bigreal::operator <<= : a = " << a  << ", n = " << n << endl;
	debug << "bigreal::operator <<= : a = ";
	a.dump(debug);
	debug << endl;	
}

	for (int i=0; i < n; i++)
	{
		aptr = a.n.rbegin();
		
 		if (a.n.size() == 1)
		{
//cout << "\naptr = " << *aptr << endl;		
			if (*aptr & BIGINT_MSB)
			{
				a.n.resize(a.n.size()+1);
				*a.n.rbegin() = 1;
//cout << "\nafter resizing aptr = " << *aptr << endl;		
			}
			*aptr <<= 1;
			*aptr &= BIGINT_BITS;
		}
		else
		{
//cout << "\naptr = " << *aptr << endl;		
			if (*aptr & BIGINT_MSB)
			{
				a.n.resize(a.n.size()+1);
				*a.n.rbegin() = 1;
//cout << "\nafter resizing aptr = " << *aptr << endl;		
			}
			// r.n.size() >= 2
			for ( ;aptr != a.n.rend()-1; aptr++)
			{
				*aptr <<= 1;
				*aptr &= BIGINT_BITS;
				if (*(aptr+1) & BIGINT_MSB)
					*aptr |= 1;
			}
			*aptr <<= 1;
			*aptr &= BIGINT_BITS;
		}
	}

if (bigreal_control::DEBUG & bigreal_control::l_shift)
{
	debug << "bigreal::operator <<= : returning a = " << a << endl;
	debug << "bigreal::operator <<= : a = ";
	a.dump(debug);
	debug << endl;
}
	return a;
}


/* The following binary gcd algorithm is Algorithm B from Knuth, The Art of
   Computer Programming vol 2, section 4.5.2.  The notation used here is
   chosen to be consistent with Knuth's.
*/
//int hcf_count;
bigreal gcd (const bigreal& a, const bigreal& b)
{
	int k = 0;
	bigreal t;
	bigreal u(a);
	bigreal v(b);
	
	/* could use abs but faster to set value directly 
	   and avoid the function call
	*/
	u.negative = false;
	v.negative = false;

	if (u == bigreal(0))
		return v;
	
	if (v == bigreal(0))
		return u;

//	hcf_count = 0;
//  cout << endl;
	
if (bigreal_control::DEBUG & bigreal_control::gcd)
	debug << "bigreal::gcd : u = " << u << " v = " << v << endl;

	while (u.even() && v.even()) //B.1
	{
		k++;
		u >>= 1;
		v >>= 1;
if (bigreal_control::DEBUG & bigreal_control::gcd)
{
	debug  << "bigreal::gcd : B.1 both u and v even, k incremented to " << k << endl;
	debug << "bigreal::gcd :   u >>= " << u << " v >>= " << v << endl;
}
	}
	
	if (u.even())               // B.2
	{

if (bigreal_control::DEBUG & bigreal_control::gcd)
	debug << "bigreal::gcd : u remains even" << endl;
		t = u;
	}
	else
		t = v * bigreal(-1);
	
if (bigreal_control::DEBUG & bigreal_control::gcd)
	debug << "bigreal::gcd : B.2 t = " << t << endl;

	do
	{

if (bigreal_control::DEBUG & bigreal_control::gcd)
{
	debug << "bigreal::gcd : t remains non-zero" << t << endl;
//	debug << "\n" << "\t" << t << flush;
}

		while (t.even())       // B.4
		{
			t >>= 1;           // B.3

if (bigreal_control::DEBUG & bigreal_control::gcd)
	debug << "bigreal::gcd : B.4   t remains even t >> = " << t << endl;
	
		}
			
		if (t > bigreal(0))     // B.5
		{
			u = t;

if (bigreal_control::DEBUG & bigreal_control::gcd)
	debug << "bigreal::gcd : B.5   u set to " << u << endl;
	
		}
		else
		{
			v = t * bigreal(-1);
if (bigreal_control::DEBUG & bigreal_control::gcd)
	debug << "bigreal::gcd : B.5   v set to " << v << endl;
	
		}

		t = u-v;              // B.6

if (bigreal_control::DEBUG & bigreal_control::gcd)
	debug << "bigreal::gcd : B.6   t set to " << t << endl;

	} while (t != bigreal(0));
	
	u <<= k;

if (bigreal_control::DEBUG & bigreal_control::gcd)
	debug << "bigreal::gcd : u right shifted " << k << " to give " << u << endl;
	
	return u;
}
#endif
