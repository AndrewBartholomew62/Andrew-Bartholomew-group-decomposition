/***********************************************************************
			   surd.cpp
			   
Original coding as surd.h A. Bartholomew 9th May 2008

Split into surd.spp and surd.h Mat 2010
	       
	       


***********************************************************************/
using namespace std;

#include <fstream>
#include <cstdlib>
#include <cmath>

extern ofstream debug;
#include <util.h>
#include <surd.h>

void newton_raphson_sqrt (bigreal& x_n, int prime);

Surd::Surd (const int p): prime(p),approx(0),coeff(1)
{
	newton_raphson_sqrt(approx,prime);
}

void Surd::dump(ostream& os) const
{
	os << "Surd(prime = " << prime << ", coeff = " << coeff << ", approx = ";
	approx.dump(os);
	os << ")";
}

istream& operator >> (istream& is, Surd& S)
{
	char ch;
	
	S.coeff = 0;
	is >> S.coeff;

if (surd_control::DEBUG & surd_control::input)	
	debug << "Surd::operator >> : read coeff " << S.coeff << endl;

	if (is.fail())
	{
		/* no coefficent characters read, this could be because there's only a sqrt to read.*/

if (surd_control::DEBUG & surd_control::input)
	debug << "Surd::operator >> : clearing failbit after reading coeff" << endl;

		is.clear();
	}
	
	ch = is.get();
	
	if (is.fail() && S.coeff==0)
	{

if (surd_control::DEBUG & surd_control::input)
	debug << "Surd::operator >> : get failed returning with failbit set" << endl;

		/* there really is nothing ro read, return with the fail bit set */
		return is;
	}
	else if (!(S.coeff==0) && (ch == '-' || ch == '+'))
	{
		/* the coefficient is all there is to read, put back the ch and return */
		S.prime = 1;
		S.approx = 0;				
		is.putback(ch);
		return is;
	}

	if (ch == '-' && is.peek() == 's')
	{
	
		S.coeff = -1;
		ch = is.get();
		
if (surd_control::DEBUG & surd_control::input)
	debug << "Surd::operator >> : minus sign in '-sqrt...' detected, set coeff to -1" << endl;
	
	}
	else if (ch == '+' && is.peek() == 's')
	{
	
		S.coeff = 1;
		ch = is.get();
		
if (surd_control::DEBUG & surd_control::input)
	debug << "Surd::operator >> : plus sign in '+sqrt...' detected, set coeff to 1" << endl;
	
	}
	else if ((is.fail() && !(S.coeff==0)) || ch != 's')
	{

if (surd_control::DEBUG & surd_control::input)
	debug << "Surd::operator >> : no sqrt detected, only a coeff" << endl;

		if (is.fail())
			is.clear();
		else
			is.putback(ch);
			
		return is;
	}

if (surd_control::DEBUG & surd_control::input)	
	debug << "Surd::operator >> : read char (which should be an 's') " << ch << endl;

	if (S.coeff == 0)
	{
		S.coeff = 1;

if (surd_control::DEBUG & surd_control::input)	
	debug << "Surd::operator >> : set coeff to 1 " << endl;
	}

	ch = is.get();

if (surd_control::DEBUG & surd_control::input)	
	debug << "Surd::operator >> : read char " << ch << endl;

	if (ch != 'q')
	{
if (surd_control::DEBUG & surd_control::input)	
	debug << "Surd::operator >> : char not 'q' setting failbit and returning" << endl;
		is.setstate(ios_base::failbit);
		return is;
	}

	ch = is.get();
if (surd_control::DEBUG & surd_control::input)	
	debug << "Surd::operator >> : read char " << ch << endl;
	if (ch != 'r')
	{
if (surd_control::DEBUG)	
	debug << "Surd::operator >> : char not 'r' setting failbit and returning" << endl;
		is.setstate(ios_base::failbit);
		return is;
	}

	ch = is.get();
if (surd_control::DEBUG & surd_control::input)	
	debug << "Surd::operator >> : read char " << ch << endl;
	if (ch != 't')
	{
if (surd_control::DEBUG & surd_control::input)	
	debug << "Surd::operator >> : char not 't' setting failbit and returning" << endl;
		is.setstate(ios_base::failbit);
		return is;
	}

	ch = is.get();
if (surd_control::DEBUG & surd_control::input)	
	debug << "Surd::operator >> : read char " << ch << endl;
	if (ch != '(')
	{
if (surd_control::DEBUG & surd_control::input)	
	debug << "Surd::operator >> : char not '(' setting failbit and returning" << endl;
		is.setstate(ios_base::failbit);
		return is;
	}

	int p;
	is >> p;

if (surd_control::DEBUG & surd_control::input)	
	debug << "Surd::operator >> : read prime " << p << endl;

	ch = is.get();
if (surd_control::DEBUG & surd_control::input)	
	debug << "Surd::operator >> : read char " << ch << endl;
	if (ch != ')')
	{
		is.setstate(ios_base::failbit);
if (surd_control::DEBUG & surd_control::input)	
	debug << "Surd::operator >> : char not ')' setting failbit and returning" << endl;
		return is;
	}

	S.prime = p;
	S.approx = 0;
	newton_raphson_sqrt(S.approx,p);

if (surd_control::DEBUG & surd_control::input)	
{
	debug << "Surd::operator >> : read Surd ";
	S.dump(debug);
	debug << endl;
}

		
	return is;
}

void surd::dump(ostream& os) const
{

	os << "surd(list length = " << terms.size() << endl;

	list<Surd>::const_iterator ptr = terms.begin();

	while (ptr != terms.end())
	{
		os << "     ";
		ptr->dump(os);
		os << endl;
		
		ptr++;
	}
	os << "     end of list)";
}

bigreal surd::estimate()
{
	bigreal result;
	list<Surd>::iterator ptr = terms.begin();

	while (ptr != terms.end())
	{
		if (ptr->prime == 1)
		{
			bigreal bigreal_coeff = bigreal(1) * ptr->coeff;
			result += bigreal_coeff;
		}
		else
			result += ptr->approx * ptr->coeff;						

		ptr++;
	};
	
	return result;
}

ostream& operator << (ostream& os, const surd& s)
{
	if (s.terms.size() == 0)
		os << "0";
	else
	{
		list<Surd>::const_iterator ptr = s.terms.begin();

		while (ptr != s.terms.end())
		{
			if (ptr != s.terms.begin() && ptr->coeff > 0)
					os << '+';
			else if (ptr->prime != 1 && ptr->coeff == rational<int>(-1))
					os << '-';

			if (ptr->prime == 1 || (ptr->coeff != rational<int>(1) && ptr->coeff != rational<int>(-1)))
				os << ptr->coeff;
			
/*			if (ptr->coeff == rational<int>(-1))
			{
				os << "-";
				if (ptr->prime == 1)
					os << "1";
			}
			else if (abs(ptr->coeff) != rational<int>(1) || ptr->prime == 1)
				os << ptr->coeff;
*/
			if (ptr->prime != 1)
				os << "sqrt(" << ptr->prime << ")";

			ptr++;
		};
	}
		
	return os;
}

surd operator + (const surd& a, const surd& b)
{

if (surd_control::DEBUG & surd_control::add)	
{
	debug << "surd::operator + : a = ";
	a.dump(debug);
	debug << endl;
	debug << "surd::operator + : b = ";
	b.dump(debug);
	debug << endl;
}

	surd result(a);
	
	list<Surd>::const_iterator ptr = b.terms.begin();
	while (ptr != b.terms.end())
	{
		result += *ptr;
		ptr++;
	};
	
	return result;
}

surd operator + (const surd& s, const Surd& t)
{

if (surd_control::DEBUG & surd_control::add)	
{
	debug << "surd::operator + : s = ";
	s.dump(debug);
	debug << endl;
	debug << "surd::operator + : t = ";
	t.dump(debug);
	debug << endl;
}

	surd result(s);
	Surd term(t);
	bool not_found = true;
	
	list<Surd>::iterator ptr = result.terms.begin();
	
	while (ptr != result.terms.end())
	{
		if (ptr->prime == term.prime)
		{
			ptr->coeff += term.coeff;

if (surd_control::DEBUG & surd_control::add)	
	debug << "surd::operator + :  t found , new coefficient = " << ptr->coeff << endl;	
		
			if (ptr->coeff == rational<int>(0))
			{
if (surd_control::DEBUG & surd_control::add)
{
	debug << "surd::operator + :  erasing ";
	if (ptr->prime == 1)
		debug << "unit";
	else
		debug << "sqrt(" << ptr->prime << ")";
	debug << " term" << endl;		
}
				ptr = result.terms.erase(ptr);
			}			
			else if (term.approx.decimal_places() > ptr->approx.decimal_places())
			{			
if (surd_control::DEBUG & surd_control::add)	
	debug << "surd::operator + :  copying approximation from t" << endl;		
	
				ptr->approx = term.approx;
			}
				
			not_found = false;
			break;
		}
		
		ptr++;
	};	
	
	if (not_found)
	{
if (surd_control::DEBUG & surd_control::add)	
	debug << "surd::operator + :  t is not in s, adding to list" << endl;		

		result.terms.push_back(term);
	}
	
	return result;
}

surd operator += (surd& s, const Surd& t)
{

if (surd_control::DEBUG & surd_control::add)	
{
	debug << "surd::operator += : s = ";
	s.dump(debug);
	debug << endl;
	debug << "surd::operator += : t = ";
	t.dump(debug);
	debug << endl;
}
	s = s + t;
	return s;
}
	
surd operator - (const surd& a, const surd& b)
{

if (surd_control::DEBUG & surd_control::subtract)	
{
	debug << "surd::operator - : a = ";
	a.dump(debug);
	debug << endl;
	debug << "surd::operator - : b = ";
	b.dump(debug);
	debug << endl;
}
	surd result(a);
	
	list<Surd>::const_iterator ptr = b.terms.begin();
	while (ptr != b.terms.end())
	{
		result -= *ptr;
		ptr++;
	};
	
	return result;
}

surd operator - (const surd& s, const Surd& t)
{
	
if (surd_control::DEBUG & surd_control::subtract)	
{
	debug << "surd::operator - : s = ";
	s.dump(debug);
	debug << endl;
	debug << "surd::operator - : t = ";
	t.dump(debug);
	debug << endl;
}

	surd result(s);
	Surd term(t);
	bool not_found = true;

	list<Surd>::iterator ptr = result.terms.begin();
	
	while (ptr != result.terms.end())
	{
		if (ptr->prime == term.prime)
		{
			ptr->coeff -= term.coeff;		

if (surd_control::DEBUG & surd_control::subtract)	
	debug << "surd::operator - :  t found , new coefficient = " << ptr->coeff << endl;		

			if (ptr->coeff == rational<int>(0))
			{
if (surd_control::DEBUG & surd_control::subtract)
{
	debug << "surd::operator - :  erasing ";
	if (ptr->prime == 1)
		debug << "unit";
	else
		debug << "sqrt(" << ptr->prime << ")";
	debug << " term" << endl;		
}
				ptr = result.terms.erase(ptr);
			}
			else if (term.approx.decimal_places() > ptr->approx.decimal_places())
			{
if (surd_control::DEBUG & surd_control::subtract)	
	debug << "surd::operator - :  copying approximation from t" << endl;		

				ptr->approx = term.approx;
			}
				
			not_found = false;
			break;
		}
		
		ptr++;
	};	
	
	if (not_found)
	{
if (surd_control::DEBUG & surd_control::subtract)	
	debug << "surd::operator - :  t is not in s, inverting and adding to list" << endl;		

		term.coeff *= rational<int>(-1);
		result.terms.push_back(term);
	}
	
	return result;
}

surd operator -= (surd& s, const Surd& t)
{

if (surd_control::DEBUG & surd_control::subtract)	
{
	debug << "surd::operator -= : s = ";
	s.dump(debug);
	debug << endl;
	debug << "surd::operator -= : t = ";
	t.dump(debug);
	debug << endl;
}
	s = s - t;
	return s;
}

surd operator += (surd& a, const surd& b)
{

if (surd_control::DEBUG & surd_control::add)	
{
	debug << "surd::operator += : a = ";
	a.dump(debug);
	debug << endl;
	debug << "surd::operator += : b = ";
	b.dump(debug);
	debug << endl;
}
    return a = a+b;
}

surd operator -= (surd& a, const surd& b)
{

if (surd_control::DEBUG & surd_control::subtract)	
{
	debug << "surd::operator - : a = ";
	a.dump(debug);
	debug << endl;
	debug << "surd::operator - : b = ";
	b.dump(debug);
	debug << endl;
}
    return a = a-b;
}

bool operator == (const surd& a, const surd& b)
{

if (surd_control::DEBUG & surd_control::equal)	
{
	debug << "surd::operator == : a = ";
	a.dump(debug);
	debug << endl;
	debug << "surd::operator == : b = ";
	b.dump(debug);
	debug << endl;
}

	if (b.terms.size() == 0)
	{
		return (!a.terms.size());
	}
	else
	{
		surd diff = a-b;
		surd zero;
		return (diff == zero);
	}
}

bool operator != (const surd& a, const surd& b)
{

if (surd_control::DEBUG & surd_control::equal)	
{
	debug << "surd::operator != : a = ";
	a.dump(debug);
	debug << endl;
	debug << "surd::operator != : b = ";
	b.dump(debug);
	debug << endl;
}
	return (!(a == b));
}

bool operator > (const surd& aa, const surd& b)
{

if (surd_control::DEBUG & surd_control::greater)	
{
	debug << "surd::operator > : aa = ";
	aa.dump(debug);
	debug << endl;
	debug << "surd::operator > : b = ";
	b.dump(debug);
	debug << endl;
}
	
	if (b.terms.size() == 0)
	{
if (surd_control::DEBUG & surd_control::greater)	
	debug << "surd::operator > : b is zero, testing sign of a" << endl;
	
		if (aa.terms.size() == 0)
		{
if (surd_control::DEBUG & surd_control::greater)	
	debug << "surd::operator > : a is zero, returning false" << endl;

			return false;
		}
		else
		{	
			/* since a includes some terms it is either > or < zero.  First check that 
			   it contains some prime root, since if it does not, it is rational and we can determine
			   its sign directly.
			*/
			bool aa_rational = true;

			list<Surd>::const_iterator cptr = aa.terms.begin();
			while (cptr != aa.terms.end())
			{
				if (cptr->prime != 1)
				{
					aa_rational = false;
					break;
				}
				cptr++;
			};
			
			if (aa_rational)
			{
if (surd_control::DEBUG & surd_control::greater)	
{
	debug << "surd::operator > : a is rational, " << "coeff = " << aa.terms.begin()->coeff 
	      << ", returning " << (aa.terms.begin()->coeff > 0 ? "true": "false") << endl;
}

				return (aa.terms.begin()->coeff > 0); // the prime is 1
			}
			
			/* So, here we know that a is irrational and non-zero
			
			   If one is trying to decide whether a non-zero linear combination 

				   r_1 p_1 +/- ... +/- r_k p_k    (1)

			   of k distinct prime roots p_i over rationals r_i is positive, and one has 
			   approximations q_i to p_i each of which is accurate to n decimal places 
			   then each approximation is accurate to +/- 10^{-n}.  This means that 

				   r_1 q_1 +/- ... +/- r_k q_k    (2)

			   is accurate to 

				   (r_1 +/- ... +/- r_k) 10^{-n}.

			   Now (r_1 +/- ... +/- r_k) comprises m = 1+[log_10(r_1 +/- ... +/- r_k)] digits, 
			   where [] is the integral part.  Hence (assuming large enough n) if the n-m th 
			   decimal place of (2) is not 0 or 9 then (2) is accurate to n-m places.

			   Thus, to decide whether (1) is positive or negative, evaluate (2) and 
			   find the first non-zero digit that is not 9 and for which the q_i are all accurate 
			   to at least m additional places (use Newton-Raphson as required to ensure this).  

			   If such a digit exists, the sign of (2) is accurate.			
			*/

			surd a(aa);
			
			bool not_found = true;
			bool result;
	
			do 
			{

if (surd_control::DEBUG & surd_control::greater)	
{
	debug << "surd::operator > : testing approximation given by: " << endl;
	a.dump(debug);
	debug << endl;
	
}
				list<Surd>::iterator ptr = a.terms.begin(); 
				int accuracy = ptr->approx.decimal_places();

if (surd_control::DEBUG & surd_control::greater)	
	debug << "surd::operator > : initial accuracy = " << accuracy << endl;

				bigreal digits_check;
				bigreal estimate;

				while (ptr != a.terms.end())
				{
					if (ptr->approx.decimal_places() != 0 && ptr->approx.decimal_places() < accuracy)
					{
						accuracy = ptr->approx.decimal_places();

if (surd_control::DEBUG & surd_control::greater)	
	debug << "surd::operator > : initial accuracy reduced to = " << accuracy << endl;
					}

					if (ptr->prime == 1)
					{
						bigreal bigreal_coeff = bigreal(1) * ptr->coeff;
						estimate += bigreal_coeff;
					}
					else
						estimate += ptr->approx * ptr->coeff;						

if (surd_control::DEBUG & surd_control::greater)	
	debug << "surd::operator > : estimate increased to = " << estimate << endl;
	
					digits_check += ptr->coeff;

if (surd_control::DEBUG & surd_control::greater)	
	debug << "surd::operator > : digits_check increased to = " << digits_check << endl;

					ptr++;
				};

if (surd_control::DEBUG & surd_control::greater)	
{
	debug << "surd::operator > : estimate = " << estimate << endl;
	debug << "surd::operator > : digits_check = " << digits_check << endl;
}	

				/* determine the number of digits in estimate */
				unsigned int offset = 1;

				while (digits_check.integral_part() != bigint(0))
				{
					digits_check.integral_part() /= bigint(10);
					offset++;
				}

if (surd_control::DEBUG & surd_control::greater)	
	debug << "surd::operator > : accuracy digits offset = " << offset << endl;

				/* check whether we have sufficient digits in the estimate to decide its sign */
				ostringstream oss;
				oss << estimate.decimal_part();

				string places = oss.str();

				if (places.length() > offset)
				{
					for (unsigned int i=0; i< places.length()-offset; i++)
					{
						if (places[i] != '0' && places[i] != '9')
						{
							not_found = false;
							result = (estimate > 0);

if (surd_control::DEBUG & surd_control::greater)	
	debug << "surd::operator > : decimal digit " << i << " (from zero) is not 0 or 9" << endl;

							break;
						}
					}
				}
				
				if (not_found)
				{

if (surd_control::DEBUG & surd_control::greater)	
	debug << "surd::operator > : no determining decimal found, extending accuracy of approximations" << endl;

					/* increase the accuracy of each approximation */
					ptr = a.terms.begin();
					while (ptr != a.terms.end())
					{
						newton_raphson_sqrt(ptr->approx,ptr->prime);
						ptr++;
					};
					
				}
				
			} while (not_found);
			
if (surd_control::DEBUG & surd_control::greater)	
	debug << "surd::operator > : result = " << (result? "true" : "false") << endl;

			return result;
		}			
	}
	else
	{
if (surd_control::DEBUG & surd_control::greater)	
	debug << "surd::operator > : b non-zero, creating difference " << endl;
	
		surd diff = aa-b;
		surd zero;
		return (diff > zero);
	}
}

bool operator < (const surd& a, const surd& b)
{

if (surd_control::DEBUG & surd_control::greater)	
{
	debug << "surd::operator < : a = ";
	a.dump(debug);
	debug << endl;
	debug << "surd::operator < : b = ";
	b.dump(debug);
	debug << endl;
}

  if (a == b || a > b)
	  return false;
  else
  	return true;
}

bool operator >= (const surd& a, const surd& b)
{

if (surd_control::DEBUG & surd_control::greater)	
{
	debug << "surd::operator >= : a = ";
	a.dump(debug);
	debug << endl;
	debug << "surd::operator >= : b = ";
	b.dump(debug);
	debug << endl;
}
	return (!(a<b));
}

bool operator <= (const surd& a, const surd& b)
{

if (surd_control::DEBUG & surd_control::greater)	
{
	debug << "surd::operator <= : a = ";
	a.dump(debug);
	debug << endl;
	debug << "surd::operator <= : b = ";
	b.dump(debug);
	debug << endl;
}
	return (!(a>b));
}

istream& operator >> (istream& is, surd& s)
{
	s.terms.clear();
	
	Surd first_term;
	is >> first_term;
	
	if (is.fail())
		return is;
	else if (first_term.coeff != rational<int>(0))
	{

if (surd_control::DEBUG & surd_control::input)
{
	debug << "surd::operator >> : read first term ";
	first_term.dump(debug);
	debug << endl;
}
		s.terms.push_back(first_term);
	}
		
	char ch;
	
	while ((ch = is.get()))
	{
		if (ch != '+' && ch != '-')
		{
			is.putback(ch);
			break;
		}
		else
		{
			Surd term;
			is.putback(ch);
			is >> term;

if (surd_control::DEBUG & surd_control::input)
{
	debug << "surd::operator >> : read term ";
	term.dump(debug);
	debug << endl;
}

			if (term.coeff != rational<int>(0))
				s.terms.push_back(term);
		}
	}
	
	return is;		
}


void newton_raphson_sqrt (bigreal& x_n, int prime)
{
	if (x_n == 0)
	{
		/* provide initial approximation */
		double root_p = sqrt(prime);

if (surd_control::DEBUG & surd_control::newton)	
	debug << "Surd::newton_raphson: root_" << prime << " = " << root_p << endl;
	
		unsigned int int_part = int(root_p);
		x_n=bigreal(int_part);
		root_p -= int_part;
		x_n.decimal_part() = int(1000000000*root_p);
		x_n.set_decimal_places(9);
	
if (surd_control::DEBUG & surd_control::newton)	
{
	debug << "Surd::newton_raphson: initial approx for sqrt(" << prime << ") = " << x_n << ": ";
	x_n.dump(debug);	
	debug << endl;
}

		newton_raphson_sqrt(x_n,prime);

if (surd_control::DEBUG & surd_control::newton)	
{
	debug << "Surd::newton_raphson: final approx for sqrt(" << prime << ") = " << x_n << ": ";
	x_n.dump(debug);	
	debug << endl;
}
	}
	else
	{

/*
	bigreal x_n_squared = x_n * x_n;

if (surd_control::DEBUG >= surd_control::SUMMARY)
{
	debug << "newton_raphson: x_n squared = ";
	x_n_squared.dump(debug);
	debug << endl;
}
	bigreal numerator = x_n_squared - bigreal(prime);
	
if (surd_control::DEBUG >= surd_control::SUMMARY)
{
	debug << "newton_raphson: numerator = ";
	numerator.dump(debug);
	debug << endl;
}
	
	bigreal denominator = bigreal(2) * x_n;
	
if (surd_control::DEBUG >= surd_control::SUMMARY)
{
	debug << "newton_raphson: denominator = ";
	denominator.dump(debug);	
	debug << endl;
}
	
	bigreal quotient = numerator / denominator;
	
if (surd_control::DEBUG >= surd_control::SUMMARY)
{
	debug << "newton_raphson: quotient = ";
	quotient.dump(debug);	
	debug << endl;
}

	bigreal result = x_n - quotient;

if (surd_control::DEBUG >= surd_control::SUMMARY)
{
	debug << "newton_raphson: result = ";
	result.dump(debug);	
	debug << endl;
}

	return result;
*/				
		x_n = x_n - (x_n * x_n - bigreal(prime))/(bigreal(2) * x_n);

if (surd_control::DEBUG & surd_control::newton)	
{
	debug << "Surd::newton_raphson: iterated approximation for sqrt(" << prime << ") = " << x_n << ": ";
	x_n.dump(debug);	
	debug << endl;
}
	}
}

bigreal operator * (const bigreal& a, rational<int> r)
{

if (surd_control::DEBUG & surd_control::times)	
{
	debug << "bigreal::operator * : a = ";
	a.dump(debug);
	debug << endl;
	debug << "bigreal::operator * : r = " << r << endl;
}
	bigreal result(a);
	result *= bigreal(r.getn());
	result /= bigreal(r.getd());
	return result;
}

bigreal operator += (bigreal& a, rational<int> r)
{

if (surd_control::DEBUG & surd_control::add)	
{
	debug << "bigreal::operator += : a = ";
	a.dump(debug);
	debug << endl;
	debug << "bigreal::operator += : r = " << r << endl;
}

	bigreal temp(r.getn());
	temp /= bigreal(r.getd());
	
	return a += temp;

}

