using namespace std;

#include <fstream>
#include <iostream>
#include <sstream>
#include <list>
#include <vector>
#include <iomanip>

extern unsigned int DEBUG;
extern ofstream		debug;

#include <matrix.h>
#include <decomp.h>
#include <n_limit.h>


n_limit::n_limit(vector<int> z)
{
	n = z.size();
	threshold t;
	t.first = z[n-1];
	
	if (n > 1)
	{
		z.pop_back();
		t.second = new n_limit(z);
	}
	else
	{
		t.second = 0;			
	}
	
	thresholds.push_back(t);
}

/* The copy constructor needs to make a copy of the n_limit
   pointed at by each threshold rather than copy the pointer
*/
n_limit::n_limit(const n_limit& l)
{
	n = l.n;
	list<threshold>::const_iterator tptr = l.thresholds.begin();
	while (tptr != l.thresholds.end())
	{
		threshold t;
		t.first = tptr->first;

		if (n == 1)
			t.second = 0;
		else
			t.second = new n_limit(*tptr->second);
		
		thresholds.push_back(t);
		tptr++;
	}	
}

/* the n_limit destructor has to delete recursively 
any n-limits created under thresholds */
n_limit::~n_limit()
{
	if (n > 1)
	{
		list<threshold>::iterator tptr = thresholds.begin();
		while (tptr != thresholds.end())
		{
			delete tptr->second;					
			tptr++;
		}
	}		
}

void n_limit::print (ostream& os, int indent,string prefix)
{
	list<threshold>::iterator tptr = thresholds.begin();
	
	string indent_string;
	
	for (int i=0; i< indent; i++)
		indent_string += "  ";
		
//	os << prefix << indent_string << "n = " << n << endl;
	
	int count = 0;
	while (tptr != thresholds.end())
	{
		count++;
		os << prefix << indent_string << "  t_" << n << "^" << count << " = " << *tptr << endl;
		if (tptr->second)
			tptr->second->print(os, indent+1,prefix);
		
		tptr++;
	}
}

ostream& operator << (ostream& os, threshold t)
{
	os << '(' << t.first << ',' << t.second << ')';
	return os;
}

ostream& operator << (ostream& os, list<threshold> thresholds)
{
	list<threshold>::iterator tptr = thresholds.begin();
	while (tptr != thresholds.end())
	{
		os << *tptr << ' ';
		tptr++;
	}
	return os;
}

bool operator > (n_limit& l, vector<int> u)
{
if (decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "n_limit::operator >: " << l.n << "-limit = " << l.thresholds << endl;
	//l.print(debug,0,"n_limit::operator >: ");
	debug << "n_limit::operator >: u = ";
	for (unsigned int i=0; i< u.size(); i++)
		debug << u[i] << ' ';
	debug << endl;
}

	if (l.n == 1)
	{
if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "n_limit::operator >: n = 1, returning " << (l.thresholds.begin()->first > u[0]? "true" : "false") << endl;

		return (l.thresholds.begin()->first > u[0]);
	}
	else
	{
			/* find the first threshold > u[l.n-1], if it exists */
			list<threshold>::iterator tptr = l.thresholds.begin();
			bool found = false;
			while (tptr != l.thresholds.end())
			{
				if (tptr->first > u[l.n-1])
				{
					found = true;
					break;
				}
				tptr++;
			}
		
			if (found)
			{
if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "n_limit::operator >: n = " << l.n << ", first threshold > u[" << l.n-1 << "] = " << u[l.n-1] << " is " << tptr->first << ", recursing" << endl;
				return (*tptr->second > u);
			}
			else
			{
if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "n_limit::operator >: n = " << l.n << ", no threshold > u[" << l.n-1 << "] = " << u[l.n-1] << ", returning false" << endl;
				return false;
			}
	}
}

bool operator == (n_limit&l1, n_limit& l2)
{
		
if (decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "n_limit::operator ==: l1: " << l1.n << "-limit = " << l1.thresholds << endl;
	debug << "n_limit::operator ==: l2: " << l2.n << "-limit = " << l2.thresholds << endl;
	//l.print(debug,0,"n_limit::operator >: ");
}

	if (l1.n != l2.n)
	{
if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "n_limit::operator ==: l1 and l2 are not of the same dimension, returning false" << endl;

		return false;
	}

	if (l1.thresholds.size() != l2.thresholds.size())
	{
if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "n_limit::operator ==: l1 and l2 contain a different number of thresholds, returning false" << endl;

		return false;
	}
	
	/* check that the threshold values in l1 and l2 are the same */
	list<threshold>::iterator tptr1 = l1.thresholds.begin();
	list<threshold>::iterator tptr2 = l2.thresholds.begin();

	int threshold_num = 1;
	while (tptr1 != l1.thresholds.end())
	{
		if (tptr1->first != tptr2->first)
		{
if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "n_limit::operator ==: l1 and l2 threshold " << threshold_num << " values differ, returning false" << endl;
			return false;
		}
		
		tptr1++;
		tptr2++;
		threshold_num++;
	}
		
	if (l1.n>1)
	{
		/* check that the threshold (n-1)-limits are the same */
		tptr1 = l1.thresholds.begin();
		tptr2 = l2.thresholds.begin();

		threshold_num = 1;
		while (tptr1 != l1.thresholds.end())
		{
			if (*tptr1->second != *tptr2->second)
			{
if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "n_limit::operator ==: l1 and l2 threshold " << threshold_num << " (n-1)-limits differ, returning false" << endl;
				return false;
			}
		
			tptr1++;
			tptr2++;
			threshold_num++;
		}
	}
	
	
if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "n_limit::operator ==: l1 and l2 equal, returning true" << endl;
	return true;
}

bool operator != (n_limit&l1, n_limit& l2) 
{
	return !(l1 == l2);
}
