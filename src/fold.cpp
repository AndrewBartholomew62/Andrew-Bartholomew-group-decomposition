/***********************************************************************
The fold function performs a folding sequence on a list of 2-cells whose
attaching maps are expressed as signed integers rather than as words.

The fold function was derived from the programme jsj-fold (that uses character
based attaching maps), incorporating the quadratic_pair stop condition.
   
					Andrew Bartholomew May 2010
	


***********************************************************************/
using namespace std;
#include <fstream>
#include <iostream>
#include <cstdlib>
#include <cstring>
#include <iomanip>
#include <map>
#include <vector>

#include <util.h>


extern ofstream debug;
extern ofstream output;

#include <matrix.h>
#include <surd.h>
#include <decomp.h>


typedef pair<list<int>,list<int> > two_cell;


int 		max_folds = 5;
int			num_folds = 0;
bool		QUADRATIC_WORD = false;
bool folding_required(two_cell& tc, vector<surd>& edge_length);
void substitute(vector<int>& word, int a, int b, bool a_away, bool b_away);

#ifdef USE_FOLD
list<two_cell_record> word_list;

void dump(ostream&, list<two_cell>&);
string inv_word(vector<int>& w);
bool quadratic_pair(two_cell& tc);
bool quadratic_word(two_cell &tc);
#endif

//void initialize_word_list(two_cell& tc, map<char,surd>& edge_length);


void fold(matrix<int>& Two_Cells, matrix<int>& track_basis)
{
#ifdef USE_FOLD	
	int n = track_basis.numcols();
	
	map<char,surd> edge_length;
	list<two_cell> two_cells;

	/* create the list of two_cells from Two_Cells */
	for (int i=o; i< Two_Cells.numrows(); i++)
	{
		list<int> top_word;
		top_word.push_back(Two_Cells[i][0]);

		list<int> bottom_word;
		bottom_word.push_back(Two_Cells[i][1]);
		bottom_word.push_back(Two_Cells[i][2]);
		
		two_cells.push_back(two_cell(top_word,bottom_word));
	}
	
	vector<surd> edge_length(n);
	
	/* calculate initial edge lengths.  These are determined by taking the elements
	   of the track_basis u1, ..., u_n and calculating \Sum alpha_i u_i where
	   \alpha_1 = 1 and \alpha_j = sqrt(p_j) for j>1, where the p_j are distinct
	   primes.  We read the primes from an input_file called "primes" so we can 
	   change or add to the prime list as required.
	*/
	
	ifstream input;
	input.open("primes");
	if (!input)
	{
		cout << "\nError cannot find the list of primes.  The programme requires a file in the current\n";
		cout << "directory containing a list of primes to calculate the initial edge lengths for the\n";
		cout << "folding sequence." << endl;
		exit(0);
	}
	else
	{
		for (int j=0; j < n; j++)
			edge_length[i] = surd(1);
	}
#endif

#ifdef STRING_VERSION	
	cout << "\ninitial two_cell list:" << endl;
	output << "\ninitial two_cell list:" << endl;
	
	lptr = two_cells.begin();
	
	while (lptr != two_cells.end())
	{
		cout << "  w = " << lptr->first << ", w' = " << lptr->second << ", \\bar w' = " << invstr(lptr->second)<< endl;
		output << "  w = " << lptr->first << ", w' = " << lptr->second << ", \\bar w' = " << invstr(lptr->second)<< endl;
		lptr++;
	};

	cout << "\ninitial edge lengths:\n";
	output << "\ninitial edge lengths:\n";

	map<char,surd>::iterator eptr(edge_length.begin());
	while (eptr != edge_length.end())
	{
		cout << "  " << eptr->first << " = " << eptr->second << " = " << eptr->second.estimate() << endl;
		output << "  " << eptr->first << " = " << eptr->second << " = " << eptr->second.estimate() << endl;
		eptr++;
	};
	cout << endl;
	output<< endl;
	
if (decomp_control::DEBUG >= decomp_control::SUMMARY)
{
	debug << "jsj::fold: initial two_cell list:" << endl;
	list<two_cell>::iterator lptr = two_cells.begin();
	
	while (lptr != two_cells.end())
	{
		debug << "             w = " << lptr->first << ", w' = " << lptr->second << endl;
		lptr++;
	};
	debug << endl;

	debug << "jsj::fold: initial edge lengths: ";
	map<char,surd>::iterator eptr(edge_length.begin());
	while (eptr != edge_length.end())
	{
		debug << eptr->first << '=' << eptr->second << ' ';
		eptr++;
	};
	debug << endl;
}	
	
	bool finished_folding;
	
	do
	{
		finished_folding = true;

		list<two_cell>::iterator lptr = two_cells.begin();
		
		while (lptr != two_cells.end())
		{
			if (lptr->first == "folded" || lptr->first.find("quadratic") != string::npos)
			{
				lptr++;
				continue;
			}

if (decomp_control::DEBUG >= decomp_control::SUMMARY)
	debug << "jsj::fold: considering folding 2-cell with w=" << lptr->first << ", w' =" << lptr->second << endl;

			num_folds = 0;
			QUADRATIC_WORD = false;
			
			while (folding_required(*lptr,edge_length))
			{

if (decomp_control::DEBUG >= decomp_control::SUMMARY)		
	debug << "jsj::fold: folding " << lptr->first << ',' << lptr->second << endl; 

				cout << "\nfolding 2-cell w=" << lptr->first << ", w' =" << lptr->second << ", \\bar w' =" << invstr(lptr->second) << endl;
				output << "\nfolding 2-cell w=" << lptr->first << ", w' =" << lptr->second << ", \\bar w' =" << invstr(lptr->second) << endl;

				finished_folding = false;
				

				/* First check that the top and bottom lengths are equal, and that there
				   are no two vertices on the same leaf of the foliation determined by the
				   lengths.  If there are, replace the 2-cell by the two 2-cells produced
				   by contracting the leaf
				*/
				string::iterator tptr = lptr->first.begin();
				string::reverse_iterator bptr = lptr->second.rbegin();
				surd top_length;
				surd bottom_length;
				string::iterator top_common;
				string::reverse_iterator bottom_common;
				bool common_found = false;
				bool leaf_contraction_required = false;
				
				do
				{
					if (*tptr == '-')
						tptr++;
					
					/* the bottom pointer is a reverse iterator, so we're moving backwards along a word
					   that might contain minus signs, therefore check for them *after* processing the edge
					*/

if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << "jsj::fold:   next top edge = " << *tptr << " bottom edge = " << *bptr << endl;

					top_length += edge_length[*tptr];
					bottom_length += edge_length[*bptr];

if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << "jsj::fold:   cumulative lengths: top = " << top_length << " bottom = " << bottom_length << endl;
					
					if (top_length == bottom_length)
					{
					
						if (!common_found)
						{
							common_found = true;
							top_common = tptr;
							bottom_common = bptr;
						}
						else
							leaf_contraction_required = true; // the first common vertex is an intermediate vertex
					}
						
					tptr++;
					bptr++;
					if (*bptr == '-')
						bptr++;

				} while (tptr != lptr->first.end() && bptr != lptr->second.rend());
							
				while (tptr != lptr->first.end())
				{
if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << "jsj::fold:   further top edges to consider" << endl;
					
					if (*tptr == '-')
						tptr++;

if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << "jsj::fold:   next top edge = " << *tptr << endl;
						
					top_length += edge_length[*tptr];					
					tptr++;

if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << "jsj::fold:   cumulative top length = " << top_length << endl;

				}
				
				while (bptr != lptr->second.rend())
				{
if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "jsj::fold:   further bottom edges to consider" << endl;
	debug << "jsj::fold:   next bottom edge = " << *bptr << endl;
}
						
					bottom_length += edge_length[*bptr];					
					bptr++;

					if (*bptr == '-')
						bptr++;

if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << "jsj::fold:   cumulative bottom length = " << bottom_length << endl;
				}
				
				if (top_length != bottom_length)
				{

if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "\nError! Folding sequence being attempted on 2-cell with different top and bottom lengths"	<< endl;
	debug << " w = " << lptr->first << ", w' = " << lptr->second << endl;
	debug << "edge lengths: ";
	map<char,surd>::iterator eptr(edge_length.begin());
	while (eptr != edge_length.end())
	{
		debug << eptr->first << '=' << eptr->second << ' ';
		eptr++;
	};	
	debug << endl;
	debug << " w length = " << top_length << " w' length = " << bottom_length << endl;
}

					cout << "\nError! Folding sequence being attempted on 2-cell with different top and bottom lengths"	<< endl;
					cout << " w = " << lptr->first << ", w' = " << lptr->second << endl;
					cout << "edge lengths: ";
					output << "\nError! Folding sequence being attempted on 2-cell with different top and bottom lengths"	<< endl;
					output << " w = " << lptr->first << ", w' = " << lptr->second << endl;
					output << "edge lengths: ";

					map<char,surd>::iterator eptr(edge_length.begin());
					while (eptr != edge_length.end())
					{
						cout << eptr->first << '=' << eptr->second << ' ';
						output << eptr->first << '=' << eptr->second << ' ';
						eptr++;
					};	
					cout << endl;
					cout << " w length = " << top_length << " w' length = " << bottom_length << endl;
					output << endl;
					output << " w length = " << top_length << " w' length = " << bottom_length << endl;

					exit(0);				
				}
				else
				{
if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "jsj::fold:   compatible lengths: top = " << top_length << " bottom = " << bottom_length << endl;
	debug << "jsj::fold:   comon_found = " << (common_found? "true" : "false") << endl;
	debug << "jsj::fold:   leaf_contraction_required = " << (leaf_contraction_required? "true" : "false") << endl;
}
															
					if (leaf_contraction_required)
					{
						/* there is at least one more edge in both the top and the bottom */
						top_common++;
						bottom_common++;					
						if (*bottom_common == '-')
							bottom_common++;

						/* having adjusted top_common and bottom_common the initial vertices of the edges they 
						   lie on the same leaf of the foliation.  We may adjust the 2-cell to contract the leaf 
						   connecting the two vertices and add the new 2-cell to the list, adjusting
						   the attaching words of the current 2-cell accordingly.  Note that since the bottom word is 
						   specified in the reverse direction, the bottom word for the new 2-cell is the first part 
						   of the current word.

						*/
				
						string new_top = lptr->first.substr(top_common-lptr->first.begin());
						string new_bottom = lptr->second.substr(0, lptr->second.length()-(bottom_common-lptr->second.rbegin()));

						cout << "\nfound a vertex on the same leaf of the foliation" << endl;
						cout << "  creating 2-cell w = " << new_top << ", w' = " << new_bottom << ", \\bar w' =" << invstr(new_bottom) << endl;
						output << "\nfound a vertex on the same leaf of the foliation" << endl;
						output << "  creating 2-cell w = " << new_top << ", w' = " << new_bottom << ", \\bar w' =" << invstr(new_bottom) << endl;

if (decomp_control::DEBUG >= decomp_control::BASIC)		
{
	debug << "jsj::fold: found common vertex in top and bottom foliation, new_top = " << new_top
		  << " new_bottom = " << new_bottom << endl; 
}

						/* check again for a Type 3 fold, in either the new words or the remainder of the old */
						if (word_length(new_top) == 1 && word_length(new_bottom) == 1)
						{
							cout << "  Type III fold detected in new 2-cell, skipping 2-cell" << endl;
							output << "  Type III fold detected in new 2-cell, skipping 2-cell" << endl;

if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << "jsj::fold:   Type III fold detected, 2-cell w=" << new_top << ", w' = " << new_bottom << " not added to list" << endl;
						}
						else
							two_cells.push_back(two_cell(new_top,new_bottom));
						
						lptr->first.erase(top_common-lptr->first.begin(),string::npos);
						lptr->second.erase(0,lptr->second.length()-(bottom_common-lptr->second.rbegin()));					

if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "jsj::fold: adjusted two_cell list:" << endl;
	list<two_cell>::iterator lptr = two_cells.begin();
	
	while (lptr != two_cells.end())
	{
		debug << "jsj::fold:   w = " << lptr->first << ", w' = " << lptr->second << ", \\bar w' =" << invstr(lptr->second) << endl;
		lptr++;
	};
}	
						cout << "  continuing to fold 2-cell w = " << lptr->first << ", w' = " << lptr->second << ", \\bar w' =" << invstr(lptr->second) << endl;
						output << "  continuing to fold 2-cell w = " << lptr->first << ", w' = " << lptr->second << ", \\bar w' =" << invstr(lptr->second) << endl;

if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << "jsj::fold: continuing to fold 2-cell w = " << lptr->first << ", w' = " << lptr->second << ", \\bar w' =" << invstr(lptr->second) << endl;

						/* Since the search above detects the first common vertex in the original words, we know
						   that the the remainder of the old words do not contain a common vertex themselves.  We
						   check here whether they represent a Type III fold, or whether they constitute a 
						   quadratic pair.  if neither of these conditions hold, we will need to perform a Type I
						   fold on these words, which is done when we leave the enclosing conditional clause.
						   
						   After checking whether our folding can stop immediately we also check whether the
						   remaining words constitute a quadratic attaching map.  If so we initialize the word_list
						   to start recording the words after each fold so we may detect the quadratic_repetition
						   stop condition.
						*/
						if (word_length(lptr->first) == 1 && word_length(lptr->second) == 1)
						{

							cout << "  contracting foliation leaf leaves a Type III fold, removing 2-cell" << endl;
							output << "  contracting foliation leaf leaves a Type III fold, removing 2-cell" << endl;

if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << "jsj::fold:   Type III fold detected, removing 2-cell w=" << lptr->first << ", w' = " << lptr->second << endl;

							lptr->first = "folded";
							lptr->second = "away";					
							break;
						}
						else if (quadratic_pair(*lptr))
						{
if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << "jsj::fold:   remaining 2-cell is a quadratic pair, no further folding required" <<  endl;

							/* mark the entry in two_cell as a quadratic pair */
							lptr->first += lptr->second;
							lptr->second = lptr->first;
							lptr->first = "quadratic pair";

							cout << "  identified as a quadratic pair, no further folding required" << endl;
							output << "  identified as a quadratic pair, no further folding required" << endl;
							break;
						}
						else if (quadratic_word(*lptr))
						{
							initialize_word_list(*lptr,edge_length);
						}
					}
					else if (word_length(lptr->first) == 1 && word_length(lptr->second) == 1) // Type III
					{
						cout << "  Type III fold detected, removing 2-cell" << endl;
							output << "  Type III fold detected, removing 2-cell" << endl;

if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << "jsj::fold:   Type III fold detected, removing 2-cell w=" << lptr->first << ", w' = " << lptr->second << endl;

						/* we'd like to remove the two_cell from the list but the iterators are not
						   valid after a call to list.erase(), list.remove() etc. , so we mark the
						   entry instead
						*/
						lptr->first = "folded";
						lptr->second = "away";					

						break;
					}													
					
					/* Perform a Type I fold.  If the edges at the corner of the 2-cell are a and b with a < b
					   we fold a onto b.  We have to deal with the cases where a and b are oriented away from or
					   towards the corner.

					   When we fold a onto b it has the effect of subdividing b into two edges, which we also 
					   call a and b, where the 'new' b has length equal to the difference between the 'old' b
					   and a.

					   After the fold, the attaching words of all the 2-cells other than the folded 2-cell have 
					   to be adjusted to replace the 'old' b with the 'new' b, as follows:

					   a oriented away from, b oriented away from the corner:  b -> ab		
					   a oriented away from, b oriented towards the corner:    b -> b-a
					   a oriented towards, b oriented away from the corner:    b -> -ab
					   a oriented towwards b oriented towards the corner:      b -> ba

					   we also replace -b as follows:

					   a oriented away from, b oriented away from the corner:  b -> -b-a		
					   a oriented away from, b oriented towards the corner:    b -> a-b
					   a oriented towards, b oriented away from the corner:    b -> -ba
					   a oriented towwards b oriented towards the corner:      b -> -a-b

					   The folded 2-cell looses folded edge a and it's image in the 'old' b, just retaining 
					   the 'new' b at the folded corner.  All other instances of the 'old' b in this 2-cell are
					   replaced as above.

					   Note that we have to remember that the bottom word is specified from the far corner towards
					   the corner being folded but we are concerned her with the orientation of the folded edges with
					   respect to the corner 'at' the fold.
					*/ 

					bool top_away;
					char top_edge;

					if (lptr->first[0] == '-')
					{
						top_away = false;
						top_edge = lptr->first[1];
					}
					else
					{
						top_away = true;
						top_edge = lptr->first[0];
					}

					bool bottom_away;
					char bottom_edge;
					int len = lptr->second.length();

					if (len > 1 && lptr->second[len-2] == '-')
						bottom_away = true; // bottom word in reverse direction to the top
					else
						bottom_away = false;

					bottom_edge = lptr->second[len-1];														
					
if (decomp_control::DEBUG >= decomp_control::BASIC)		
{
	debug << "jsj::fold:   performing Type I fold, top_edge = " << top_edge << ", bottom_edge = "  << bottom_edge << endl; 
	debug << "jsj::fold:     top_edge oriented " << (top_away? "away from" : "towards") <<" the folded corner" << endl; 
	debug << "jsj::fold:     bottom_edge oriented " << (bottom_away? "away from" : "towards") <<" the folded corner" << endl; 
}

					 /* to perform the fold we take the top_edge and bottom edge out of w and w', make the 
						substitutions described above, then put back the long edge.  This avoids having to 
						check we only substitue 'old' b and not 'new' b.
					 */
					 if (top_away)
						 lptr->first = lptr->first.substr(1,string::npos);
					 else
						 lptr->first = lptr->first.substr(2,string::npos);

					 if (bottom_away)
						 lptr->second = lptr->second.substr(0,len-2);
					 else
						 lptr->second = lptr->second.substr(0,len-1);

if (decomp_control::DEBUG >= decomp_control::BASIC)		
	debug << "jsj::fold:     after removong top and bottom edge w = "  << lptr->first << ", w' = " << lptr->second << endl; 
	
					 /* substitute 'old' b for 'new' b */
					 char a,b;
					 bool a_away, b_away;

					 if (edge_length[top_edge] < edge_length[bottom_edge])
					 {						
						 a = top_edge;
						 a_away = top_away;

						 b = bottom_edge;
						 b_away = bottom_away;
					 }
					 else
					 {
						 a = bottom_edge;
						 a_away = bottom_away;

						 b = top_edge;
						 b_away = top_away;
					}
					cout << "\n" << num_folds << ". folding edge " << a << ", onto edge "  << b << endl; 
					output << "\n" << num_folds << ". folding edge " << a << ", onto edge "  << b << endl; 
							
//					cout << "  folding edge " << a << ", length " << edge_length[a]
//				         << ", onto edge "  << b << ", length " << edge_length[b] << endl; 
//					output << "  folding edge " << a << ", length " << edge_length[a]
//				         << ", onto edge "  << b << ", length " << edge_length[b] << endl; 

if (decomp_control::DEBUG >= decomp_control::BASIC)		
{
	debug << "jsj::fold:     folding edge " << a << ", length " << edge_length[a]
	      << ", onto edge "  << b << ", length " << edge_length[b] << endl; 
}							
					 list<two_cell>::iterator sptr = two_cells.begin(); // substitute pointer

					 while (sptr != two_cells.end())
					 {
						if (sptr->first == "folded" || sptr->first.find("quadratic") != string::npos)
						{
							sptr++;
							continue;
						}							
								
if (decomp_control::DEBUG >= decomp_control::BASIC)		
	debug << "jsj::fold:     substituting " << b << " in w = "  << sptr->first << endl; 

						substitute(sptr->first,a,b,a_away,b_away);

if (decomp_control::DEBUG >= decomp_control::BASIC)		
	debug << "jsj::fold:     substituting " << b << " in w' = "  << sptr->second << endl; 

						substitute(sptr->second,a,b,a_away,b_away);

						sptr++;
					};

					/* put back 'new' b */
					if (edge_length[top_edge] < edge_length[bottom_edge])
					{						

						if (bottom_away)
							lptr->second = lptr->second + "-" + string(1,bottom_edge);
						else
							lptr->second = lptr->second + string(1,bottom_edge);
if (decomp_control::DEBUG >= decomp_control::BASIC)		
	debug << "jsj::fold:     restored 'new' " << bottom_edge << " to bottom, w' = "  << lptr->second << endl; 

					}
					else
					{
						if (top_away)
							lptr->first = string(1,top_edge) + lptr->first;
						else
							lptr->first = "-" + string(1,top_edge) + lptr->first;							
								
if (decomp_control::DEBUG >= decomp_control::BASIC)		
	debug << "jsj::fold:     restored 'new' " << top_edge << " to top, w = "  << lptr->first << endl; 

					}
							
					/* adjust the length of b */
					edge_length[b] -= edge_length[a];

					cout << "\n  revised two_cell list:" << endl;
					output << "\n  revised two_cell list:" << endl;
					list<two_cell>::iterator optr = two_cells.begin(); //output pointer

					while (optr != two_cells.end())
					{
						if (optr->first != "folded")
						{
							if (optr->first == "quadratic pair")
							{
								cout << "  quadratic pair, attaching map = " << optr->second << endl;
								output << "  quadratic pair, attaching map = " << optr->second << endl;
							}
							else if (optr->first == "quadratic repeating word")
							{
								cout << "  quadratic repeating word, attaching map = " << optr->second << endl;
								output << "  quadratic repeating word, attaching map = " << optr->second << endl;
							}
							else
							{
								cout << "  w = " << optr->first << ", w' = " << optr->second << ", \\bar w' = " << invstr(optr->second)<< endl;
								output << "  w = " << optr->first << ", w' = " << optr->second << ", \\bar w' = " << invstr(optr->second)<< endl;
							}
						}

						optr++;
					};

					cout << "\n  revised edge lengths:\n";
					output << "\n  revised edge lengths:\n";
					map<char,surd>::iterator eptr(edge_length.begin());
					while (eptr != edge_length.end())
					{
						cout << "  " << eptr->first << " = " << eptr->second << " = " << eptr->second.estimate() << endl;
						output << "  " << eptr->first << " = " << eptr->second << " = " << eptr->second.estimate() << endl;
						eptr++;
					};


if (decomp_control::DEBUG >= decomp_control::BASIC)
{
	debug << "jsj::fold:     revised edge lengths: ";
	map<char,surd>::iterator eptr(edge_length.begin());
	while (eptr != edge_length.end())
	{
		debug << eptr->first << '=' << eptr->second << ' ';
		eptr++;
	};
	debug << endl;
}	

				}
			};
			
if (decomp_control::DEBUG >= decomp_control::SUMMARY)
	debug << "jsj::fold: finished folding 2-cell" << endl;

			lptr++;				
		};
	
		/*  If we have reached the end of the list of two-cells and carried out
		    the maximum number of folds on the last two-cell that required folding 
		    this indicates we have traversed the list and performed at most max_folds
		    on each 2-cell, so break.
		*/
		if (num_folds == max_folds)
			break;
		
	} while (!finished_folding);	
	
	cout << "\nfinal two_cell list:" << endl;
	output << "\nfinal two_cell list:" << endl;
	lptr = two_cells.begin();
	bool empty_list = true;
	while (lptr != two_cells.end())
	{
		if (lptr->first != "folded")
		{
			empty_list = false;
			if (lptr->first == "quadratic pair")
			{
				cout << "  quadratic pair, attaching map = " << lptr->second << endl;
				output << "  quadratic pair, attaching map = " << lptr->second << endl;
			}
			if (lptr->first == "quadratic repeating word")
			{
				cout << "  quadratic repeating word, attaching map = " << lptr->second << endl;
				output << "  quadratic repeating word, attaching map = " << lptr->second << endl;
			}
			else
			{
				cout << "  w = " << lptr->first << ", w' = " << lptr->second << ", \\bar w' = " << invstr(lptr->second)<< endl;
				output << "  w = " << lptr->first << ", w' = " << lptr->second << ", \\bar w' = " << invstr(lptr->second)<< endl;
			}
		}
		
		lptr++;
	};
	
	if (empty_list)
	{
		cout << "empty" << endl;
		output << "empty" << endl;
	}
	
	cout << endl;
	output << endl;

#endif	
}

#ifdef USE_FOLD
int word_length (string& word)
{
	int count=0;
	for (unsigned int i=0; i< word.length(); i++)
	{
		if (word[i] != '-')
			count++;
	}
	
	return count;
}

void dump(ostream& os, list<two_cell>& l)
{
	list<two_cell>::iterator lptr = l.begin();
	
	while (lptr != l.end())
	{
		os << "list<two_cell>::iterator = " << &(*lptr) << ": &first = " << &(lptr->first);
		os << " value = " << lptr->first << ": &second = " << &(lptr->second) << " value = " << lptr->second << endl;
		lptr++;
	};
	os << endl;

}

bool folding_required(two_cell& tc,map<char,surd>& edge_length)
{
	if (++num_folds == max_folds)
	{
		cout << "\nMaximum fold threshold (" << max_folds << ") reached" << endl;
		output << "\nMaximum fold threshold (" << max_folds << ") reached" << endl;

if (decomp_control::DEBUG >= decomp_control::BASIC)
	debug << "jsj::fold: maximum fold threshold (" << max_folds << ") reached" << endl;
		return false;
	}
	else if (quadratic_pair(tc))
	{
		/* mark the 2-cell as a quadratic pair */
		tc.first += tc.second;
		tc.second = tc.first;
		tc.first = "quadratic pair";
		
		return false;
	}
	else if (QUADRATIC_WORD)
	{
		/* we already know this to be a quadratic word so check the word_list to see if we have a repeat */
		bool found = false;
		


		/* mark the 2-cell as a quadratic pair */
		tc.first += tc.second;
		tc.second = tc.first;
		tc.first = "quadratic repeating word";
			
		return !found; // if we've found a repeat we don't need to continue folding
	}
	else if (quadratic_word(tc))
	{
		initialize_word_list(tc,edge_length);		
		return true;
	}
	else
		return true;
}

/* quadratic_pair returns if the two-cell is of the form w=ab w'=-a-b */
bool quadratic_pair(two_cell &tc)
{
	string::iterator tptr = tc.first.begin();
	char char_1;
	char char_2;
	bool char_1_negative = false;
	bool char_2_negative = false;
	
    /* first isolate the first two characters and signs in the top word. */	
	if (*tptr == '-')
	{
		char_1_negative = true;
		tptr++;
	}
	
	char_1 = *tptr;
	tptr++;

	if (*tptr == '-')
	{
		char_2_negative = true;
		tptr++;
	}
	
	char_2 = *tptr;
	tptr++;
	
	if (tptr != tc.first.end())
	{
		/* top word contains more than two characters, 
		   so not of form w=ab, return false		   
		*/
		return false;
	}
	
	/* now check that the bottom word is the same as the top
	   but with the two characters inverted.
	*/
	string::iterator bptr = tc.second.begin();
	if( *bptr == '-' )
	{
		bptr++;
		
		if (char_1_negative || *bptr != char_1)
			return false;
	}
	else if (!char_1_negative || *bptr != char_1)
			return false;

	bptr++; // move over first character

	if( *bptr == '-' )
	{
		bptr++;
		
		if (char_2_negative || *bptr != char_2)
			return false;
	}
	else if (!char_2_negative || *bptr != char_2)
			return false;
	
	bptr++; // move over second character

	if (bptr != tc.second.end())
	{
		/* bottom word contains more than two characters, 
		   so not of form w=-a-b, return false
		   
		*/
		return false;
	}
	else
		return true;
}

bool quadratic_word(two_cell &tc)
{
	return true;
}

void initialize_word_list(two_cell& tc,map<char,surd>& edge_length)
{
	QUADRATIC_WORD = true;
	
	word_list.clear();
	two_cell_record initial_record(tc,edge_length);
	word_list.push_back(initial_record);
}


/* The function substitute replaces (old)b in word with the result of a folded onto (old)b.
   The exact substitution depends on the orientation of a and (old)b, as given by a_away and 
   b_away, but results in expressing word in terms of a and (new)b
*/
void substitute(string& word, char a, char b, bool a_away, bool b_away)
{
	string positive_replacement;
	string negative_replacement;

	string as = string(1,a);
	string bs = string(1,b);

if (decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "jsj::substitute: word = " << word << ", a = " << a << ", b = " << b 
	      << ", a_away = " << (a_away? "true": "false") << ", b_away = " << (b_away? "true": "false") << endl;
}

	if (a_away)
	{
		if (b_away)
			positive_replacement = as + bs;
		else
			positive_replacement = bs + "-" + as;		
	}
	else
	{
		if (b_away)
			positive_replacement = "-" + as + bs;
		else
			positive_replacement = bs + as;
	}

	if (a_away)
	{
		if (b_away)
			negative_replacement = "-" + bs + "-" + as;
		else
			negative_replacement = as + "-" + bs;		
	}
	else
	{
		if (b_away)
			negative_replacement = "-" + bs + as;
		else
			negative_replacement = "-" + as + "-" + bs;
	}
	

if (decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "jsj::substitute:     positive_replacement = " << positive_replacement
	      << ", negative_replacement = " << negative_replacement << endl;
}	

	bool negative;
	
	/* in the following loop, p is incremented within the body of the loop as appropriate */
	for (string::size_type p=0;p<word.length();) 
	{
		if(word[p] == '-')
		{
			negative = true;
			p++;

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "jsj::substitute:     incremented p to " << p << " word[p] = " << word[p] << endl;
		}
		else
			negative = false;
		
		if (word[p] == b)
		{
			if (negative)
			{

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "jsj::substitute:     negative_replacement at position " << p << endl;

				/* remove the '-' since it's taken care of by the replacement */
				word.erase(p-1,1);
				p--;
				
				word.replace(p,1,negative_replacement);
				p += negative_replacement.length();

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "jsj::substitute:     word becomes " << word << ", new p = " << p << endl;
			}			
			else
			{

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "jsj::substitute:     positive_replacement at position " << p << endl;

				word.replace(p,1,positive_replacement);
				p += positive_replacement.length();

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "jsj::substitute:     word becomes " << word << ", new p = " << p << endl;
			}
		}
		else 
			p++;		
	}

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "jsj::substitute:     word change to " << word << endl;
}


/* invstr creates the inverse of src and returns it to the call. */
string invstr(string& str)
{
	string inverse;
	
    /* start at the end of the string and work back */
    string::reverse_iterator rptr = str.rbegin();
    
    while (rptr != str.rend())
    {
		char c = *rptr;
		rptr++;
		if (rptr != str.rend())
		{
			if (*rptr == '-')
				rptr++;
			else
				inverse += "-";
		}
		else
			inverse += "-";
		
		inverse += string(1,c);		
	}
	
    return inverse;
}
#endif
