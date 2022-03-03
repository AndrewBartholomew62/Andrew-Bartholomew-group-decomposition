using namespace std;

#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <list>
#include <cstring>


extern unsigned int DEBUG;

extern ofstream	debug;


#include <matrix.h>
#include <util.h>
#include <decomp.h>

//void print(matrix<int> m, ostream& s, int n, string prefix="");

/* The function inv takes a string s and inverts it, as an element of G.
   It assumes that there is sufficient memory at the location passed to 
   store the inverse word
*/
string G_inv(string s)
{
	string inverse;
    unsigned int s_length = s.length();

    if (s_length == 0)
		return inverse; /* inverse of null string is null */

	bool skip = false;
    for (unsigned int i = s_length - 1; i>0; i--)
    {
		if (skip)
	    	skip = false;
		else if ( s[i-1] == '-')
		{
	    	inverse += string(1,s[i]);
	    	skip = true;
		}
		else
	    	inverse += "-" + string(1,s[i]);
    }

    if (!skip)
    	inverse += "-" + string(1,s[0]);

	return inverse;	
}

void find(int label,int number_of_1_cells_in_rosette, matrix<int>** labels_for_1_cell, vector<int>& position)
{
    /* This procedure looks for 'label' in the t or b-label matrices 
    and writes its position to 'position' as (1-cell,row,column). */

    for (int i=0; i< number_of_1_cells_in_rosette; i++)
    {
		matrix<int>& labels = *labels_for_1_cell[i];

		for (unsigned int j = 0; j < labels.numrows(); j++)
		{
	    	for (unsigned int k = 0 ; k < labels.numcols(); k++)
			if ( labels[j][k] == label )
			{
		    	position[0] = i;
		    	position[1] = j;
		    	position[2] = k;
				return;
			}
		}
    }
}

int bvertex(int label, int number_of_0_cells, int numbcols, vector<int>** blabels_for_0_cell)
{
    /* This procedure evaluates the b-vertex in the re-labelled btrees
       to which the 0-cell b-label 'label' is attached. */
    int on_0_cell = numbcols + 1;
    bool not_found = true;

    for (int i=0; i < number_of_0_cells; i++)
    {
		for (unsigned int j = 0; j < blabels_for_0_cell[i]->size(); j++)
		{
	    	if ( (*blabels_for_0_cell[i])[j] == label )
	    	{
				not_found = false;
				break;
	    	}
		}

		if ( not_found )
	    	on_0_cell++;
		else
	    	break;
    }
    return on_0_cell;
}

/*********************************************************************
 This procedure gives the word in G corresponding to the arc l1-->l2 on
 some relator disc.  It deals with both t and b cases, and with arcs in
 trees or not.
 
	max_1_cell_labels will be either mbtl or m1cbl,  
	one_cell_labels will be either tlabels_for_1_cell or blabels_for_1_cell

 
*********************************************************************/
string disc_word
(
	int 				l1, 
	int 				l2, 
	bool 				t_labels, 
	bool 				in_a_tree,
    vector<int>& 		max_1_cell_labels, 
	matrix<int>** 		one_cell_labels,
	vector<string>&	flag_on_1_cell, 
	int* 				n_tpts_on_1_cell, 
	matrix<int>& 		ext0cbl, 
	int 				max1cbl,
	int 				number_of_1_cells_in_rosette, 
	vector<int> 				length_of_attaching_map, 
	matrix<int>& 		attaching_map
)
{
	string word;

    vector<int> minpos(3);
	vector<int> maxpos(3);
    /* "min(max)imum label position" will be used with 'find'. */

    int 	max;
    int 	min;
    if ( l1 < l2 )
    {
		max = l2;
		min = l1;
    }
    else
    {
		max = l1;
		min = l2;
    }

    if ( t_labels || max <= max1cbl )
    {
		/* Evaluate which disc the labels lie on. */
		int disc = -1;
		do disc++; while ( min > max_1_cell_labels[disc]);
		
		/* Determine which 1-cells the labels 'min' and 'max' lie on and the
		   place in 'map' corresponding to these 1-cells. */


		find(min, number_of_1_cells_in_rosette, one_cell_labels, minpos);
		find(max, number_of_1_cells_in_rosette, one_cell_labels, maxpos);

		/* 'min(max)pos[1]' now denotes the respective ocurrence we are
		   looking for.  We count how many occurrences of the 'minpos[0]'th
		   and 'maxpos[0]'th 1-cells there are in the discs 1,...,disc-1 and
		   then use 'min(max)pos[2]' to determine which occurrence in the attaching 
		   map carries 'min' and 'max'.
		*/

		int minocc = -1;
		int maxocc = -1;
		int minp = -1;
		int maxp = -1;
		
		for (int  i = 0; i < disc; i++) 
		{
	    	for (int j = 0; j < length_of_attaching_map[i]; j++)
	    	{
				if ( abs(attaching_map[i][j]) -1 == minpos[0])
					minocc++;

				if ( abs(attaching_map[i][j]) -1 == maxpos[0])
					maxocc++;
	    	}
		}

		/* First find the "place" of 'min' in 'map'. */
		for (int i = 0; i < length_of_attaching_map[disc]; i++)
		{
	    	if ( abs(attaching_map[disc][i]) - 1 == minpos[0] )
	    	{
				minocc++;
				minp++;
				if (minocc == minpos[1])
			    	break;
	    	}
	    	else
				minp++;
		}    

		/* Now 'max'. */
		for (int i = 0; i < length_of_attaching_map[disc]; i++)
		{
	    	if ( abs(attaching_map[disc][i]) - 1 == maxpos[0] )
	    	{
				maxocc++;
				maxp++;

				if ( maxocc == maxpos[1] )
			    	break;
	    	}
	    	else
				maxp++;
		}

		/* We may now determine the corresponding word between 'min' and
		   'max' by moving through the attaching map. */

		Slice_iter<int> map = attaching_map.row(disc);

		if ( minp != maxp ) /* Otherwise the word is trivial */
		{    
	    	if ( map[minp] > 0  && flag_on_1_cell[map[minp]-1] != "" )
				word += flag_on_1_cell[map[minp]-1];

	    	for (int i= minp+1 ; i< maxp ; i++)
	    	{
				if ( map[i] > 0  && flag_on_1_cell[map[i]-1] != "" )
		    		word += flag_on_1_cell[map[i]-1];
					
				if ( map[i] < 0  && flag_on_1_cell[abs(map[i])-1] != "" )
				{
		    		word += "-";
		    		word += flag_on_1_cell[abs(map[i])-1];
				}
	    	}

	    	if ( map[maxp] < 0 && flag_on_1_cell[abs(map[maxp])-1] != "" )
	    	{
	    		word += "-";
	    		word += flag_on_1_cell[abs(map[maxp])-1];
	    	}
		}

		if ( l1 == max )
	    	word = G_inv(word);
    }
    else if ( min <= max1cbl )
    {
		/* we are dealing with b-labels and the max is a 0-cell b-label */
		find(min, number_of_1_cells_in_rosette, one_cell_labels, minpos);

		if (minpos[2] == 2*(n_tpts_on_1_cell[minpos[0]])-1 && flag_on_1_cell[minpos[0]] != "" )
	    	word += flag_on_1_cell[minpos[0]];

		if ( l1 == max )
	    	word = G_inv(word);
    }
    else
    {
		/* If l1-->l2 lies in a b-tree (necessarily the first - since it must 
		   contain O ) then the arc is "carried", in R(L), by the 
		   first 1-cell corresponding to a 2-torsion generator, so there is 
		   no flag on this 1-cell.
		*/
		if (!in_a_tree)
		{
	    	/* Calculate which 1-cell l1-->l2 corresponds to. */
	    	int disc = 0;
	    	int check = 0;
	    	int ocbl=min-max1cbl; // "0-cell b-label"=l, where min is the lth 0-cell b-label.

	    	do
	    	{
				disc++;
				check += length_of_attaching_map[disc-1];
	    	} while (ocbl > check);
	    	disc--;

	    	/* 'disc' now gives the relator disc "carrying" 'min'.
	    	   The 1-cell can now be calculated. */

	    	if ( min == ext0cbl[disc][0] && max == ext0cbl[disc][1] )
	    	{
				int one_cell = attaching_map[disc][length_of_attaching_map[disc]-1];
				word += flag_on_1_cell[abs(one_cell)-1];

				if ( (min == l1 && one_cell > 0) || (min == l2 && one_cell < 0) )
		    		word = G_inv(word);  // Notice we are moving anticlockwise round the disc in this case.
	    	}
	    	else
	    	{
				int one_cell = attaching_map[disc][length_of_attaching_map[disc]-1-check+ocbl];
				word += flag_on_1_cell[abs(one_cell)-1];

				if ( (min == l1 && one_cell < 0) || (min == l2 && one_cell > 0) )
		    		word = G_inv(word); 
	    	}
		}
    }
    
	return word;
}


/* This procedure determines a path from the initial point ("ip"), which lies
   in some maximal tree, to the chosen base point in that tree.  It 
   returns the word corresponding to moving along that path by writing 
   it to the char* output, which is assumed to have sufficient memory
   allocated to it. 
   
   tree will be either ttree or btree
   word will be either tvw or bvw
   
*/
string base_point_path(int v, bool t, bool twisted, int tbp, int bbp1, int bbp2, int* tree, string* word)
{
	string w;
    do
    {
		if ( (t && v == tbp) || (!t && twisted && v == bbp1)   || (!t && !twisted && (v == bbp1 || v == bbp2) ) )
	    	break;

		w += word[v-1]; /* Add word to output. */
		v = tree[v-1];  /* Move "up" tree. */
    } while (true);
	return w;
}


/* This procedure evaluates the image of the generator (given by "initial"-->"terminal"), 
   of some fundamental group, in G. 

   max_1_cell_labels will be either mbtl or m1cbl,  
   labels_for_1_cell will be either tlabels_for_1_cell or blabels_for_1_cell
   tree will be either ttree or btree
   word will be either tvw or bvw
   
*/
string image
(
	int 				in, 
	int 				ter, 
	bool 				t, 
	int* 				tree, 
	string* 			word, 
	vector<int>& 		max_1_cell_labels, 
	matrix<int>** 		labels_for_1_cell,
	bool 				twisted, 
	int 				tbp, 
	int 				bbp1, 
	int 				bbp2, 
	int* 				tccc, 
	int* 				bccc, 
	int 				number_of_0_cells, 
	int 				numbcols, 
	vector<int>** 		blabels_for_0_cell,            
	vector<string>&	flag_on_1_cell, 
	int* 				n_tpts_on_1_cell, 
	matrix<int>& 		ext0cbl, 
	int 				max1cbl,
	int 				number_of_1_cells_in_rosette, 
	vector<int> 				length_of_attaching_map, 
	matrix<int>& 		attaching_map
)
{
    /* This procedure evaluates the image of the generator (given by
    "initial"-->"terminal"), of some fundamental group, in G. */

    vector<int> position(3);
    int ini_v,ter_v; /* "initial/terminal vertex" */

    /* Evaluate the vertices to which 'in' and 'ter' are attached. */
    if (t)
    {
		find(in,number_of_1_cells_in_rosette, labels_for_1_cell, position);
		ini_v = tccc[position[0]] + position[2] + 1;
		find(ter,number_of_1_cells_in_rosette, labels_for_1_cell, position);
		ter_v = tccc[position[0]] + position[2] +1;
    }
    else
    {
		if (in <=max1cbl)
		{
			find(in,number_of_1_cells_in_rosette, labels_for_1_cell, position);
		    ini_v = bccc[position[0]] + position[2] + 1;
		}
		else
		    ini_v = bvertex(in, number_of_0_cells, numbcols, blabels_for_0_cell);
	
		if (ter <=max1cbl)
		{
			find(ter,number_of_1_cells_in_rosette, labels_for_1_cell, position);
		    ter_v = bccc[position[0]] + position[2] + 1;
		}
		else
		    ter_v = bvertex(ter, number_of_0_cells, numbcols, blabels_for_0_cell);
    }
	
    /* The required word is now easy to determine */
    string w = base_point_path(ini_v, t, twisted, tbp, bbp1, bbp2, tree, word);
    w = G_inv(w);
			 	 
    w += disc_word(in, ter, t, false, max_1_cell_labels, labels_for_1_cell, flag_on_1_cell, n_tpts_on_1_cell, 
	               ext0cbl, max1cbl, number_of_1_cells_in_rosette, length_of_attaching_map, attaching_map);
				   
    w+= base_point_path(ter_v, t, twisted, tbp, bbp1, bbp2, tree, word);
	return w;
}

vector<int>* find_row (int label, int number_of_0_cells, vector<int>** blabels_for_0_cell)
{    
    /* This function isolates the labels on a 0-cell,
       containing the given label. */

    for (int i =0; i< number_of_0_cells; i++)
    {
		vector<int>& labels = *blabels_for_0_cell[i];
	
		for (unsigned int j=0; j < labels.size(); j++)
		{
		    if (labels[j] == label)
				return blabels_for_0_cell[i];
		}
    }
	return 0; // to stop control reaching end of non-void function warning
}


void set_negative
(
	int 			l1, 
	int 			l2, 
	matrix<int>& 	jbba,
	int 			max1cbl, 
	int 			number_of_1_cells_in_rosette, 
	matrix<int>**	blabels_for_1_cell, 
	matrix<int>**	b0cbl_for_1_cell, 
	int* 			occurrences_of_1_cell, 
	vector<int> 			length_of_attaching_map,	
	matrix<int>& 	attaching_map,
	matrix<int>& 	ext0cbl
)
{
    /**********************************************************************
    This procedure sets all occurrences of an arc identified with l1-->l2 to
    be negative in 'jbba'.  There are only three possibilities we shall need
    to consider, by the way we removed repetitions.
    ***********************************************************************/

	vector<int> position(3);
    if ( l1 <= max1cbl && l2 <= max1cbl )
    {
		/* The two labels must be adjacent in some 1-cell 
		   on the boundary of some relator disc. */
		
		find (l1,number_of_1_cells_in_rosette,blabels_for_1_cell,position);
		int one_cell_1 = position[0];
		int col1 = position[2];
		matrix<int>& matptr1 = *blabels_for_1_cell[one_cell_1];

		find (l2,number_of_1_cells_in_rosette,blabels_for_1_cell,position);
		int one_cell_2 = position[0];
		int col2 = position[2];
		matrix<int>& matptr2 = *blabels_for_1_cell[one_cell_2];

		for (int i = 0; i < occurrences_of_1_cell[one_cell_1]; i++ )
		{
	    	int p1 = matptr1[i][col1];
	    	int p2 = matptr2[i][col2];

	    	if ( p1 < p2 )
			{
				if (jbba[p1-1][0] == p2)
			    	jbba[p1-1][0] *= -1;
				else
			    	jbba[p1-1][2] *= -1;
			}
	    	else
			{
				if (jbba[p2-1][0] == p1)
			    	jbba[p2-1][0] *= -1;
				else
			    	jbba[p2-1][2] *= -1;
			}
		}
    }
    else if ( l1 <= max1cbl && l2 > max1cbl )
    {
		/* All occurrences of this arc will be recorded in the rows of 'jbba'
		   determined by the column of the b-label matrices containing l1.  
		   An occurrence will be recorded in each of these rows, either in 
		   column 1 or 3 (but not both, since the relators are reduced words).
		   The column not denoting the occurrence will contain zero.
		*/
		find (l1,number_of_1_cells_in_rosette,blabels_for_1_cell,position);
		int one_cell_1 = position[0];
		int col1 = position[2];
		matrix<int>& matptr1 = *blabels_for_1_cell[one_cell_1];

		for (int i = 0; i < occurrences_of_1_cell[one_cell_1]; i++ )
		{
	    	int p1 = matptr1[i][col1];

			if (jbba[p1-1][0] != 0)
		    	jbba[p1-1][0] *= -1;
			else
		    	jbba[p1-1][2] *= -1;
		}
    }
    else
    {
		/* Both l1 and l2 > max1cbl.  Notice that l1<l2 in this case, by the 
		   way we removed repetitions.  First we find the 1-cell in the 
		   rosette carrying the arc l1-->l2.
		*/
		int disc = 0;
		int check = 0;
		int ocbl = l1-max1cbl;
		 /* "0-cell b-label"=l, where l1 is the lth 0-cell b-label. */

		do
		{
	    	check += length_of_attaching_map[disc++];
		} while ( ocbl>check );
		disc--;

		/* 'disc' now gives the relator disc "carrying" l1.  The 1-cell can
		   now be calculated. */

		int one_cell_1;
		if ( l1 == ext0cbl[disc][0] && l2 == ext0cbl[disc][1] )
	    	one_cell_1 = abs(attaching_map[disc][length_of_attaching_map[disc]-1])-1;
		else
	    	one_cell_1 = abs(attaching_map[disc][length_of_attaching_map[disc]-1-check+ocbl])-1;

		/* Now use 'b0cbl for 1 cell' to negate the required occurrences. */
		matrix<int>& matptr1 = *b0cbl_for_1_cell[one_cell_1];
		
		for (unsigned int i = 0; i < matptr1.numrows(); i++ )
		{
	    	int p1 = matptr1[i][0];

			if (jbba[p1-1][0] == matptr1[i][1])
		    	jbba[p1-1][0] *= -1;
			else
		    	jbba[p1-1][2] *= -1;
		}
    }
}

/* find_next occurrence of s2 in s1, starting from p */
string::size_type find_next( string& s1, string::size_type p, string& s2)
{
	string ss = s1.substr(p);
	string::size_type pp = ss.find(s2);
	
	if (pp == string::npos)
		return pp;
	else
		return pp+p;
}


string simple_reduce(string str)
{
    /* This procedure removes adjacent cancelling pairs of generators in the
      string s.*/

	int length = str.length();
	
    if (length < 3)
		return str;

	char* s = c_string(str);
    int place = 0;
    bool skip_next_char = false;
    bool cancelling_pair_found = false;
	
    do
    {    
		if ( skip_next_char )
	    	skip_next_char = false;
		else if ( *(s+place) == '-' )
		{
	    	if ( *(s+place+1) == *(s+place+2) )
				cancelling_pair_found = true;
	    	else
				skip_next_char = true;
		}
		else
		{
	    	if ( *(s+place+1) == '-' && *(s+place+2) == *(s+place))
				cancelling_pair_found = true;
		}

		if ( !cancelling_pair_found )
	    	place++;

    } while (!cancelling_pair_found && place <= length-3);

    if (cancelling_pair_found)
    {
		*(s+place) = '\0';
		strcat ( s, s+place+3);
		return simple_reduce(s);
    }

	string result = s;	
	delete [] s;
	return result;
}

void reduce(string* generators, int first_gen, int num_gens, vector<string>& relator, int number_of_relators)
{
    bool reduced_generators = false;

if (decomp_control::REDUCE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "reduce: reduce called with first_gen = " << first_gen << " num_gens = " << num_gens << endl;		

	/* first, simple reduce generators to remove adjacent cancelling pairs */
    for (int i=0; i< num_gens; i++)
		generators[first_gen+i] = simple_reduce(generators[first_gen+i]);


    /* now remove duplicates and inverse duplicates */
    for (int i=0; i < num_gens-1; i++)
    {
		string str1 = generators[first_gen+i];
		if (str1.length())
		{
	    	for (int j=i+1; j < num_gens; j++)
	    	{
				string str2 = generators[first_gen+j];
				if (str2.length())
				{
				    if (str1 == str2)
					{
if (decomp_control::REDUCE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "reduce: duplicate found: i=" << i+1 << " j=" << j+1 << endl;		
						generators[first_gen+j] = "";
					}
				    else
				    {
						string inv = G_inv(str2);
						generators[first_gen+j] = inv;
						if (str1 == inv)
						{
if (decomp_control::REDUCE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "reduce: duplicate inverse found: i=" << i+1 << " j=" << j+1 << endl;		
					    	generators[first_gen+j] = "";
						}
				    }
				}
	    	}
        }
    }

if (decomp_control::DECOMPOSE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "reduce: generators after initial reduction" << endl;
	for (int i=0; i<num_gens; i++)
	{
	    debug << "reduce:   generator " << i+1 << " ";
	    
	    if (generators[first_gen+i] == "")
			debug << "<>";
		else
			debug << generators[first_gen+i];
			
		debug << endl;
	}
}

    /******************************************
     Remove the remainder of this function to 
     retrieve original 1987 output format.			    
     *****************************************/

    /******************************************************************
     Reduce the generators as follows: for each generator g(i), check
     all generators g(j), j != i, to see whether g(j) or inv(g(j)) appears
	 as a subword of g(i). If so which case the subword can be removed to 
	 simplify the generator g(i).  In fact we check for subwords in both 
     g(i) and inv(g(i)).
     *****************************************************************/

    for (int i=0; i<num_gens; i++)
    {
		string str_i = generators[first_gen+i];
    	string::size_type li = str_i.length();

if (decomp_control::REDUCE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "reduce: check for presence of generator " << i+1 << " = ";
	if (li == 0)
		debug << "<>";
	else
		debug << str_i;
	debug << endl;
}

		if (li != 0)
		{
			for (int j=0; j <num_gens; j++)
			{
				if (j==i)
					continue;
					
				string str_j = generators[first_gen+j];
	    		string::size_type lj = str_j.length();
	    		if (lj != 0)
	    		{
					if (li < lj)
					{
		    			/* check for str_i at the end of str_j */
						if (str_j.rfind(str_i) == lj-li && str_j[lj-li-1] != '-' )
						{
if (decomp_control::REDUCE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "reduce:   found at end of generator " << j+1 << " = " << str_j << endl;
							str_j.erase(lj-li, li);
							lj -= li;
							generators[first_gen+j] = str_j;
							reduced_generators = true;
						}						

						if ( li < lj ) /* lj may have changed */
						{
			    			/* check for str_i at the beginning of str_j */
							if (str_j.find(str_i) == 0)
							{
if (decomp_control::REDUCE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "reduce:   found at start of generator " << j+1 << " = " << str_j << endl;
								str_j.erase(0, li);
								generators[first_gen+j] = str_j;
								lj -= li;
								reduced_generators = true;						
							}
						}
					}

					str_i = G_inv(str_i);
					generators[first_gen+i] = str_i;
					li = str_i.size();

					if ( li < lj )
					{
		    			/* check for inv(str_i) at the end of str_j */		
						if (str_j.rfind(str_i) == lj-li && str_j[lj-li-1] != '-' )
						{
if (decomp_control::REDUCE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "reduce:   inverse = " << str_i << " found at end of generator " << j+1 << " = " << str_j << endl;
							str_j.erase(lj-li, li);
							generators[first_gen+j] = str_j;
							lj -= li;
							reduced_generators = true;
						}						


						if ( li < lj ) /* lj may have changed */
						{
							/* check for inv(str_i) at the beginning of str_j */
							if (str_j.find(str_i) == 0)
							{
if (decomp_control::REDUCE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "reduce:   inverse = " << str_i << " found at start of generator " << j+1 << " = " << str_j << endl;
								str_j.erase(0, li);
								generators[first_gen+j] = str_j;
								lj -= li;
								reduced_generators = true;						
							}
						}
					}
	    		}
			}
		}
    }

    /*************************************************************
     Check to see whether any of the relators or their inverses are 
     contained as substrings in any of the generators.  If so, remove 
     them, since the substrings are equal to the identity.
     *************************************************************/
    for (int i=0; i<num_gens; i++)
    {
		string str_i = generators[first_gen+i];
		if ( str_i.length() != 0 )
		{

if (decomp_control::REDUCE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "reduce: looking for relators in generator " << i+1 << " = " << str_i << endl;
	
	    	for (int j=0; j < number_of_relators; j++)
	    	{
				string::size_type pj = 0;			
				string rel_j = relator[j];
if (decomp_control::REDUCE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "reduce:   checking relator " << rel_j << ", pj = " << pj << endl;

				do
				{
					pj = find_next(str_i, pj, rel_j);
if (decomp_control::REDUCE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "reduce:     find_next returned ";
	if (pj == string::npos)
		debug << "string::npos";
	else
		debug << pj;
	debug << endl;
}					
		    		if (pj != string::npos)
		    		{
						if (pj == 0 || str_i[pj-1] != '-')
						{
if (decomp_control::REDUCE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "reduce: relator " << rel_j << " found in generator " << i+1 << ": " << str_i << endl;
							str_i.erase(pj,rel_j.length());
if (decomp_control::REDUCE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "reduce:   reduced to " << str_i << endl;
							generators[first_gen+i] = str_i;
				    		reduced_generators = true;
						}
						else
						{
							pj += rel_j.length();
							if (pj == str_i.length())
								break;
						}
		    		}
					
				} while(pj != string::npos);

				/* now check the inverse of the relator */
				string rel_inv = G_inv(rel_j);

				pj = 0;
				
if (decomp_control::REDUCE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "reduce:   checking relator inverse " << rel_inv << ", pj = " << pj << endl;

				do
				{
					pj = find_next(str_i, pj, rel_inv);

if (decomp_control::REDUCE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "reduce:     find_next returned ";
	if (pj == string::npos)
		debug << "string::npos";
	else
		debug << pj;
	debug << endl;
}					
					
					
		    		if (pj != string::npos)
		    		{
						if (pj == 0 || str_i[pj-1] != '-')
						{
if (decomp_control::REDUCE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "reduce: relator inverse " << rel_inv << " found in generator " << i+1 << ": " << str_i << endl;
							str_i.erase(pj,rel_inv.length());
if (decomp_control::REDUCE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "reduce:   reduced to " << str_i << endl;
							generators[first_gen+i] = str_i;
				    		reduced_generators = true;
						}
						else
						{
							pj += rel_inv.length();
							if (pj == str_i.length())
								break;
						}
		    		}
					
				} while(pj != string::npos);	
	    	}
		}
    }

    /*	We now ensure that generators currently described as '-a' are 
	replaced with the corresponding 'a'
    */
    for (int i=0; i<num_gens; i++)
    {
		string& str = generators[first_gen+i];
		if ( str.size() == 2 && str[0] == '-' )
		{
	    	/* invert this generator */
	    	str = G_inv(str);
		}
    }

    /* Now use recursion to check for cancelling pairs again. */
    if (reduced_generators)
	{
if (decomp_control::REDUCE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "reduce: recursing " << endl;
		reduce(generators, first_gen, num_gens, relator, number_of_relators);   
	}

if (decomp_control::REDUCE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "reduce: complete" << endl;
}


/*************************************************************************
This procedure decomposes the pattern P (:=P('pattern')) into its connected
components and returns the N-tuples corresponding to each component.  The
strategy for the decomposition is as follows.  For each relator disc we
label the points of intersection of P with the 1-skeleton of the disc with
consecutive integers (without repetitions) starting at O and labelling
clockwise around the boundary, then for each interior 1-cell ("clockwise")
starting at O.  We then use these labels to build up equivalence classes of
labels on each disc for the relator lbl1~lbl2 iff lbl1 and lbl2 are joined
by an edge of P.  This enables us to discover how the pattern crosses each
disc. We then build classes of "end" labels of the arcs across each disc by
considering how the discs are identified under the attaching maps (we shall
have one class for each component of P), finally returning to the crossing
arcs to determine the tuples corresponding to each component.

If short_split is true, split returns true if the pattern is a track and 
false otherwise.
************************************************************************/
bool split (vector<int>& pattern, K_complex K, bool short_split, ofstream& txt, bool return_components, matrix<int>** component_ptr)     
{
if (decomp_control::SPLIT_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "split: pattern provided: ";
	for (unsigned int i=0; i<pattern.size(); i++)
		debug << pattern[i] << ' ';
	debug << endl;
}

    /**********  read variables required by split from the K_complex **********/
	int number_of_1_cells = K.number_of_1_cells;
	int number_of_2_cells = K.number_of_2_cells;
	int number_of_relators = K.number_of_relators;
	int number_of_1_cells_in_rosette = K.number_of_1_cells_in_rosette;
	int* occurrences_of_1_cell = K.occurrences_of_1_cell;
	vector<int> length_of_attaching_map = K.length_of_attaching_map;
	matrix<int>& two_cell = *K.two_cell;
	matrix<int>& attaching_map = *K.attaching_map;
	


    /***********************************************************************
    For the 1-cells of the rosette, the total number of occurrences of 1-cell
    i in 'two cell' is the number of occurrences of i in the boundary of some
    relator disc.  We therefore set out, for each such 1-cell, a matrix with
    r rows and c columns, where r is the number of occurrences of the 1-cell
    in the boundary of the relator discs and c is the number of points of P
    intersecting that 1-cell.  For the interior 1-cells, which occur exactly
    once in a unique relator disc, we set out a 1xc matrix, where c is
    as above.  We store in these matrices the labels attached to the 1-cells
    corresponding to their rows (in accordance with the canonical order
    through the relator discs, ignoring identifications) in such a way that
    the columns of each matrix corresponding to a 1-cell in the rosette
    determines how those labels become identified in L.
    ***********************************************************************/

    matrix<int>* labels_for_1_cell[number_of_1_cells];
    /* labels_for_1_cell[i] will reference the "label matrix" 
       for the ith 1-cell.
    */

    int n_pts_on_1_cell[number_of_1_cells];
    /* "number of points on 1-cell" will be used later. */

    /*  First deal with the 1-cells in the rosette. */
    for (int i = 0; i< number_of_1_cells_in_rosette; i++)
    {
		bool one_cell_found = false;

		for (int j =0; j< number_of_2_cells; j++)
		{       
	    	for ( int k = 0; k<3; k++)
	    	{
				if ( abs(two_cell[j][k]) == i+1)
				{
		    		if (k == 2)	    
						n_pts_on_1_cell[i] = pattern[3*j+k] + pattern[3*j];
		    		else 
						n_pts_on_1_cell[i] = pattern[3*j+k] + pattern[3*j+k+1];

		    		one_cell_found = true;
		    		break;
				}  	     
	    	} 
			
	    	if (one_cell_found)
			break;
		}

		labels_for_1_cell[i] = new matrix<int> (occurrences_of_1_cell[i],n_pts_on_1_cell[i]);	      
    }

    /* Now deal with the interior 1-cells, for which purpose we shall need a
    different marker for 'labels for 1 cell'.  We need to know the first
    variable on each relator disc, so evaluate these first. */

    int start[number_of_relators];
    start[0] = 0;   

    for (int i = 1; i< number_of_relators; i++)
		start[i] = start[i-1] + 3 * (length_of_attaching_map[i-1] - 2);  

    int interior_1_cell = number_of_1_cells_in_rosette;
    for (int i =0; i< number_of_relators; i++)
    {
		int s = start[i];
		for (int j = s + 3 ; j < s + 3 * (length_of_attaching_map[i] - 2) - 2 ; j += 3 )
		{
	    	n_pts_on_1_cell[interior_1_cell] = pattern[j] + pattern[j+1];

		   labels_for_1_cell[interior_1_cell] = new matrix<int> (1,n_pts_on_1_cell[interior_1_cell]);   
		   interior_1_cell ++;
		}
    }

if (decomp_control::SPLIT_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "split: number of points on 1-cells: ";
	for (int i=0; i<number_of_1_cells; i++)
		debug << n_pts_on_1_cell[i] << ' ';
	debug << endl;
}

    /********************************************************************
    Assign the labels to the space just allocated.  We take the convention
    that, for the 1-cells in the rosette, the orientation of the labels in
    the label matrix is left to right if the 1-cell is "positive in the
    attaching map" and right to left otherwise.  This ensures the
    identification matching mentioned above.  We also record maximum and
    minimum label on each disc and the maximum label on the boundary of
    each disc, all for future use.
    ********************************************************************/
    int marker[number_of_1_cells_in_rosette];
    /* These markers will mark the rows in the label matrices. */

    for (int i=0; i< number_of_1_cells_in_rosette; i++)
		marker[i] = 0;

    int label = 1;
    /* 'one cell' will denote a 1-cell of the rosette; label will 
      be used to allocate labels. */

    interior_1_cell = number_of_1_cells_in_rosette;
	
    int minlb[number_of_relators];
    int maxlb[number_of_relators];   
    int mbl[number_of_relators];
    /*"minimum label","maximum label","maximum boundary label"*/

    for ( int i = 0; i < number_of_relators ; i++)
    {

		/* First assign to 'minlb[i]'. */
		minlb[i] = label;
		for (int j = 0; j < length_of_attaching_map[i]; j++)
		{
	    	int one_cell = attaching_map[i][j];
	    	if ( one_cell > 0 )
	    	{
				Slice_iter<int> eltptr = labels_for_1_cell[one_cell-1]->row(marker[one_cell-1]);
								 
				for (int k=0; k< n_pts_on_1_cell[one_cell-1]; k++)
		    		*eltptr++ = label++;

				marker[one_cell-1]++;
	    	}
	    	else
	    	{
				one_cell *= -1;
				
				Slice_iter<int> eltptr = labels_for_1_cell[one_cell-1]->row(marker[one_cell-1]);

				/* new code to allocate labels in reverse order 18/2/06 */
				label += n_pts_on_1_cell[one_cell-1]-1;
				
				for (int k=0; k< n_pts_on_1_cell[one_cell-1]; k++)
	    			*eltptr++ = label--;
				
				label += n_pts_on_1_cell[one_cell-1]+1;
				marker[one_cell-1]++;
	    	}
		}

		/* Now assign the maximum label on the boundary of
		 this disc to 'mbl[i]'.*/
		mbl[i] = label-1;

		for (int j=0; j< length_of_attaching_map[i]-3; j++)
		{
			Slice_iter<int> eltptr = labels_for_1_cell[interior_1_cell]->row(0);

	    	for( int k =0; k< n_pts_on_1_cell[interior_1_cell]; k++)
			*eltptr++ = label++;
			
	    	interior_1_cell++;
		}

		/* Now assign to 'maxlb[i]'. */
		maxlb[i] = label-1;

    }

    if (decomp_control::SPLIT_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
    {
		debug << "split: minlb = ";
		for (int i=0; i< number_of_relators; i++)
			debug << setw(4) <<minlb[i] << ' ';
		debug << endl;
			
		debug << "split: maxlb = ";
		for (int i=0; i< number_of_relators; i++)
			debug << setw(4) <<maxlb[i] << ' ';
		debug << endl;

		debug << "split: mbl = ";
		for (int i=0; i< number_of_relators; i++)
			debug << setw(4) <<mbl[i] << ' ';
		debug << endl;

		debug << "split: label matracies\n";
		for (int i = 0; i< number_of_1_cells; i++)
		{
	    	debug << "split: one cell " << i+1 << endl;
			if (n_pts_on_1_cell[i])
				print(*labels_for_1_cell[i], debug, 0, "split:   ");
			else
				debug << "split:   void intersection with track" << endl;
		}

    }

    /*****************************************************************
    We next determine which labels are joined by an edge of P on each
    relator disc; we store the information as pairs, forming the rows of a
    matrix.  We have one pair for each edge of P.
                                                    
         Given a variable and the two sets of labels the edges of P
    determined by that variable "join", we make the pairing by moving "into
    the vertex", ie if 'pattern[variable]'=x we must have x labels in each
    list paired by the corresponding edges of t, thus the pairings must be of
    the x labels closest to that vertex on the relator disc.

         Each variable lies in a two cell of one of three types, A, B, C.
    Type A 2-cells are those that appear first in a relator disc (with the
    canonical ordering), type B 2-cells are the remaining 2-cells of that
    disc, except the last, and type C 2-cells are those that appear last in a
    relator disc.  Within any 2-cell a variable appears in place 1, 2, or 3,
    determined, again, by the canonical order.
   
    We store the pairs corresponding to each relator disc separately, and
    with each pair shall store the corresponding variable,
    to assist later calculations.
    ********************************************************************/

    matrix<int>* pair[number_of_relators];
    int partial_sum[number_of_relators];

    /* 'partial sum' will indicate how many 
       pairs arise from each relator disc. */

    for (int i = 0 ; i < number_of_relators ; i++)
		partial_sum[i] = 0;
   
    /* Assign partial sums. */
    int variable = 0;

    for ( int i = 0; i < number_of_relators ; i++)
	for (int j = 0; j < 3*(length_of_attaching_map[i]-2) ; j++)
	    partial_sum[i] += pattern[variable++];
   
    /* Declare space for 'pair'. */
    for (int i=0; i<number_of_relators; i++)
		pair[i] = new matrix<int>(partial_sum[i], 3);   

    /****************************************************************
    We traverse the relator discs in the canonical order, thereby moving
    through 'pattern'.  We make use of 'marker' to determine which set of
    labels we should consider for the boundary 1-cells and 'interior 1 cell'
    to take care of the rest.
    ***************************************************************/

    for (int i=0; i< number_of_1_cells_in_rosette; i++)
		marker[i] = 0;

    interior_1_cell = number_of_1_cells_in_rosette;
    variable = 0; /* Reset 'variable' for assignation to 'pair'.*/
       
    for ( int i = 0; i< number_of_relators; i++)
    {
		int pair_marker = 0; /* reset pair_marker for next relator disk */

		if (length_of_attaching_map[i] == 3)
		{
	    	/* The relator must be of the form 'aaa' for some generator. */
	    	int one_cell = attaching_map[i][0]; 

	    	int mark = marker[one_cell-1];
	    	int pv = pattern[variable];      

	    	Slice_iter<int> r1 = labels_for_1_cell[one_cell-1]->row(mark);
	    	Slice_iter<int> r2 = labels_for_1_cell[one_cell-1]->row(mark+1);
	    	Slice_iter<int> r3 = labels_for_1_cell[one_cell-1]->row(mark+2);

	    	matrix<int>& matrixptr = *pair[i]; /* matrixptr now points at the target
						 of the asignation */
	    	/* Place 1. */
	    	r1 += pv-1;
	    	r3 += pv;

	    	for (int j=0; j< pv; j++)
			{
				matrixptr[pair_marker][0] = *r1--;
				matrixptr[pair_marker][1] = *r3++;
				matrixptr[pair_marker][2] = variable;
				pair_marker++;
			}

            variable++;

	    	/* Place 2. */
	    	r1 += pv+1;
	    	r2 += pv-1;

	    	for (int j=0; j< pv; j++)
			{
				matrixptr[pair_marker][0] = *r1++;
				matrixptr[pair_marker][1] = *r2--;
				matrixptr[pair_marker][2] = variable;
				pair_marker++;
			}

	    	variable++;

	    	/* Place 3. */
	    	r2 += pv+1;
	    	r3 -= pv+1; /* i.e. the pv th label in the row */

	    	for (int j=0; j< pv; j++)
			{
				matrixptr[pair_marker][0] = *r2++;
				matrixptr[pair_marker][1] = *r3--;
				matrixptr[pair_marker][2] = variable;
				pair_marker++;
			}

	    	variable++;

	    	marker[one_cell-1] += 3;
		}
		else
		{
	    	matrix<int>& matrixptr = *pair[i]; /* matrixptr now points at the target
						 of the asignation */
						 
	    	/* Type A 2-cell. */
	    	int one_cell1 = attaching_map[i][0];
	    	int one_cell2 = attaching_map[i][1];
			bool s1;
			bool s2;

	    	if (one_cell1 < 0 )
	    	{
				s1 = false;
				one_cell1 *= -1;
	    	}
	    	else
				s1 = true;

	    	if (one_cell2 < 0 )
	    	{
				s2 = false;
				one_cell2 *= -1;
	    	}
	    	else
				s2 = true;

	    	Slice_iter<int> r1 = labels_for_1_cell[one_cell1-1]->row(marker[one_cell1-1]);

			Slice_iter<int> r2 = r1; // assignation just for the declaration, since there's no default constructor 
			
	    	if (one_cell1 == one_cell2)
				r2 = labels_for_1_cell[one_cell1-1]->row(marker[one_cell1-1]+1); 						
	    	else
				r2 = labels_for_1_cell[one_cell2-1]->row(marker[one_cell2-1]); 

	    	Slice_iter<int> r3 = labels_for_1_cell[interior_1_cell]->row(0); 

	    	int n1 = labels_for_1_cell[one_cell1-1]->numcols();
	    	int n2 = labels_for_1_cell[one_cell2-1]->numcols();
	    	int n3 = labels_for_1_cell[interior_1_cell]->numcols();    

	    	/* Place 1. */
	    	int pv = pattern[variable];      

	    	if ( s1 ) r1 += pv-1; else r1 += n1-pv;
	    	r3 += pv-1;

	    	for (int j=0; j< pv; j++)
	    	{

				matrixptr[pair_marker][0] = *r1;
				matrixptr[pair_marker][1] = *r3--;
				matrixptr[pair_marker][2] = variable;
				pair_marker++;
				
				if ( s1 ) r1--; else r1++;
	    	}
	    	variable++;

	    	/* Place 2. */
	    	pv = pattern[variable];      

	    	if ( s1 ) r1 += n1-pv+1; else r1 -= n1-pv+1;
	    	if ( s2 ) r2 += pv-1; else r2 += n2-pv;

	    	for (int j=0; j< pv; j++)
	    	{
				matrixptr[pair_marker][0] = *r1;
				matrixptr[pair_marker][1] = *r2;
				matrixptr[pair_marker][2] = variable;
				pair_marker++;

				if ( s1 ) r1++; else r1--;
				if ( s2 ) r2--; else r2++;
	    	}
	    	variable++;

	    	/* Place 3. */
	    	pv = pattern[variable];      

	    	if ( s2 ) r2 += n2-pv+1; else r2 -= n2-pv+1;
	    	r3 += n3-pv+1; 

	    	for (int j=0; j< pv; j++)
	    	{
				matrixptr[pair_marker][0] = *r2;
				matrixptr[pair_marker][1] = *r3++;
				matrixptr[pair_marker][2] = variable;
				pair_marker++;

				if ( s2 ) r2++; else r2--;
	    	}
	    	variable++;

	    	marker[one_cell1-1] += 1;
	    	marker[one_cell2-1] += 1;
	    	/* Notice that if one_cell1=one_cell2 then this code is
	    	 still correct. */

	    	/* Type B 2-cells. */
	    	int a_map_place = 2;
	    	for ( int j =0; j< length_of_attaching_map[i]-4; j++)
	    	{
				int in_rosette = attaching_map[i][a_map_place];

				if (in_rosette < 0 )
				{
		    		s2 = false;
		    		in_rosette *= -1;
				}
				else
		    		s2 = true;

				Slice_iter<int> r1 = labels_for_1_cell[interior_1_cell]->row(0); 
				Slice_iter<int> r2 = labels_for_1_cell[in_rosette-1]->row(marker[in_rosette-1]);
				Slice_iter<int> r3 = labels_for_1_cell[interior_1_cell+1]->row(0); 
				
				int n1 = labels_for_1_cell[interior_1_cell]->numcols();
				int n2 = labels_for_1_cell[in_rosette-1]->numcols();
				int n3 = labels_for_1_cell[interior_1_cell+1]->numcols();

				/* Place 1. */
				pv = pattern[variable];      
				r1 += pv-1;
				r3 += pv-1;
				
				for (int k=0; k< pv; k++)
				{
					matrixptr[pair_marker][0] = *r1--;
					matrixptr[pair_marker][1] = *r3--;
					matrixptr[pair_marker][2] = variable;
					pair_marker++;
				}
					
				variable++;

				/* Place 2. */
				pv = pattern[variable];      
				r1 += n1-pv+1;
				if ( s2 ) r2 += pv-1; else r2 += n2-pv;

				for (int k=0; k< pv; k++)		
				{
					matrixptr[pair_marker][0] = *r1++;
					matrixptr[pair_marker][1] = *r2;
					matrixptr[pair_marker][2] = variable;
					pair_marker++;

		    		if (s2) r2--; else r2++;
				}
				variable++;

				/* Place 3. */
				pv = pattern[variable];      
				if (s2) r2 += n2-pv+1; else r2 -= n2-pv+1;
				r3 += n3-pv+1;

				for (int k=0; k< pv; k++)		
				{
					matrixptr[pair_marker][0] = *r2;
					matrixptr[pair_marker][1] = *r3++;
					matrixptr[pair_marker][2] = variable;
					pair_marker++;
		    		if (s2) r2++; else r2--;
				}
				variable++;
				interior_1_cell++;
				marker[in_rosette-1] += 1;
				a_map_place ++;
	    	}

	    	/* Type C 2-cell. */
	    	one_cell1 = attaching_map[i][a_map_place];
	    	one_cell2 = attaching_map[i][a_map_place+1];

			bool s3;
	    	if (one_cell1 < 0 )
	    	{
				s2 = false;
				one_cell1 *= -1;
	    	}
	    	else
				s2 = true;

	    	if (one_cell2 < 0 )
	    	{
				s3 = false;
				one_cell2 *= -1;
	    	}
	    	else
				s3 = true;

	    	r1 = labels_for_1_cell[interior_1_cell]->row(0); 
	    	r2 = labels_for_1_cell[one_cell1-1]->row(marker[one_cell1-1]); 
		
	    	if (one_cell1 == one_cell2)
				r3 = labels_for_1_cell[one_cell1-1]->row(marker[one_cell1-1]+1); 
	    	else
				r3 = labels_for_1_cell[one_cell2-1]->row(marker[one_cell2-1]); 

           	n1 = labels_for_1_cell[interior_1_cell]->numcols();
	    	n2 = labels_for_1_cell[one_cell1-1]->numcols();
	    	n3 = labels_for_1_cell[one_cell2-1]->numcols();

	    	/* Place 1. */
	    	pv = pattern[variable];          
	    	r1 += pv-1;
	    	if ( s3 ) r3 += n3-pv; else r3 += pv-1;

	    	for (int j=0; j< pv; j++)
	    	{
				matrixptr[pair_marker][0] = *r1--;
				matrixptr[pair_marker][1] = *r3;
				matrixptr[pair_marker][2] = variable;
				pair_marker++;
				if ( s3 ) r3++; else r3--;
	    	}
	    	variable++;

	    	/* Place 2. */
	    	pv = pattern[variable];             
	    	r1 += n1-pv+1;
	    	if (s2) r2 += pv-1; else r2 += n2-pv;

	    	for (int j=0; j< pv; j++)
	    	{
				matrixptr[pair_marker][0] = *r1++;
				matrixptr[pair_marker][1] = *r2;
				matrixptr[pair_marker][2] = variable;
				pair_marker++;
				if ( s2 ) r2--; else r2++;
	    	}
	    	variable++;

	    	/* Place 3. */
	    	pv = pattern[variable];             
	    	if (s2) r2 += n2-pv+1; else r2 -= n2-pv+1;
	    	if (s3) r3 -= n3-pv+1; else r3 += n3-pv+1;

	    	for (int j=0; j< pv; j++)
	    	{
				matrixptr[pair_marker][0] = *r2;
				matrixptr[pair_marker][1] = *r3;
				matrixptr[pair_marker][2] = variable;
				pair_marker++;
				if (s2) r2++; else r2--;
				if (s3) r3--; else r3++;
	    	}
	    	variable++;
	    	interior_1_cell += 1;
	    	marker[one_cell1-1] += 1;
	    	marker[one_cell2-1] += 1;
	    	/* Again, this code is correct for all cases. */
		}
    }


if (decomp_control::SPLIT_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "split: pair matracies\n";
	for (int i=0; i< number_of_relators; i++)
	{
	    debug << "split:   relator " << i+1;
		print(*pair[i],debug,0,"split:    ");
	}
}

    /*******************************************************************
    We now determine the component arcs on each relator disc.
   
    Declare space for the disc permutations.  For each disc, i, we set up
    a row with maxlb[i]-minlb[i] elements.  Note that each label occurs
    in "pair[i]" exactly once iff the label is on the boundary of the disc 
    and exactly twice iff it is on an interior 1-cell.
    *********************************************************************/

    vector<int>* dperm[number_of_relators]; /*"disc permutations" */
    for( int i =0; i< number_of_relators ; i++)
		dperm[i] = new vector<int> (maxlb[i] - minlb[i] + 1);
	
    /*********************************************************************
    Assign disc permutations.  By moving successively through the rows of 
    "pair[i]" we deal with the "variable corners" in the
    canonical order and, by nature of our construction of "pair[i]", have the 
    following properties satisfied by each row:
   
    1.  If a label in the first column is <='mbl[i]', it is the initial
	label of a component arc on the disc.  If it is
	>'mbl[i]', the pair in that row determines an edge of P whose
	initial point lies on an interior 1-cell.  In this case the label
	at this interior point will have occured exactly once, in the
	second column of a row of "pair[i]", BEFORE the current pair.
   
    2.  If a label in the second column is <='mbl[i]', it is the terminal 
	point of a component arc on the disc.  If it is >'mbl[i]', the 
	pair in that row determines an edge of P whose terminal point
	lies on an interior 1-cell.  In this case the label
	at this terminal point will occur exactly once, in the first
	column of a row of "pair[i]", AFTER the current pair.
   
        The disc permutation will allow us to find the component arcs.  A
    label in the first column of "pair[i]" determines the initial point of
    some edge of P and therefore shall be made to index the label at the
    terminal point of that edge in the permutation.  If a label in the second
    column is <='mbl[i]' it terminates a component arc and, so, shall be made
    to index zero, to signify this.  If the label in the second column is
    >'mbl[i]' it will occur, also, in the first column in a subsequent row and
    so no assignation will be made at this first encounter.  By properties 1
    and 2 above, it can be seen that this procedure completely determines the 
    permutation.
    *********************************************************************/

    for (int i=0;i<number_of_relators; i++)
    {
		matrix<int>& pairs = *pair[i];
		vector<int>& perm = *dperm[i];
		int min = minlb[i];
		int bdry_max = mbl[i];

		for ( int j = 0; j< partial_sum[i]; j++)
		{
		    perm[pairs[j][0] - min] = pairs[j][1] - min;

		    if ( pairs[j][1] <= bdry_max )
				perm[pairs[j][1] - min] = 0;
		}
    }

if (decomp_control::SPLIT_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "split: disc permutations\n";
	for (int i=0; i< number_of_relators; i++)
	{
	    debug << "split:   relator " << i+1 << ": ";
		vector<int>& perm = *dperm[i];

		for (unsigned int j=0; j < perm.size(); j++)
			debug << perm[j] << ' ';
	    debug << endl;
	}
}
    /********************************************************************
    Evaluate the number of components on each disc.  This is simply a matter 
    of counting the terminating points of arcs on the disc, which are
    identified by a unique zero amongst the first 'mbl[i]' elements of the 
    corresponding disc permutation.
    *************************************************************/
    int ndc[number_of_relators]; /* "number of disc components" */
    for ( int i=0; i< number_of_relators ;i++)
		ndc[i] = 0;

    for ( int i =0; i< number_of_relators; i++)
    {
		vector<int>& perm = *dperm[i];
		for ( int j = 0; j < mbl[i] - minlb[i] + 1 ; j++)
		{
		    if ( perm[j] == 0 )
				ndc[i] += 1;
		}
    }
   
    /***************************************************************
    Declare space for disc arcs.  The maximum number of edges an arc may 
    contain is the number, N, of 2-cells in the disc under consideration. 
    The number of labels on an arc is thus <=N+1 and, by storing the 
    terminating zero, we must record a maximum of N+2 integers for each arc.  These
    integers will be held as the rows of a matrix.
    *******************************************************************/
    matrix<int>* darcs[number_of_relators]; /*"disc arcs"*/
    for (int i=0; i < number_of_relators ; i++)
		darcs[i] = new matrix<int>(ndc[i], length_of_attaching_map[i]);
   
    /********************************************************************
    Assign disc arcs.  The initial and terminal labels of each arc may be
    recovered from the indexes of the first 'mbl[i]' components of the 
    corresponding permutation (this number must therefore be even).  
    The terminal labels are indicated by zero and so the labels on each arc may be obtained by 
    "tracing through" the permutation until a zero is encountered.  
    This zero is stored to indicate the termination of that arc.
    *******************************************************************/
    
    for (int i=0; i< number_of_relators; i++)
    {
		vector<int>& perm = *dperm[i];
		matrix<int>& arcs = *darcs[i];
		int min = minlb[i];
		int arcn = 0; /* arc number */

		for ( int j  = 0 ; j < mbl[i] - minlb[i] + 1; j++)
		{
		    if ( perm[j] != 0)
		    {
				arcs[arcn][0] = j + min;
				int arc_cpt = 1; /* "arc component"*/
				bool arc_not_complete = true;
				do
				{
				    arcs[arcn][arc_cpt] = perm[arcs[arcn][arc_cpt-1]-min] + min;
				    arc_cpt++;
	
				    if (arcs[arcn][arc_cpt-1] == min)
				    {
						arcs[arcn][arc_cpt-1] = 0;
						arc_not_complete = false;
				    }
	
				} while (arc_not_complete);
				arcn++;
		    }
		}
	}

if (decomp_control::SPLIT_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "split: number of disc components: ";
	for (int i=0; i< number_of_relators; i++)
		debug << ndc[i] << ' ';
	debug << endl;

	debug << "split: disc arcs" << endl;
	for (int i=0; i< number_of_relators; i++)
	{
	    debug << "split:   relator " << i+1 << ":" << endl;
	    print(*darcs[i],debug,0,"split:     ");
	}
}

       
    /********************************************************************
    It remains to determine how the disc arc components are identified when 
    the relator discs are attached to the rosette.  To accomplish this we 
    utilize a matrix whose rows are made up of the end labels of each disc 
    component.
    ********************************************************************/
    int tndc = 0; //"total number of disc components" 
    for ( int i =0; i<number_of_relators ; i++)
		tndc += ndc[i];

    matrix<int> ael(tndc,2); /* "arc end labels" */
    int elm = 0; /* "end label marker" */

    for ( int i =0; i< number_of_relators ; i++)
    {
		matrix<int>& arcs = *darcs[i];
		for (int j=0; j< ndc[i]; j++)
		{
		    ael[elm][0] = arcs[j][0];
		    int k = 1;
		    do  
				k++; 
		    while (arcs[j][k] != 0);
			
		    /* We know there are at least two labels on each arc.  When 
		    this loop terminates, k indicates the terminating zero of the
		    current  arc. */

		   ael[elm][1] = arcs[j][k-1];
		   elm++;
		}
    }
   

    /********************************************************************
    We build up components of P as classes of boundary labels, by copying 
    required columns from the label matrices into a row of a component
    matrix.  The  labels in these columns are identified by the attaching
    maps, as discussed above. Given one column, we may find the other 
    columns that contain labels, identified in L, belonging to the same 
    component of P by moving  along an arc in some relator disc from each 
    label in turn, ie examining the corresponding place in 'ael'.  We shall 
    set elements of 'ael' to be negative to prevent us obtaining the same 
    column twice thus enabling a more efficient method of constructing the 
    classes.
   
	 We have 'tndc' as an upper bound for the number of components and 
    the number of boundary labels as an upper bound for the size of a 
    component.  Note the labels on interior 1-cells are irrelevant to this 
    procedure, since no identification of  discs occurs along interior 
    1-cells.
    *********************************************************************/
    /* calculate the number of boundary labels */
    label = 0;
    for (int i=0; i< number_of_relators; i++)
    	label += mbl[i] - minlb[i] + 1;

    matrix<int> cptm(tndc,label+1); /* "component matrix" We need to add 1
				    to this total to ensure we have a zero
				    at the end of each row of cptm, for 
				    later use */
 
    /*Build up component classes.*/
    int cpt = 0;
  	bool all_components_found;
	
    do /* while !all_components_found */
    {
		/* Find a boundary label not already recorded. */
		all_components_found = true;
		int new_label;
		for (int i = 0; i< tndc; i++)
		{
			for (int j = 0; j < 2; j++)
			{
		    	if (ael[i][j] > 0)
		    	{
					all_components_found = false;
					new_label = ael[i][j];
					goto new_label_found;
	    		}
			}
		}
		new_label_found:
		
		if (!all_components_found)
		{
	    	int celt = 0; // component element
	    	Slice_iter<int> pointer = cptm[cpt]; /* pointer will reference a label in the current comonent */
	    	bool component_complete = false;
	    	do /* While !component_complete */
	    	{
				int one_cell;
				int column;
				
				/* Look for this new label in the label matrices. */
				for (int i=0; i< number_of_1_cells_in_rosette; i++)
				{
		    		matrix<int>& matrixptr = *labels_for_1_cell[i];
					
		    		for (unsigned int j=0; j < matrixptr.numrows(); j++)
		    		{
						for (unsigned int k=0; k< matrixptr.numcols(); k++)
						{
				    		if ( matrixptr[j][k] == new_label)
				    		{
								one_cell = i;
								column = k;
								goto found_in_label_matrices;
				    		}
						}
		    		}
				}
				found_in_label_matrices:

				 /************************************************************************				   
				 Set those integers in 'ael' that appear in the column just found to be 
				 negative, so we do not look for this column again.  These integers will 
				 currently be positive, since '*new_label' was positive and each label 
				 occurs in a unique place amongst the label matrices.  The integers may be 
				 found by the symmetry of 'jba'; at the same time we copy this column 
				 into the current component.
				 ************************************************************************/

	    		matrix<int>& matrixptr = *labels_for_1_cell[one_cell];
				for (unsigned int i=0; i< matrixptr.numrows(); i++)
				{
		    		int label = matrixptr[i][column];
					for (int i = 0; i< tndc; i++)
					{
						for (int j = 0; j < 2; j++)
						{
					    	if (ael[i][j] == label)
					    	{
								ael[i][j] *= -1;
								goto found_in_arc_end_labels;
				    		}
						}
					}
					found_in_arc_end_labels:
		    		
					cptm[cpt][celt++] = label;
				}

				/* Now find the next label to look for. */
				bool label_found = false;
				do
				{
		    		if ( *pointer == 0 )
					component_complete = true;
		    		else
		    		{
						/* Look for *pointer in ael.  When we find it,
						   we check to see if we already have the labels
						   at the other end of this arc.  If we have them,
						   increment pointer and look again */
						   
						for (int i=0; i< tndc; i++)
						{
							for (int j=0; j< 2; j++)				
							{
				    			if ( abs(ael[i][j]) == *pointer )
					    		{
									if ( j == 0 && ael[i][1] > 0)
									{
							    		new_label = ael[i][1];
							    		label_found = true;
									}
									else
									if ( j == 1 && ael[i][0] > 0) 
									{
							    		new_label = ael[i][0];
							    		label_found = true;
									}
									break;
					    		}
							}
						}

						if (!label_found)
			    			pointer++;
		    		}
				} while(!component_complete && !label_found);
	    	} while (!component_complete);
	    	cpt++;
		}
    } while (!all_components_found);

if (decomp_control::SPLIT_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "split: total number of disc arcs = " << tndc << endl;
	debug << "split: arc end labels:" << endl;
	print(ael,debug,0,"split:   ");
	debug << "split: component matrix:" << endl;
	print(cptm,debug,0,"split:  ");
}
	
	/* we can now do some clearing up of memory*/
	for (int i=0; i< number_of_1_cells; i++)
		delete labels_for_1_cell[i];

    /*********************************************************************
    We may now use the rows of 'cptm' to determine the tuples corresponding to 
    each component of P.  Each label of a disc arc component occurs exactly 
    once in 'cptm', clearly with both end labels of a given arc in the same 
    component.  We therefore determine which component each "initial" label 
    lies in and then use that disc arc to build up the tuple for the   
    corresponding component.  Notice 'cpt' now indicates how many components 
    P has.
    **********************************************************************/

    /* Check whether the components are required by call to function. */
    if ( short_split )
		return (cpt == 1? true: false);
		 
	if (!return_components)
	{
	    cout << "\n";
		for (unsigned int i=0; i<pattern.size(); i++)
			cout << pattern[i] << ' ';
		cout << endl;

    	txt << "\n";
		for (unsigned int i=0; i<pattern.size(); i++)
			txt << pattern[i] << ' ';
		txt << endl;
	}

	matrix<int> ctracks(cpt, pattern.size()); /* "component tracks" */

	int b = cptm.numcols();
	for (int i=0; i< number_of_relators; i++)
	{
		matrix<int>& arcs = *darcs[i];
		matrix<int>& pairs = *pair[i];

		for (int j=0;  j < ndc[i]; j++)
		{
			/* Determine in which component the initial label of 'arcs[j,]' lies. */
			bool label_found = false;
			int mark = arcs[j][0];
			int disc_cpt; // "disc_component"
			int pattern_cpt; // "pattern component"

			for (int k=0; k < cpt; k++)
			{
				for (int l=0; l < b; l++)
		    	{
					label = cptm[k][l];
					if ( label == mark )
					{
						disc_cpt = j;
						pattern_cpt = k;
						label_found = true;
						break;
					}
					else
					if ( label == 0 )
						break;
				} 

				if (label_found)
					break;
			}

			 /**********************************************************************
			 We now traverse the arc 'arcs[disc_cpt,]', at each stage determining with which 
			 variable the corresponding edge is associated. In the loop below, 
			 'arcs[disc_cpt,k]' is the initial label of some edge on disc i and so it appears 
			 exactly once in the first column of 'pairs'.  We now determine in which 
			 row, then, the third element of that row determines the variable 
			 associated with the pair, so we add 1 to the corresponding variable in 
			 'ctracks[in cpt,]'.
			 ***********************************************************************/

			for (unsigned int k = 0; k < arcs.numcols(); k++)
			{
				label = arcs[disc_cpt][k];
				if (label == 0)
					break;

		    	label_found = false;
				for (int l = 0; l< partial_sum[i]; l++)
				{
					if ( pairs[l][0] == label )
					{
						ctracks[pattern_cpt][pairs[l][2]] += 1;
						label_found = true;
						break;
					}
				}
			}
	    }
	}
		
if (decomp_control::SPLIT_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "split: component tracks" << endl;           
	print(ctracks,debug,0,"split:   ");
}

	if (return_components)
		*component_ptr = new matrix<int>(ctracks);
	else
	{
	    if ( cpt == 1 )
		{
			cout << "which is a track\n";
			txt << "which is a track\n";
		}
		else
		{
			cout <<	"which has components\n";           
			txt <<	"which has components\n";           
			cout << ctracks << endl;
			txt << ctracks << endl;
		}
	}

    /* This case is 'short_split = false', we just need to return 
       a bool which will be discarded at the call. */
	
	for (int i=0; i< number_of_relators; i++)
	{
		delete pair[i];
		delete dperm[i];
		delete darcs[i];
	}

    return true;
}

int find_generator(vector<char> generators, char c)
{

if (decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "find_generator: generators: ";
	for (unsigned int i=0; i< generators.size(); i++)
		debug << generators[i] << ' ';
	debug << endl;
	debug << "find_generator: c = " << c << endl;
}	
	int place = -1;
	for (unsigned int i=0; i < generators.size(); i++)
	{
		if (generators[i] == c)
		{

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "find_generator: found c at place " << i << endl;
			place = i;
			break;
		}
	}

if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "find_generator: returning " << place << endl;

	return place;
}

/*************************************************************************
This procedure extends the current tree over the component of the boundary 
of B containing a b-label assigned to 'new label' at the call.  The 
algorithm is essentially that used for the t-tree. 
**************************************************************************/
void band_extend
(
	int new_label, 
	int number_of_1_cells_in_rosette, 
	int* occurrences_of_1_cell, 
	matrix<int>** blabels_for_1_cell, 
	int& lpointer, 
	matrix<int>& bposition, 
	unsigned int& ppointer, 
	int& pmarker, 
	matrix<int>& jbba, 
	matrix<int>& btrarcs, 
	int& amarker, 
	bool& twisted
)
{
    bool tree_complete = false;
    do
    {
		int one_cell;
		int column;
		/* Look for new_label in the blabel matrices. */
		for (int i=0; i< number_of_1_cells_in_rosette; i++)
		{
	    	matrix<int>& matrixptr = *blabels_for_1_cell[i];
	    	for (unsigned int j=0; j < matrixptr.numrows(); j++)
	    	{
				for (unsigned int k=0; k< matrixptr.numcols(); k++)
				{
		    		if ( matrixptr[j][k] == new_label)
		    		{
						bposition[pmarker][0] = i;
						bposition[pmarker][1] = k;
						pmarker++;
						
						one_cell = i;
						column = k;
						goto new_label_found_in_blabel_matrices;
		    		}
				}
	    	}
		}
		new_label_found_in_blabel_matrices:	
	
	    /************************************************************************				   
	     Set those elements in 'jbba' column 2 that appear in the column just found 
	 	 to be negative, so we do not look for this column again.  These integers 
	     will currently be positive, since '*new_label' was positive and each label 
	     occurs in a unique place amongst the label matrices.  We also look out for 
	     the b-label 2.
	 ************************************************************************/

    	matrix<int>& matrixptr = *blabels_for_1_cell[one_cell];
		for (unsigned int i=0; i< matrixptr.numrows(); i++)
		{
	    	int blabel = matrixptr[i][column];
			jbba[abs(jbba[blabel-1][1])-1][1] *= -1;
	    	if (blabel == 2)
				twisted = true;
		}

		/* We now see how we can extend the tree, ie find a "new label" */
		bool label_found = false;  
		do
		{
			int one_cell = bposition[ppointer][0]; // a 1-cell in the rosette
			int column = bposition[ppointer][1];
			
			lpointer++;
	    	if ( lpointer > occurrences_of_1_cell[bposition[ppointer][0]])
	    	/* ie if lpointer goes off the end of its current column. */
	    	{
				ppointer++;
				lpointer = 1;
	    	}

	    	if (ppointer >= bposition.numrows() || bposition[ppointer][0] == -1)
			{
				tree_complete = true;
			}
			else
			{
				one_cell = bposition[ppointer][0]; // a 1-cell in the rosette
				column = bposition[ppointer][1];
			}

	    	if (!tree_complete)
	    	{
				matrix<int>& blabels = *(blabels_for_1_cell[one_cell]);

				int il = blabels[lpointer - 1][column];

				/*******************************************************
				"initial label" is the label referred to by 'lpointer' 
				in the current column.
				*****************************************************/
				if (jbba[il-1][1] > 0 )
				{
		    		new_label = jbba[il-1][1];
		    		label_found = true;
		    		btrarcs[amarker][0] = il;
		    		btrarcs[amarker][1] = new_label;
					amarker++;
				}
	    	}
		} while(!label_found && !tree_complete);
    } while(!tree_complete);
}

/***********************************************************************
			    add_0_cells
This procedure adds as many 0-cells as possible to the current tree.  We
add to 'btrarcs1(2)' only one occurrence of an edge of EX (joining a
0-cell b-label to the tree) in the boundary of some relator disc but
assign 'jbba[b][1]=-1' for each 0-cell b-label, b, attached to a 0-cell
added to the tree, to signify this fact.  Notice, if t separates, on the
second call to 'add 0 cells', this will indicate that this 0-cell has
been added to the first tree.  (No confusion with the 1-cell b-label 1
will arise.)  

The first 0-cell of L lies in the same component of L-t as the b-label 1, 
so we shall always be able to add the first 0-cell to the
first tree.  By considering the rosette, R(L), we see that there is
always an edge of EX joining some 0-cell b-label (attached to the first
0-cell) to a 1-cell b-label in the first tree.  

2006 comment: let c be the 1-cell containing b-label 1.  If the generator corresponding 
to c is not a 2-torsion generator, or if c=a1 for some 2-torsion generator a whose 1-cells 
are a1,a2, then the above paragraph is clear.  If c=a2 for the 2-torsion generator then from
the relator disc for a^2=1 there is an arc joining b-label 0 to jbba[0][1], the blabel at the
other end of the disc-band-arc from blanel 1, so again blabel 0 lies in the same component of
L-t as blabel 1 (t does not intersect a1 in this last case).


Notice that it is possible
that a 0-cell (not the first) may only be added to the (first) tree by
joining it, by an edge of EX, to the first 0-cell.  This only occurrs if
there is a 2-torsion generator with corresponding relator disc having
void intersection with t (as can be seen from R(L), since we have
constructed a maximal tree in the component of the boundary of B
containing the b-label 1).  In this case we shall obtain the first 1-cell
"carried" by this generator as the joining edge of EX, since we are
insisting that 2-torsion generators are written in the positive at all
times and 2-torsion relators appear first in R.
*************************************************************************/
void add_0_cells(int number_of_0_cells, vector<int>** blabels_for_0_cell, matrix<int>& jbba, 
				 int& amarker, matrix<int>& btrarcs)
{

    for ( int i=0; i< number_of_0_cells; i++)
    {
		/* Isolate the b-labels on 0-cell i. */
		vector<int> blabels = *blabels_for_0_cell[i];

		/* check if this 0-cell has been not been added to the tree. */
		if  ( jbba[blabels[0]-1][1] == 0 )
		{
	    	/* Check whether any of the b-labels referenced by 'blabels
	    	   are joined to the tree by an edge of EX. */
	    	for (unsigned int j=0; j < blabels.size(); j++)
	    	{
				if ( jbba[jbba[blabels[j]-1][0]-1][1] < 0 )

/************************************************************************
In 1987, the above condition was recorded as :

		   IF     IF   jbba[row[j],1]>max1cbl
			       {ie, if 'jbba[row[j],1]' is a 0-cell b-label.}
			  THEN jbba[jbba[row[j],1],2]<0 {1}
			  ELSE jbba[ABS jbba[jbba[row[j],1],2],2]<0 {2}
                          FI

1. The 0-cell b-label, b, is in the tree iff 'jbba[b][1]' has been
   assigned -1.

2. If b:=jbba[blabels[j]-1][0] is a 1-cell b-label, 'ABS jbba[b-1][1]' is the 1-cell
   b-label joined to b by a 'band-arc component' of the boundary of B.  We
   inspect the 2nd (index 1) column of this row in 'jbba' (initially assigned b) to
   determine whether b is in the tree.


However, in 2.,if b is a 1-cell b-label, then b is in the tree iff 
jbba[b-1][1] is in the tree, by the way the tree is constructed, so we only 
need to check whether jbba[b-1][1] < 0 to determine whether b is in the tree.  
This observation would have produced the Algol code

			  ELSE jbba[jbba[row[j],1],2]<0 {2}

which is the same condition as {1}; the C has been written accordingly.

If in 2 b is a 0-cell b-label, then it is in the tree iff jbba[b-1][1] = -1, so the same 
test suffices.
************************************************************************/
				{
		    		/* blabels[j] is joined to jbba[blabels[j]-1][0], which is in the tree. */
					btrarcs[amarker][0] = blabels[j];
					btrarcs[amarker][1] = jbba[blabels[j]-1][0];
					amarker++;

		    		for (unsigned int k = 0; k < blabels.size(); k++)
						jbba[blabels[k]-1][1] = -1;
						
		    		break;
	    		}
				else
				if ( jbba[jbba[blabels[j]-1][2]-1][1] < 0 )
				{
		    		/* blabels[j] is joined to jbba[blabels[j]-1][2], which is in the tree. */
					btrarcs[amarker][0] = blabels[j];
					btrarcs[amarker][1] = jbba[blabels[j]-1][2];
					amarker++;

		    		for (unsigned int k = 0; k < blabels.size(); k++)
						jbba[blabels[k]-1][1] = -1;
						
		    		break;
				}
			}
		}
    }
}


/* find_cpt is a version of the find algorithm for finding tracks
   in pattern component lists.  It was implemented explicitly 
   as a test of MJD MAC compiler */
list<ray>::iterator find_cpt(list<ray>& cpt_list, ray cpt)
{
	unsigned int n = cpt.r.size();
	list<ray>::iterator rptr = cpt_list.begin();

if (decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "find_cpt: ray size = " << n << endl;
	debug << "find_cpt: looking for component ";
	for (unsigned int i=0; i< n; i++)
		debug << cpt.r[i] << ' ';
	debug << endl;
	debug << "find_cpt: cpt_list:" << endl;
	
	while (rptr != cpt_list.end())
	{
		debug << "find_cpt:   ";
		for (unsigned int i=0; i< n; i++)
			debug << rptr->r[i] << ' ';
		debug << endl;

		rptr++;
	}
}	
	
	rptr = cpt_list.begin();
	
	while (rptr != cpt_list.end())
	{

if (decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "find_cpt: checking component ";
	for (unsigned int i=0; i< n; i++)
		debug << rptr->r[i] << ' ';
	debug << endl;
}

		bool equal = true;
		
		for (unsigned int i=0; i< n; i++)
		{
			if (rptr->r[i] != cpt.r[i])
			{
				equal = false;			
if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "find_cpt:   i = " << i << " component[i] = " << rptr->r[i] << " ray[i] = " << cpt.r[i] << (equal? " equal": " different") << endl;
				
				break;
			}
if (decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "find_cpt:   i = " << i << " component[i] = " << rptr->r[i] << " ray[i] = " << cpt.r[i] << (equal? " equal": " different") << endl;

		}
		
		if (equal)
			break;
		else
			rptr++;
	}
	
if (decomp_control::DEBUG >= decomp_control::DETAIL)
{
	if (rptr == cpt_list.end())
		debug << "find_cpt: not found, returning end() iterator" << endl;
	else
		debug << "find_cpt: found, returning iterator pointing to component" << endl;
}
	return rptr;
}

/* Two patterns P1 and P2 are compatible if the components of P1+P2 are the
   disjoint union of the components of P1 and P2 (see Dicks & Dunwoody p 235)
   That is, each component of P1+P2 is a component of either P1 or P2, or both,
   conversely, if a component of P1+P2 is neither a component of P1 nor of P2,
   then P1 and P2 are incompatible.
*/
bool compatible (ray& P1, ray& P2, K_complex K)
{
	bool result = true;
	int n = P1.r.size();
	
if (decomp_control::DEBUG >= decomp_control::INTERMEDIATE)
{
	debug << "compatible: P1 = ";
	for (int i=0; i< n; i++)
		debug << P1.r[i] << ' ';
	debug << endl << "compatible: P2 = ";
	for (int i=0; i< n; i++)
		debug << P2.r[i] << ' ';
	debug << endl;
}

	/* Break P1 into its components and create a list of them */
	matrix<int>* P1_Components;
	ofstream dummy_ofstream;

	split (P1.r, K, false, dummy_ofstream, true, &P1_Components);

	list<ray> P1_list;
	matrix<int>& P1_components = *P1_Components;
	
/*if (decomp_control::DEBUG >= decomp_control::INTERMEDIATE)
{
	debug << "P1_components =\n\t";
	for (unsigned int i=0; i< P1_components.numrows(); i++)
	{
		debug << "\n\t";
		for (unsigned int j=0; j< P1_components.numcols(); j++)
			debug << P1_components[i][j] << ' ';
	}
	debug << endl;
}*/

	for (unsigned int i = 0; i< P1_components.numrows(); i++)
	{
		vector<int> ray(n);
		
		for (int j = 0; j<n; j++)
			ray[j] = P1_components[i][j];
		
		P1_list.push_back(ray);
	}		

if (decomp_control::DEBUG >= decomp_control::INTERMEDIATE)
{
	debug << "compatible: P1_list:" << endl;
	list<ray>::iterator P1_ptr = P1_list.begin();
	
	while (P1_ptr != P1_list.end())
	{
		ray& R = *P1_ptr;
		
		debug << "compatible:   ";
		for (unsigned int j=0; j< R.r.size(); j++)
			debug << R.r[j] << ' ';
		debug << endl;
		P1_ptr++;
	}
}

	matrix<int>* P2_Components;

	split (P2.r, K, false, dummy_ofstream, true, &P2_Components);

	list<ray> P2_list;
	matrix<int>& P2_components = *P2_Components;
	
/*if (decomp_control::DEBUG >= decomp_control::INTERMEDIATE)
{
	debug << "P2_components =\n\t";
	for (unsigned int i=0; i< P2_components.numrows(); i++)
	{
		debug << "\n\t";
		for (unsigned int j=0; j< P2_components.numcols(); j++)
			debug << P2_components[i][j] << ' ';
	}
	debug << endl;
}*/

	for (unsigned int i = 0; i< P2_components.numrows(); i++)
	{
		vector<int> ray(n);
		
		for (int j = 0; j<n; j++)
			ray[j] = P2_components[i][j];
		
		P2_list.push_back(ray);
	}		

if (decomp_control::DEBUG >= decomp_control::INTERMEDIATE)
{
	debug << "compatible: P2_list:" << endl;
	list<ray>::iterator P2_ptr = P2_list.begin();
	
	while (P2_ptr != P2_list.end())
	{
		ray& R = *P2_ptr;
		
		debug << "compatible:   ";
		for (unsigned int j=0; j< R.r.size(); j++)
			debug << R.r[j] << ' ';
		debug << endl;		
		P2_ptr++;
	}
}

	/* create the sum P = P1 + P2 and split that to form a list of components as above */
	vector<int> P = P1.r;
	for (int i=0; i < n; i++)
		P[i] += P2.r[i];

if (decomp_control::DEBUG >= decomp_control::INTERMEDIATE)
{
	debug << "compatible: sum = ";
	for (int i=0; i< n; i++)
		debug << P[i] << ' ';
	debug << endl;
}

	matrix<int>* P_Components;
	
	split (P, K, false, dummy_ofstream, true, &P_Components);
	
	list<ray> P_list;
	matrix<int>& P_components = *P_Components;
	
/*if (decomp_control::DEBUG >= decomp_control::INTERMEDIATE)
{
	debug << "P_components =\n\t";
	for (unsigned int i=0; i< P_components.numrows(); i++)
	{
		debug << "\n\t";
		for (unsigned int j=0; j< P_components.numcols(); j++)
			debug << P_components[i][j] << ' ';
	}
	debug << endl;
}*/

	for (unsigned int i = 0; i< P_components.numrows(); i++)
	{
		vector<int> ray(n);
		
		for (int j = 0; j<n; j++)
			ray[j] = P_components[i][j];
		
		P_list.push_back(ray);
	}		


if (decomp_control::DEBUG >= decomp_control::INTERMEDIATE)
{
	debug << "compatible: P_list:" << endl;
	list<ray>::iterator P_ptr = P_list.begin();
	
	while (P_ptr != P_list.end())
	{
		ray& R = *P_ptr;
		
		debug << "compatible:   ";
		for (unsigned int j=0; j< R.r.size(); j++)
			debug << R.r[j] << ' ';
		debug << endl;		
		P_ptr++;
	}
}

	if (P_list.size() != P1_list.size() + P2_list.size())
	{
if (decomp_control::DEBUG >= decomp_control::INTERMEDIATE)
	debug << "compatible: number of components of sum is not the sum of number of components of P1 and P2" << endl;
		result = false;
	}
	else
	{
		/* for each component c of P, we check that c is a components of either P1 or P2, 
		   and remove it from the corresponding list if it is found.  If c is not found in 
		   either P1 or P2 then the patterns are incompatible.*/
		list<ray>::iterator P1_ptr;
		list<ray>::iterator P2_ptr;
		list<ray>::iterator P_ptr = P_list.begin();
		  
		while (P_ptr != P_list.end())
		{
if (decomp_control::DEBUG >= decomp_control::INTERMEDIATE)
{
	debug << "compatible: looking for component ";
	for (unsigned int i=0; i< P_ptr->r.size(); i++)
		debug << P_ptr->r[i] << ' ';
	debug << endl;
}
//			P1_ptr = find(P1_list.begin(), P1_list.end(), *P_ptr);
//			P2_ptr = find(P2_list.begin(), P2_list.end(), *P_ptr);

if (decomp_control::DEBUG >= decomp_control::INTERMEDIATE)
	debug << "compatible: checking P1" << endl;
			P1_ptr = find_cpt(P1_list, *P_ptr);

if (decomp_control::DEBUG >= decomp_control::INTERMEDIATE)
	debug << "compatible: checking P1" << endl;
			P2_ptr = find_cpt(P2_list, *P_ptr);
			
			if (P1_ptr == P1_list.end() && P2_ptr == P2_list.end())
			{

				/* component cannot be found in either pattern, so they're not compatible */
if (decomp_control::DEBUG >= decomp_control::INTERMEDIATE)
	debug << "compatible: not found in either pattern, so they are incompatible" << endl;
			
				result =  false;
				break;
			}
			else if (P1_ptr != P1_list.end())
			{
if (decomp_control::DEBUG >= decomp_control::INTERMEDIATE)
	debug << "compatible: found in pattern P1" << endl;
				/* remove the component from P1_list */
				P1_list.erase(P1_ptr);
			}
			else // component must be in P2
			{
if (decomp_control::DEBUG >= decomp_control::INTERMEDIATE)
	debug << "compatible: found in pattern P2" << endl;
				/* remove the component from P2_list */
				P2_list.erase(P2_ptr);

			}
			
			P_ptr++;
		}	  
	}
		
	delete P1_Components;
	delete P2_Components;
	delete P_Components;
	
if (decomp_control::DEBUG >= decomp_control::INTERMEDIATE)
{
	debug << "compatible: patterns are ";
	
	if (!result)
		debug << "not ";
	debug << "compatibile" << endl;
}

	return result;
}
