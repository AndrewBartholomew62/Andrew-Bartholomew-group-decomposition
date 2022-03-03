using namespace std;

#include <fstream>
#include <iostream>
#include <cstdlib>
#include <list>
#include <string>
#include <ctype.h>
#include <stdio.h>
#include <iomanip>
#include <vector>

extern unsigned int	DEBUG;

extern bool		TWISTED_ALLOWED;
extern bool		FLAGS_ON_FIRST_2_TORSION_1_CELL;
extern ofstream	debug;

#define NOT_FOUND -1

#include <util.h>
#include <matrix.h>
#include <decomp.h>

/********** Function prototypes *************************/
void band_extend(int new_label, int number_of_1_cells_in_rosette, int* occurrences_of_1_cell, 
				 matrix<int>** blabels_for_1_cell, int& lpointer, matrix<int>& bposition, unsigned int& ppointer, 
				 int& pmarker, matrix<int>& jbba, matrix<int>& btrarcs, int& amarker, bool& twisted);

void add_0_cells(int number_of_0_cells, vector<int>** blabels_for_0_cell, matrix<int>& jbba, 
				 int& amarker, matrix<int>& btrarcs);
void find(int label,int number_of_1_cells_in_rosette, matrix<int>** labels_for_1_cell, vector<int>& position);
int bvertex(int label, int number_of_0_cells, int numbcols, vector<int>** blabels_for_0_cell);
vector<int>* find_row (int label, int number_of_0_cells, vector<int>** blabels_for_0_cell);
//void print(matrix<int> m, ostream& s, int n, string prefix="");
int find_generator(vector<char> generators, char c);

string disc_word
(
	int 				l1, 
	int 				l2, 
	bool 				t_labels, 
	bool 				in_a_tree,
    vector<int>& 		max_1_cell_labels, 
	matrix<int>** 		one_cell_labels,
	vector<string>&		flag_on_1_cell, 
	int* 				n_tpts_on_1_cell, 
	matrix<int>& 		ext0cbl, 
	int 				max1cbl,
	int 				number_of_1_cells_in_rosette, 
	vector<int> 				length_of_attaching_map, 
	matrix<int>& 		attaching_map
);

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
	vector<string>&		flag_on_1_cell, 
	int* 				n_tpts_on_1_cell, 
	matrix<int>& 		ext0cbl, 
	int 				max1cbl,
	int 				number_of_1_cells_in_rosette, 
	vector<int>			length_of_attaching_map, 
	matrix<int>& 		attaching_map
);

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
	vector<int> 	length_of_attaching_map,	
	matrix<int>& 	attaching_map,
	matrix<int>& 	ext0cbl
);

string simple_reduce(string str);
void reduce(string* generators, int first_gen, int num_gens, vector<string>& relator, int number_of_relators);

/********************************************************************
decompose determines the decomposition of G given by the track t 
(:=P('track')), writing the results to the file 'txt'.  The strategy is as 
follows.  We construct two sets of labels on the relator discs: firstly 
"t-labels".  For each relator disc, D (taken in the canonical order) we
label the points of intersection of t with the 1-skeleton of D
with consecutive integers, starting with 1 and labelling clockwise around
the boundary of D, beginning at O(D), then continuing for each interior
1-cell of D (in order) starting at O(D).  The other labels, "b-labels" are
constructed only on the boundary of the relator discs.  Thicken t to a band
B, such that t=t(B) and label, clockwise round the boundary of the discs,
taken in order and starting at O(D) on each disc, the points of intersection
of the boundary of B and the boundary of the discs with consecutive integers,
starting with 1.  Finally, extend the b-labelling to the 0-cells in the
boundary of each disc, again via the canonical order.  The t-labels mark the
vertices of t, regarded as a graph with Vt the intersection of t with the
1-skeleton of L, and the b-labels mark the 0-cells of the cell decomposition
described in chapter 5 for L-t.  The t-labels are used to build a
representation of t, regarded as a graph with Vt the intersection of t with
R(L), which allows the construction of a maximal tree in this "t-graph".
From the b-labels and the representation of  t we are able to build a
complete representation of the 1-skeleton (X, say) of our cell decomposition
for L-t.  This enables us to construct a  maximal tree in each component of
X, the decision on connectivity actually being made during the construction.

      The decomposition is given by the induced images of the fundamental 
groups of t and the components of L-t in G (induced by inclusion).  We are 
able to calculate generators for these fundamental groups, since we now have 
the required maximal trees:  moreover, the relator discs allow us to 
calculate the induced image of these generators thus giving a generating set 
for the induced image, which is what is required from the procedure.

The boolean full_decomp was added 30/8/09 to provide reuse of the decompose code to 
report on the status of a track, being twisted, non-separating or separating.  If 
full_decomp is false, the decompose function returns true after reaching the point
where these characteristics are known and reported.
*************************************************************************/
unsigned int decompose (vector<int>& track, K_complex& K, bool echo, bool full_decomp, ofstream& txt)     
{
	if (!full_decomp)
		TWISTED_ALLOWED = true; // otherwise twisted tracks will be doubled
	
	unsigned int status=0; // used for the return value
		
    /**********  read variables required by decompose from the K_complex **********/
	int number_of_0_cells = K.number_of_0_cells;
	int number_of_1_cells = K.number_of_1_cells;
	int number_of_2_cells = K.number_of_2_cells;
	int number_of_generators = K.number_of_generators;
	int number_of_relators = K.number_of_relators;
	vector<string> relator = K.relator;
	int occ_of_cv = K.occ_of_cv;
	vector<char> generator_corresponding_to_1_cell = K.generator_corresponding_to_1_cell;
	int number_of_1_cells_in_rosette = K.number_of_1_cells_in_rosette;
	vector<char> generators_of_order_2 = K.generators_of_order_2;
	vector<int> occ_of_order_2_gen = K.occ_of_order_2_gen;
	int* occurrences_of_1_cell = K.occurrences_of_1_cell;
	vector<int> length_of_attaching_map = K.length_of_attaching_map;
	matrix<int>& two_cell = *K.two_cell;
	matrix<int>& attaching_map = *K.attaching_map;

	
    /************************************************************************
    For the 1-cells of the rosette, the total number of occurrences of 1-cell 
    i in 'two cell' is the number of occurrences of i in the boundary of some 
    relator disc.  For the t-labels, we therefore set out, for each such 
    1-cell, a matrix with r rows and c columns, where r is the number of 
    occurrences of the 1-cell in the boundary of the relator discs and c is 
    the number of points of t intersecting that 1-cell.  For the interior 
    1-cells, which occur exactly once in a unique relator disc, we set out a 
    1xc matrix, where c is as above.  For the b-labels, we set out, for each 
    1-cell in the rosette, a matrix with r rows and 2c columns (r and c as 
    above) and for the 0-cells a row of n integers, where n is the number of 
    occurrences of the 0-cell in the boundary of the relator discs.  We
    store, in  these matrices, the labels attached to the 1-cells
    corresponding to their rows (in accordance with the canonical order
    through the relator discs, ignoring identifications) in such a way that
    the columns of each matrix corresponding to a 1-cell in the rosette
    determines how those labels become identified in L (of course, the row is
    identified with the 0-cells).
    ************************************************************************/

if (decomp_control::DECOMPOSE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "decompose: decomposing the track: ";
	for (unsigned int i=0; i< track.size(); i++)
		debug << track[i] << ' ';
	debug << endl;
}

    /*  First deal with the 1-cells in the rosette. */

	matrix<int>* tlabels_for_1_cell[number_of_1_cells];
	matrix<int>* blabels_for_1_cell[number_of_1_cells_in_rosette];
    int n_tpts_on_1_cell[number_of_1_cells];
    /* "number of track points on 1-cell" will be used later. */

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
						n_tpts_on_1_cell[i] = track[3*j+k] + track[3*j];
		    		else 
						n_tpts_on_1_cell[i] = track[3*j+k] + track[3*j+k+1];

		    		one_cell_found = true;
		    		break;
				}  	     
	    	} 
			
	    	if (one_cell_found)
			break;
		}

		tlabels_for_1_cell[i] = new matrix<int> (occurrences_of_1_cell[i],n_tpts_on_1_cell[i]);	      
		blabels_for_1_cell[i] = new matrix<int> (occurrences_of_1_cell[i],2*n_tpts_on_1_cell[i]);	      
    }


    /* Now deal with the t-labels on the interior 1-cells, for which purpose 
    we shall need a different marker for 'labels for 1 cell'.  We need to 
    know the first variable on each relator disc, so evaluate these first. */


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
	    	n_tpts_on_1_cell[interior_1_cell] = track[j] + track[j+1];

		   tlabels_for_1_cell[interior_1_cell] = new matrix<int> (1,n_tpts_on_1_cell[interior_1_cell]);   
		   interior_1_cell ++;
		}
    }

if (decomp_control::DECOMPOSE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "decompose: number of track points on 1-cells: ";
	for (int i=0; i<number_of_1_cells; i++)
		debug << n_tpts_on_1_cell[i] << ' ';
	debug << endl;
}


    /**********************************************************************
    Assign the t-(b-)labels to the space allocated.  We take the convention 
    that, for the 1-cells in the rosette, the orientation of the labels in the
    label matrix is left to right if the 1-cell is "positive in the attaching 
    map" and right to left otherwise.  This ensures the identification
    matching mentioned above.  We also record maximum and minimum t-label on
    each disc, the maximum t-label on the boundary of each disc, and the
    maximum 1-cell b-label on each disc, all for future use.
    ************************************************************************/

    int bmarker[number_of_1_cells_in_rosette];
    /* These markers will be used to mark the rows in the b0cbl_for_1_cell matracies. */

    int marker[number_of_1_cells_in_rosette];
    /* These markers will mark the rows in the label matrices. */

    for (int i=0; i< number_of_1_cells_in_rosette; i++)
    {
		marker[i] = 0;
		bmarker[i] = 0;
    }

    int tlabel = 1;
	int blabel = 1;
	
    /* 'one cell' will denote a 1-cell of the rosette; tlabel and blabel will 
      be used to allocate labels. */

    interior_1_cell = number_of_1_cells_in_rosette;
    int mintlb[number_of_relators];
    int maxtlb[number_of_relators];   
    vector<int> mbtl(number_of_relators);
    vector<int> m1cbl(number_of_relators);
    /*"minimum t-label","maximum t-label","maximum boundary t-label",
       "maximum 1-cell b-label".*/

    
	
	
	
    for ( int i = 0; i < number_of_relators ; i++)
    {

		/* First assign to 'minlb[i]'. */
		mintlb[i] = tlabel;
		for (int j = 0; j < length_of_attaching_map[i]; j++)
		{
	    	int one_cell = attaching_map[i][j];
	    	if ( one_cell > 0 )
	    	{
				Slice_iter<int> teltptr = tlabels_for_1_cell[one_cell-1]->row(marker[one_cell-1]);
				Slice_iter<int> beltptr = blabels_for_1_cell[one_cell-1]->row(marker[one_cell-1]);
								 
				for (int k=0; k< n_tpts_on_1_cell[one_cell-1]; k++)
		    		*teltptr++ = tlabel++;

				for (int k=0; k< 2*n_tpts_on_1_cell[one_cell-1]; k++)
		    		*beltptr++ = blabel++;

				marker[one_cell-1]++;
	    	}
	    	else
	    	{
				one_cell *= -1;
				
				Slice_iter<int> teltptr = tlabels_for_1_cell[one_cell-1]->row(marker[one_cell-1]);
				Slice_iter<int> beltptr = blabels_for_1_cell[one_cell-1]->row(marker[one_cell-1]);

				/* new code to allocate labels in reverse order 18/2/06 */
				tlabel += n_tpts_on_1_cell[one_cell-1]-1;
				
				for (int k=0; k< n_tpts_on_1_cell[one_cell-1]; k++)
	    			*teltptr++ = tlabel--;
				
				tlabel += n_tpts_on_1_cell[one_cell-1]+1;
				blabel += 2*n_tpts_on_1_cell[one_cell-1]-1;
				
				for (int k=0; k< 2*n_tpts_on_1_cell[one_cell-1]; k++)
	    			*beltptr++ = blabel--;
				
				blabel += 2*n_tpts_on_1_cell[one_cell-1]+1;
				marker[one_cell-1]++;
	    	}
		}

		/* Now assign the maximum t-label on the boundary of this disc to 
		'mbtl[i]' and the maximum 1-cell b-label to 'm1cbl[i]'. */
		mbtl[i] = tlabel-1;
		m1cbl[i] = blabel-1;

		/* Assign t-labels to the interior 1-cells */
		for (int j=0; j< length_of_attaching_map[i]-3; j++)
		{
			Slice_iter<int> teltptr = tlabels_for_1_cell[interior_1_cell]->row(0);

	    	for( int k =0; k< n_tpts_on_1_cell[interior_1_cell]; k++)
			*teltptr++ = tlabel++;
			
	    	interior_1_cell++;
		}

		/* Now assign 'maxtlb[i]'. */
		maxtlb[i] = tlabel-1;

    }

    if (decomp_control::DECOMPOSE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
    {
		debug << "decompose: mintlb = ";
		for (int i=0; i< number_of_relators; i++)
			debug << setw(4) <<mintlb[i] << ' ';
		debug << endl;
			
		debug << "decompose: maxtlb = ";
		for (int i=0; i< number_of_relators; i++)
			debug << setw(4) <<maxtlb[i] << ' ';
		debug << endl;

		debug << "decompose: mbtl =   ";
		for (int i=0; i< number_of_relators; i++)
			debug << setw(4) <<mbtl[i] << ' ';
		debug << endl;

		debug << "decompose: m1cbl =  ";
		for (int i=0; i< number_of_relators; i++)
			debug << setw(4) <<m1cbl[i] << ' ';
		debug << endl;

		debug << "decompose: t-label matracies" << endl;
		for (int i = 0; i< number_of_1_cells; i++)
		{
	    	debug << "decompose: one cell " << i+1 << endl;

			if (n_tpts_on_1_cell[i])
				print(*tlabels_for_1_cell[i], debug, 4, "decompose:   ");
			else
				debug << "decompose:   void intersection with track" << endl;
		}

		debug << "decompose: b-label matracies" << endl;
		for (int i = 0; i< number_of_1_cells_in_rosette; i++)
		{
	    	debug << "decompose: one cell " << i+1 << endl;

			if (n_tpts_on_1_cell[i])
				print(*blabels_for_1_cell[i],debug,4, "decompose: ");
			else
				debug << "decompose:   void intersection with track" << endl;
		}
    }


	
    /**********************************************************************
    It remains to assign the 0-cell labels.  Since the #(0-cells)<=#(1-cells 
    in the rosette) we may utilize 'marker' for this purpose.  We shall note 
    the min/max 0-cell blabels attached to each disc, recording the
    information in  "extocbl".
    
    If there are added 0-cells to consider, ie. if there are 2-torsion 
	generators given, we require the number of occurrences of these added 
	0-cells in the relator discs.  However such a 0-cell lies at the mid point 
	of a petal of R(L) corresponding to a 2-torsion generator and by our insistance 
	that 2-torsion generators and relators be entered first, and in the same
    order, we may use 'occ of order 2 gen' to determine the number of  occurrences.
    ***********************************************************************/
    vector<int>* blabels_for_0_cell[number_of_0_cells];

    blabels_for_0_cell[0] = new vector<int> (occ_of_cv);   

    for (int i=1; i<number_of_0_cells; i++)
		blabels_for_0_cell[i] = new vector<int> (occ_of_order_2_gen[i-1]);

    matrix<int> ext0cbl(number_of_relators,2);
			   /* "extreme 0-cell b-labels[i]" will hold (min,max) 
			       0-cell blabels for the ith disc. */

    for (int i=0; i< number_of_0_cells; i++)
		marker[i] = 0;


    /********************************************************************
    For future use, we shall need to know, for each 1-cell in R(L) that has 
    void intersection with t, the 0-cell b-labels attached to the boundary of 
    each occurrence of that 1-cell in the boundary of some relator disc.  For 
    each such 1-cell, i, in R(L) we set up a 'occurrences of 1 cell[i]'x2 
    matrix to hold these boundary labels, referenced by "b0cbl for 1 cell", 
    "boundary 0-cell b-labels for 1-cell".
    **********************************************************************/
    matrix<int>* b0cbl_for_1_cell[number_of_1_cells_in_rosette];

    for (int i=0; i<number_of_1_cells_in_rosette; i++)
	{
		if ( n_tpts_on_1_cell[i] == 0 )
	    	b0cbl_for_1_cell[i] = new matrix<int> (occurrences_of_1_cell[i],2);
		else
			b0cbl_for_1_cell[i]	=0;
	}
	

    for (int i = 0 ; i < number_of_relators; i++)
    {
		if (relator[i].length() == 2)
		{
	    	/* There are exactly 4 0-cells on the disc. */

	    	/* 1st 0-cell. */
			vector<int>& blabels0 = *blabels_for_0_cell[0];
	    	blabels0[marker[0]] = blabel;

	    	/*This is the min 0-cell blabel on the disk.*/
	    	ext0cbl[i][0] = blabel;

	    	blabel++;
	    	marker[0]++;

	    	int map_elt = abs(attaching_map[i][0])-1;
	    	if (n_tpts_on_1_cell[map_elt] == 0)
	    	{
				matrix<int>& b0cbl = *b0cbl_for_1_cell[map_elt];
				b0cbl[bmarker[map_elt]][0] = blabel-1;
				b0cbl[bmarker[map_elt]][1] = blabel;
				bmarker[map_elt]++;
	    	}

	    	/* 2nd 0-cell. */

	    	/* By our requirements on the input information to SOLVE, we must
	    	   be dealing with the (i+1)st 0-cell. */
			vector<int>& blabels2 = *blabels_for_0_cell[i+1]; //blabels for 0-cell in 2-torsion generator
	    	blabels2[marker[i+1]] = blabel;

	    	blabel++;
	    	marker[i+1]++;

	    	map_elt = abs(attaching_map[i][1])-1;
	    	if (n_tpts_on_1_cell[map_elt] == 0)
	    	{
				matrix<int>& b0cbl = *b0cbl_for_1_cell[map_elt];
				b0cbl[bmarker[map_elt]][0] = blabel-1;
				b0cbl[bmarker[map_elt]][1] = blabel;
				bmarker[map_elt]++;
	    	}

	    	/* 3rd 0-cell. */
	    	blabels0[marker[0]] = blabel;

	    	blabel++;
	    	marker[0]++;

	    	map_elt = abs(attaching_map[i][2])-1;
	    	if (n_tpts_on_1_cell[map_elt] == 0)
	    	{
				matrix<int>& b0cbl = *b0cbl_for_1_cell[map_elt];
				b0cbl[bmarker[map_elt]][0] = blabel-1;
				b0cbl[bmarker[map_elt]][1] = blabel;
				bmarker[map_elt]++;
	    	}


	    	/* 4th 0-cell. */
	    	blabels2[marker[i+1]] = blabel;

	    	/* This is the max 0-cell blabel on the disk. */
	    	ext0cbl[i][1] = blabel;

	    	blabel++;
	    	marker[i+1]++;

	    	map_elt = abs(attaching_map[i][3])-1;
	    	if (n_tpts_on_1_cell[map_elt] == 0)
	    	{
				matrix<int>& b0cbl = *b0cbl_for_1_cell[map_elt];
				b0cbl[bmarker[map_elt]][0] = blabel-4;
				b0cbl[bmarker[map_elt]][1] = blabel-1;
				bmarker[map_elt]++;
	    	}
		}   
		else
		{
			vector<int>& blabels0 = *blabels_for_0_cell[0];
	    	blabels0[marker[0]] = blabel;

	    	/*This is the min 0-cell blabel on the disk.*/
	    	ext0cbl[i][0] = blabel;

	    	blabel++;
	    	marker[0]++;

	    	int map_elt = abs(attaching_map[i][0])-1;
	    	if (n_tpts_on_1_cell[map_elt] == 0)
	    	{
				matrix<int>& b0cbl = *b0cbl_for_1_cell[map_elt];
				b0cbl[bmarker[map_elt]][0] = blabel-1;
				b0cbl[bmarker[map_elt]][1] = blabel;
				bmarker[map_elt]++;
	    	}

	    	for (int j = 1; j < length_of_attaching_map[i] ; j++)
	    	{
				char gen = generator_corresponding_to_1_cell[abs(attaching_map[i][j])-1];

				if ( gen == generator_corresponding_to_1_cell[abs(attaching_map[i][j-1])-1]
				     && find_generator(generators_of_order_2, gen) != NOT_FOUND
				   )
				{
		    		/**************************************************
		    		Since 2-torsion generators are always written in the
		    		positive, the j th 0-cell on the i th disc must be
		    		the (chptr-generators_of_order_2)+1 st 0-cell.
		    		****************************************************/
		    		string::size_type a0c = find_generator(generators_of_order_2, gen) + 1; // "added 0 cell"
					vector<int>& blabels2 = *blabels_for_0_cell[a0c];

	    			blabels2[marker[a0c]] = blabel;
				
		    		blabel++;
		    		marker[a0c]++;
				}
				else
				{ 
	    			blabels0[marker[0]] = blabel;

		    		blabel++;
		    		marker[0]++;
				}

				map_elt = abs(attaching_map[i][j])-1;
				if (n_tpts_on_1_cell[map_elt] == 0)
				{
		    		if (j == length_of_attaching_map[i]-1)
					{
						matrix<int>& b0cbl = *b0cbl_for_1_cell[map_elt];
						b0cbl[bmarker[map_elt]][0] = blabel - length_of_attaching_map[i];
						b0cbl[bmarker[map_elt]][1] = blabel-1;
		    		}
					else
					{
						matrix<int>& b0cbl = *b0cbl_for_1_cell[map_elt];
						b0cbl[bmarker[map_elt]][0] = blabel-1;
						b0cbl[bmarker[map_elt]][1] = blabel;
		    		}

					bmarker[map_elt]++;

				}

				if( j == length_of_attaching_map[i] - 1 )
		    		ext0cbl[i][1] = blabel-1;
	    	}
		}
    }

if (decomp_control::DECOMPOSE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "decompose: blabels for 0-cells:" << endl;
	
    debug << "decompose:   0-cell " << 1 << "\ndecompose:   ";
	vector<int> blabels = *blabels_for_0_cell[0];
	for (int j=0; j< occ_of_cv; j++)
		debug << setw(4) << blabels[j];
	debug << endl;
		
	for (int i = 1; i< number_of_0_cells; i++)
	{
	    debug << "decompose:   0-cell " << i+1 << "\ndecompose:   ";		
		blabels = *blabels_for_0_cell[i];
		for (int j=0; j< occ_of_order_2_gen[i-1]; j++)
			debug << setw(4) << blabels[j];			
	    debug << endl;
	}

	debug << "decompose: boundary 0-cell b-labels" << endl;
	for (int i = 0; i< number_of_1_cells_in_rosette; i++)
	{
	    if ( n_tpts_on_1_cell[i] == 0 )
	    {
			debug << "decompose:   1-cell " << i+1 << endl;
			print (*b0cbl_for_1_cell[i],debug, 0, "decompose:     ");
//			debug << *b0cbl_for_1_cell[i] << endl;
	    }
	    else
			debug << "decompose:   1-cell " << i+1 << " intersects track" << endl;
	}

	debug << "decompose: extreme 0cbl:" << endl;
	print (ext0cbl,debug,0,"decompose:   ");
//	debug << ext0cbl << endl;
}

    /*********************************************************************
    NOTICE: For any non empty matrix referenced by 'b0cbl for 1 cell', for 
    each row, the minimum b-label of that row occurs in the first column.
    **********************************************************************/

    /* At this point we calculate, for each disc Di, the total number of 
      t-labels in the interior of the discs D1,...,Di-1, again for future use. */
    int citls[number_of_relators]; /* "cumulative interior t-label sum" */
   citls[0] = 0;

   for (int i=1; i< number_of_relators; i++)
       citls[i] = citls[i-1] + maxtlb[i-1] - mbtl[i-1];

    /***********************************************************************
    We next determine which t-labels are joined by an edge of t (regarding t 
    as a graph with Vt the intersection of t with the 1-skeleton of L) on 
    each relator disc; we store the information as pairs, forming the rows of 
    a matrix.  We have one pair for each edge of t.

         Given a variable, and the two sets of t-labels the edges of t 
    determined by that variable "join", we make the pairing by moving "into
    the vertex";  ie if 'track[variable]'=x we must have x t-labels in each
    list paired by the corresponding edges of t, thus the pairings must be of
    the x t-labels closest to that vertex on the relator disc.

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

    int partial_sum[number_of_relators];

    /* 'partial sum' will indicate how many 
       pairs arise from each relator disc. */

    for (int i = 0 ; i < number_of_relators ; i++)
	partial_sum[i] = 0;
   
    /* Assign partial sums. */
    int variable = 0;

    for ( int i = 0; i < number_of_relators ; i++)
	{
		for (int j = 0; j < 3*(length_of_attaching_map[i]-2) ; j++)
		    partial_sum[i] += track[variable++];
    }
	 
    matrix<int>* pair[number_of_relators];

if (decomp_control::DECOMPOSE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "decompose: partial_sum = ";
	for (int i=0; i<number_of_relators; i++)
		debug << partial_sum[i] << ' ';
	debug << endl;
}

    /* Declare space for 'pair'. */
    for (int i=0; i<number_of_relators; i++)
		pair[i] = new matrix<int> (partial_sum[i], 3);   

    /****************************************************************
    We traverse the relator discs in the canonical order, thereby moving
    through 'track'.  We make use of 'marker' to determine which set of
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
	    	int tv = track[variable];
			      
//	    	matrix<int>& matrixptr = *tlabels_for_1_cell[one_cell-1]; 

	    	Slice_iter<int> r1 = tlabels_for_1_cell[one_cell-1]->row(mark);
	    	Slice_iter<int> r2 = tlabels_for_1_cell[one_cell-1]->row(mark+1);
	    	Slice_iter<int> r3 = tlabels_for_1_cell[one_cell-1]->row(mark+2);

	    	matrix<int>& matrixptr = *pair[i]; /* matrixptr now points at the target
						 of the asignation */
	    	/* Place 1. */
	    	r1 += tv-1;
	    	r3 += tv;

	    	for (int j=0; j< tv; j++)
			{
				matrixptr[pair_marker][0] = *r1--;
				matrixptr[pair_marker][1] = *r3++;
				matrixptr[pair_marker][2] = variable;
				pair_marker++;
			}

           	variable++;

	    	/* Place 2. */
	    	r1 += tv+1;
	    	r2 += tv-1;

	    	for (int j=0; j< tv; j++)
			{
				matrixptr[pair_marker][0] = *r1++;
				matrixptr[pair_marker][1] = *r2--;
				matrixptr[pair_marker][2] = variable;
				pair_marker++;
			}

	    	variable++;

	    	/* Place 3. */
	    	r2 += tv+1;
	    	r3 -= tv+1; /* i.e. the tv th label in the row */

	    	for (int j=0; j< tv; j++)
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

	    	Slice_iter<int> r1 = tlabels_for_1_cell[one_cell1-1]->row(marker[one_cell1-1]);

			Slice_iter<int> r2 = r1; // assignation just for the declaration, since there's no default constructor 

	    	if (one_cell1 == one_cell2)
				r2 = tlabels_for_1_cell[one_cell1-1]->row(marker[one_cell1-1]+1); 						
	    	else
				r2 = tlabels_for_1_cell[one_cell2-1]->row(marker[one_cell2-1]); 

	    	Slice_iter<int> r3 = tlabels_for_1_cell[interior_1_cell]->row(0); 

	    	int n1 = tlabels_for_1_cell[one_cell1-1]->numcols();
	    	int n2 = tlabels_for_1_cell[one_cell2-1]->numcols();
	    	int n3 = tlabels_for_1_cell[interior_1_cell]->numcols();    

	    	/* Place 1. */
	    	int tv = track[variable];      

	    	if ( s1 ) r1 += tv-1; else r1 += n1-tv;
	    	r3 += tv-1;

	    	for (int j=0; j< tv; j++)
	    	{

				matrixptr[pair_marker][0] = *r1;
				matrixptr[pair_marker][1] = *r3--;
				matrixptr[pair_marker][2] = variable;
				pair_marker++;
				
				if ( s1 ) r1--; else r1++;
	    	}
	    	variable++;

	    	/* Place 2. */
	    	tv = track[variable];      

	    	if ( s1 ) r1 += n1-tv+1; else r1 -= n1-tv+1;
	    	if ( s2 ) r2 += tv-1; else r2 += n2-tv;

	    	for (int j=0; j< tv; j++)
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
	    	tv = track[variable];      

	    	if ( s2 ) r2 += n2-tv+1; else r2 -= n2-tv+1;
	    	r3 += n3-tv+1; 

	    	for (int j=0; j< tv; j++)
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
	    	/* Notice that if one_cell1=one_cell2 then this code is still correct. */

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

				Slice_iter<int> r1 = tlabels_for_1_cell[interior_1_cell]->row(0); 
				Slice_iter<int> r2 = tlabels_for_1_cell[in_rosette-1]->row(marker[in_rosette-1]);
				Slice_iter<int> r3 = tlabels_for_1_cell[interior_1_cell+1]->row(0); 
				
				int n1 = tlabels_for_1_cell[interior_1_cell]->numcols();
				int n2 = tlabels_for_1_cell[in_rosette-1]->numcols();
				int n3 = tlabels_for_1_cell[interior_1_cell+1]->numcols();
				
				/* Place 1. */
				tv = track[variable];      
				r1 += tv-1;
				r3 += tv-1;
				
				for (int k=0; k< tv; k++)
				{
					matrixptr[pair_marker][0] = *r1--;
					matrixptr[pair_marker][1] = *r3--;
					matrixptr[pair_marker][2] = variable;
					pair_marker++;
				}
					
				variable++;

				/* Place 2. */
				tv = track[variable];      
				r1 += n1-tv+1;
				if ( s2 ) r2 += tv-1; else r2 += n2-tv;

				for (int k=0; k< tv; k++)		
				{
					matrixptr[pair_marker][0] = *r1++;
					matrixptr[pair_marker][1] = *r2;
					matrixptr[pair_marker][2] = variable;
					pair_marker++;

		    		if (s2) r2--; else r2++;
				}
				variable++;

				/* Place 3. */
				tv = track[variable];      
				if (s2) r2 += n2-tv+1; else r2 -= n2-tv+1;
				r3 += n3-tv+1;

				for (int k=0; k< tv; k++)		
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

	    	r1 = tlabels_for_1_cell[interior_1_cell]->row(0); 
	    	r2 = tlabels_for_1_cell[one_cell1-1]->row(marker[one_cell1-1]); 
		
	    	if (one_cell1 == one_cell2)
				r3 = tlabels_for_1_cell[one_cell1-1]->row(marker[one_cell1-1]+1); 
	    	else
				r3 = tlabels_for_1_cell[one_cell2-1]->row(marker[one_cell2-1]); 

           	n1 = tlabels_for_1_cell[interior_1_cell]->numcols();
	    	n2 = tlabels_for_1_cell[one_cell1-1]->numcols();
	    	n3 = tlabels_for_1_cell[one_cell2-1]->numcols();
			
	    	/* Place 1. */
	    	tv = track[variable];          
	    	r1 += tv-1;
	    	if ( s3 ) r3 += n3-tv; else r3 += tv-1;

	    	for (int j=0; j< tv; j++)
	    	{
				matrixptr[pair_marker][0] = *r1--;
				matrixptr[pair_marker][1] = *r3;
				matrixptr[pair_marker][2] = variable;
				pair_marker++;

				if ( s3 ) r3++; else r3--;
	    	}
	    	variable++;
			
	    	/* Place 2. */
	    	tv = track[variable];             
	    	r1 += n1-tv+1;
	    	if (s2) r2 += tv-1; else r2 += n2-tv;

	    	for (int j=0; j< tv; j++)
	    	{
				matrixptr[pair_marker][0] = *r1++;
				matrixptr[pair_marker][1] = *r2;
				matrixptr[pair_marker][2] = variable;
				pair_marker++;

				if ( s2 ) r2--; else r2++;
	    	}
	    	variable++;


	    	/* Place 3. */
	    	tv = track[variable];             
	    	if (s2) r2 += n2-tv+1; else r2 -= n2-tv+1;
	    	if (s3) r3 -= n3-tv+1; else r3 += n3-tv+1;

	    	for (int j=0; j< tv; j++)
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


if (decomp_control::DECOMPOSE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "decompose: t-label pair matracies" << endl;
	for (int i=0; i< number_of_relators; i++)
	{
	    debug << "decompose:   relator " << i+1 << endl;
	    if (partial_sum[i])
			print(*pair[i],debug,0,"decompose:     ");
		else
			debug << "decompose:     void intersection with track" << endl;
	}
}


    /*******************************************************************
    We now determine the component arcs of t on each relator disc.
   
    Declare space for the disc permutations.  For each disc, i, we set up
    a row with maxtlb[i]-mintlb[i] elements.  Note that each label occurs
    in "pair[i]" exactly once iff the label is on the boundary of the disc 
    and exactly twice iff it is on an interior 1-cell.
    *********************************************************************/

    vector<int>* dperm[number_of_relators]; /*"disc permutations" */
    for( int i =0; i< number_of_relators ; i++)
		dperm[i] = new vector<int> (maxtlb[i] - mintlb[i] + 1);

    /*********************************************************************
    Assign disc permutations.  By moving successively through the rows of 
    "pair[i]" we deal with the "variable corners" in the
    canonical order and, by nature of our construction of "pair[i]", have the 
    following properties satisfied by each row:
   
    1.  If a t-label in the first column is <='mbtl[i]', it is the initial
	t-label of a component arc on the disc.  If it is
	>'mbtl[i]', the pair in that row determines an edge of t whose
	initial point lies on an interior 1-cell.  In this case the t-label
	at this interior point will have occured exactly once, in the
	second column of a row of "pair[i]", BEFORE the current pair.
   
    2.  If a t-label in the second column is <='mbtl[i]', it is the terminal 
	point of a component arc on the disc.  If it is >'mbtl[i]', the 
	pair in that row determines an edge of t whose terminal point
	lies on an interior 1-cell.  In this case the t-label
	at this terminal point will occur exactly once, in the first
	column of a row of "pair[i]", AFTER the current pair.
   
        The disc permutation will allow us to find the component arcs.  A
    t-label in the first column of "pair[i]" determines the initial point of
    some edge of P and therefore shall be made to index the t-label at the
    terminal point of that edge in the permutation.  If a t-label in the second
    column is <='mbl[i]' it terminates a component arc and, so, shall be made
    to index zero, to signify this.  If the t-label in the second column is
    >'mbtl[i]' it will occur, also, in the first column in a subsequent row and
    so no assignation will be made at this first encounter.  By properties 1
    and 2 above, it can be seen that this procedure completely determines the 
    permutation.
    *********************************************************************/

    for (int i=0;i<number_of_relators; i++)
    {
		matrix<int>& pairs = *pair[i];
		vector<int>& perm = *dperm[i];
		int min = mintlb[i];
		int bdry_max = mbtl[i];

		for ( int j = 0; j< partial_sum[i]; j++)
		{
		    perm[pairs[j][0] - min] = pairs[j][1] - min;

		    if ( pairs[j][1] <= bdry_max )
				perm[pairs[j][1] - min] = 0;
		}
    }

if (decomp_control::DECOMPOSE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "decompose: disc t-label permutations\n";
	for (int i=0; i< number_of_relators; i++)
	{
	    debug << "decompose:   relator " << i+1 << ": ";
		vector<int>& perm = *dperm[i];

		for (unsigned int j=0; j < perm.size(); j++)
			debug << perm[j] << ' ';
	    debug << endl;
	}
}

     
    /********************************************************************
    Evaluate the number of components on each disc.  This is simply a matter 
    of counting the terminating points of arcs on the disc, which are
    identified by a unique zero amongst the first 'mbtl[i]' elements of the 
    corresponding disc permutation.
    *************************************************************/
    int ndtc[number_of_relators]; /* "number of disc track-components" */
    for ( int i=0; i< number_of_relators ;i++)
		ndtc[i] = 0;

    for ( int i =0; i< number_of_relators; i++)
    {
		vector<int>& perm = *dperm[i];
		for ( int j = 0; j < mbtl[i] - mintlb[i] + 1 ; j++)
		{
		    if ( perm[j] == 0 )
				ndtc[i] += 1;
		}
    }
   
  
    /***************************************************************
    Declare space for disc arcs.  The maximum number of edges an arc may 
    contain is the number, N, of 2-cells in the disc under consideration. 
    The number of labels on an arc is thus <=N+1 and, by storing the 
    terminating zero, we must record a maximum of N+2 integers for each arc.  
	These integers will be held as the rows of a matrix.
    *******************************************************************/
   
    matrix<int>* darcs[number_of_relators]; /*"disc arcs"*/
    for (int i=0; i < number_of_relators ; i++)
		darcs[i] = new matrix<int>(ndtc[i], length_of_attaching_map[i]);
		
    /********************************************************************
    Assign disc arcs.  The initial and terminal labels of each arc may be
    recovered from the indexes of the first 'mbl[i]' components of the 
    corresponding permutation (this number must therefore be even).  
    The terminal labels are indicated by zero and so the labels on each arc 
	may be obtained by "tracing through" the permutation until a zero is 
	encountered.  This zero is stored to indicate the termination of that arc.
    *******************************************************************/
    
	for (int i=0; i< number_of_relators; i++)
    {
		vector<int>& perm = *dperm[i];
		matrix<int>& arcs = *darcs[i];
		int min = mintlb[i];
		int arcn = 0; /* arc number */

		for ( int j  = 0 ; j < mbtl[i] - mintlb[i] + 1; j++)
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

if (decomp_control::DECOMPOSE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "decompose: number of disc components of track: ";
	for (int i=0; i< number_of_relators; i++)
		debug << ndtc[i] << ' ';
	debug << endl;

	debug << "decompose: disc t-arcs" << endl;
	for (int i=0; i< number_of_relators; i++)
	{		
	    debug << "decompose:   relator " << i+1 << ":" << endl;
	    if (partial_sum[i])
	    	print(*darcs[i],debug,0,"decompose:     ");
	    else
			debug << "decompose:     void intersection with track" << endl;
	}
}
     
    /********************************************************************
    Next we determine how the disc arc components of t cross the relator 
    discs.  To accomplish this we utilize a matrix whose rows are made up 
    of the end labels of each disc component.  This matrix will also show 
    us how theses components become connected when attached to the rosette.
    ********************************************************************/

    int tndtc = 0; //"total number of disc track components" 
    for ( int i =0; i<number_of_relators ; i++)
		tndtc += ndtc[i];

    matrix<int> tael(tndtc,2); /* "track arc end labels" */
    int elm = 0; /* "end label marker" */

    for ( int i =0; i< number_of_relators ; i++)
    {
		matrix<int>& arcs = *darcs[i];
		for (int j=0; j< ndtc[i]; j++)
		{
		    tael[elm][0] = arcs[j][0];
		    int k = 1;
		    do  
				k++; 
		    while (arcs[j][k] != 0);
			
		    /* We know there are at least two labels on each arc.  When 
		    this loop terminates, k indicates the terminating zero of the
		    current  arc. */

		   tael[elm][1] = arcs[j][k-1];
		   elm++;
		}
    }


if (decomp_control::DECOMPOSE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "decompose: track arc end labels" << endl;
	print (tael,debug, 0, "decompose:   ");
}

    /***************************************************************
    'tael' and the t-label matrices give us sufficient information to 
    determine a maximal tree in t, regarding t as a graph with Vt the  
    intersection of t with R(L).  We now turn to the b-labels; however, 
    we are now dealing with a graph, X, having some edges in the boundary 
    of the relator discs and so for each b-label, b, there is more than 
    one other label on the same disc that is joined to b by some edge in EX.
    In fact, it is not difficult to see that there are exactly two other 
    such labels but their position on the disc, in relation to b, may be 
    in one of three places.  Suppose b1 is a b-label, joined to b on some 
    disc by an oriented edge e in EX, then one of the following hold.

    0.  The orientation of e is clockwise in the boundary of the disc;
    1.  e has non empty intersection with the interior of the disc;
    2.  The orientation of e is anticlockwise in the boundary of the disc.

    Moreover, for each b-label, b, in the disc D, there is a unique b-label in 
    D satisfying one of 0,1 or 2, for exactly two of the cases 0,1,2.  If b 
    labels a 0-cell then cases 1 and 3 are satisfied and otherwise case 2 and 
    one of cases 1 and 3 are satisfied.  We therefore set out, for each 
    b-label, a row of 3 integers, stored in a matrix, and shall assign the 
    b-label satisfying case i (i=1,2,3) in position i, if case i is satisfied and shall store zero 
    in position i if case i is not satisfied.
    ***********************************************************************/
    blabel--; /*'blabel now denotes the number of b-labels used.*/
    matrix<int> jbba(blabel,3); /* "joined by a b-arc" */

    /********************************************************************
    By construction, the first 'm1cbl[number of relators]' b-labels are all 
    joined to a b-label satisfying 2 above.  These joining arcs are the disc 
    arc components of the boundary of the band B so we may use 'jbta',
    'm1cbl',  and 'citls' to determine b1.
    **********************************************************************/

    /* We make the following declaration for later use. */
    int max1cbl = m1cbl[number_of_relators-1]; /* "maximum 1-cell b-label"*/

    int bl = 1;

    for (int i=0; i< number_of_relators; i++)
    {
		int mark = m1cbl[i];
		do
		{
	    	if (bl > mark)
			break;

			int tlabel;
			
	    	/* determine the t-label adjacent to bl */
	    	if ( bl % 2 == 1)
				tlabel = (bl + 1)/2 + citls[i]; 
	    	else
				tlabel = bl/2 + citls[i];

	    	/* now determine the t-label at the other 
	    	   end of the disc t-arc from tl */
			for (int j=0; j< tndtc; j++)
			{
				for (int k=0; k< 2; k++)
				{
					if (tael[j][k] == tlabel)
					{
			    		/* reassign tlabel to be the other 
			    		   end of this disc t-arc */
			    		if ( k == 0 )
							tlabel = tael[j][1];
		    			else
							tlabel = tael[j][0];
						
						goto tlabel_reassigned;
					}
				}
			}
			tlabel_reassigned:
			
	    	/* tlabel now indicates the other end of the t-arc,
	    	   we assign to the bl-1th row of jbba since the array
	    	   is indexed from 0.
	    	*/
	    	if (bl % 2 == 1)
				jbba[bl-1][1] = 2 * (tlabel - citls[i]);
	    	else
				jbba[bl-1][1] = 2 * (tlabel - citls[i]) -1;
	    	bl++;
		} while (bl <= mark);
    }

    /* No other b-label is joined to a b-label satisfying 1, above, so we 
       may complete the assignations to 'jbba[,1]'. */
    for (int i=bl; i< blabel; i++)
		jbba[i][1] = 0;

    /***********************************************************************
    To assign b-labels to 'jbba[,0]' and 'jbba[,2]' we notice that, if b is a 
    b-label in the interior if a 1-cell, if b is odd there is a b-label b1 
    satisfying 1 above and if b is even there is a b-label satisfying 3
    above.  The first 'max1cbl' b-labels are of this type.
    **********************************************************************/
    int bl0 = max1cbl + 1; /* "b-label on 0-cell" now denotes the first b-label on
			 a 0-cell.*/
    bl = 1;
    for (int i=0; i< number_of_relators; i++)
    {
		/* Deal with the first 'length of attaching map[i]'-1 
		   1-cells in the boundary. */
		for (int j=0; j< length_of_attaching_map[i] - 1; j++)
		{
	    	int map_elt = abs(attaching_map[i][j]);
	    	int number_of_points=n_tpts_on_1_cell[map_elt - 1];

	    	if (number_of_points == 0)
			{
				bl0++;
			}
	    	else
	    	{
				jbba[bl-1][0] = bl0++;
				jbba[bl-1][2] = 0;
				bl++;

				for (int k=0; k< 2*number_of_points-2; k++)
				{
		    		if ( bl%2 == 1)
		    		{
						jbba[bl-1][0] = bl-1;
						jbba[bl-1][2] = 0;
						bl++;
		    		}
		    		else
		    		{
						jbba[bl-1][0] = 0;
						jbba[bl-1][2] = bl+1;
						bl++;
		    		}
				}
				jbba[bl-1][0] = 0;
				jbba[bl-1][2] = bl0;
				bl++;
	    	}
		}

		/* Now deal with the last 1-cell in the boundary. */
		int map_elt = abs(attaching_map[i][length_of_attaching_map[i]-1]);
		int number_of_points = n_tpts_on_1_cell[map_elt-1];

		if( number_of_points != 0)
		{
	    	jbba[bl-1][0] = bl0++;
	    	jbba[bl-1][2] = 0;
	    	bl++;

	    	for (int k=0; k< 2*number_of_points-2; k++)
	    	{
				if ( bl%2 == 1)
				{
		    		jbba[bl-1][0] = bl-1;
		    		jbba[bl-1][2] = 0;
		    		bl++;
				}
				else
				{
		    		jbba[bl-1][0] = 0;
		    		jbba[bl-1][2] = bl+1;
		    		bl++;
				}
	    	}
	    	jbba[bl-1][0] = 0;
	    	jbba[bl-1][2] = bl0 - length_of_attaching_map[i];
	    	bl++;
		}    
		else
	    	bl0++;
    }

    /* Finally, we deal with the b-labels on the 0-cells.  
       We re-initialize 'bl0' and 'bl'. */
    bl0 = bl; /*'bl0' now contains the first b-label on a 0-cell. */
    bl = 1;

    for (int i =0; i< number_of_relators; i++)
    {
		/* Deal with the first 0-cell of the disc. */
		bool no_tpoints_yet = true;
		if (n_tpts_on_1_cell[abs(attaching_map[i][length_of_attaching_map[i]-1])-1] == 0)
	    	jbba[bl0-1][0] = bl0+length_of_attaching_map[i] - 1;
		else
	    	jbba[bl0-1][0] = m1cbl[i];

		if (n_tpts_on_1_cell[abs(attaching_map[i][0])-1] == 0)
	    	jbba[bl0-1][2] = bl0 + 1;
		else
	    	jbba[bl0-1][2] = bl;


		bl0++;

		if (n_tpts_on_1_cell[abs(attaching_map[i][0])-1] != 0)
		{
	    	bl+= 2 * n_tpts_on_1_cell[abs(attaching_map[i][0])-1]-1;
	    	no_tpoints_yet = false;
		}

		/* Now the next 'length of attaching map[i]'-2 0-cells. */
		for (int j = 1; j< length_of_attaching_map[i] - 1; j++)
		{
	    	if (n_tpts_on_1_cell[abs(attaching_map[i][j-1])-1] == 0)
				jbba[bl0-1][0] = bl0 - 1;
	    	else
				jbba[bl0-1][0] = bl;

	    	if (n_tpts_on_1_cell[abs(attaching_map[i][j])-1] == 0)
				jbba[bl0-1][2] = bl0 + 1;
	    	else if (no_tpoints_yet)
				jbba[bl0-1][2] = bl;
	    	else
				jbba[bl0-1][2] = bl+1;

	    	bl0++;
	    	if (no_tpoints_yet && n_tpts_on_1_cell[abs(attaching_map[i][j])-1] != 0)
	    	{
				no_tpoints_yet = false;
				bl+= 2 * n_tpts_on_1_cell[abs(attaching_map[i][j])-1]-1;
	    	}
	    	else
				bl+= 2 * n_tpts_on_1_cell[abs(attaching_map[i][j])-1];
		}

		/* Finally, the last 0-cell. */
		if (n_tpts_on_1_cell[abs(attaching_map[i][length_of_attaching_map[i]-2])-1] == 0)
	    	jbba[bl0-1][0] = bl0 - 1;
		else
	    	jbba[bl0-1][0] = bl;

		if (n_tpts_on_1_cell[abs(attaching_map[i][length_of_attaching_map[i]-1])-1] == 0)
	    	jbba[bl0-1][2] = bl0 - length_of_attaching_map[i] + 1;
		else
	    	jbba[bl0-1][2] = bl + 1;

		/*Assign to 'bl0' the first b-label on a 0-cell in the next disc.*/
		bl0++;

		/*Assign to 'bl' the first b-label in the interior 
		  of a 1-cell on the next disc. */
    	bl = m1cbl[i] + 1;
    }

if (decomp_control::DECOMPOSE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "decompose: jbba" << endl;
	print (jbba,debug,0,"decompose:   ");
}

    /*********************************************************************
    The assignation to 'jbba' is now complete and we turn to the problem of 
    constructing maximal trees in the "t-graph" and "b-graph", then 
    completing the decomposition
    **********************************************************************/

    /**********************************************************************
    The t-graph.

    We store the t-labels at the ends of each disc arc of t in our maximal 
    tree (regarding t as a graph with Vt the intersection of t with R(L)) as 
    pairs forming the rows of a matrix, "ttrarcs" ("t-tree arcs").  We have 
    'tndtc' as an upper bound for the number of such arcs.  To construct the 
    tree we utilize 'jbta' and the label matrices referenced by 'tlabels for 
    1 cell'.  The labels on interior 1-cells are irrelevant to the
    construction, since we are not considering these points as elements of Vt.
    The labels in the columns of the t-label matrices are identified in R(L) 
    via the various attaching maps, as discussed above.  Given one column, 
    identified to the point p in R(L), we may find an arc joining p to some 
    other point q in R(L), extending the tree, by dereferencing a place in 
    'jbta'.  We shall set components of 'jbta' to be negative to prevent us 
    creating circuits in the t-graph.  Notice also that each t-label in the
    t-label matrices occurs in a unique matrix and column.  We store this 
    position (to control the construction) as a pair (m,c),
    ("matrix","column"),  of integers, forming the rows of a matrix    
    "tposition".  We must visit each vertex in Vt exactly once, ie
    shall need to reference each column of the t-label matrices (not interior
    1-cells) exactly once, thus giving an upper bound for 'tposition'.  We use 
    the first t-label to initialize the tree construction, and a Boolean 
    variable, "trosette", to record whether or not the t-tree consists of 
    just this single vertex.
    ************************************************************************/
    matrix<int> ttrarcs(tndtc,2); /* t-tree arcs */

    /* Evaluate the upper bound for 'tposition' */
    int n_tpts_on_rosette = 0;
    for (int i = 0; i < number_of_1_cells_in_rosette; i++)
		n_tpts_on_rosette += n_tpts_on_1_cell[i];

    matrix<int> tposition(n_tpts_on_rosette,2);

    int pmarker = 0;
    int amarker = 0;
    /* 'new label will be a label from a column in some t-label matrix,
     "label pointer", "position pointer","position marker" and "arc marker"
     will act as control aids. */

    bool tree_complete = false;
    int new_label = 1;

    do
    {
		int one_cell;
		int column;
		
		/* We have to ensure we don't come back to this vertex of the track again, so 
		   we look for set of t-labels that are identified with new_label, these will be
		   a column of labels in the tlabel matrices. */
		for (int i=0; i< number_of_1_cells_in_rosette; i++)
		{
	    	matrix<int>& matrixptr = *tlabels_for_1_cell[i];
	    	for (unsigned int j=0; j < matrixptr.numrows(); j++)
	    	{
				for (unsigned int k=0; k< matrixptr.numcols(); k++)
				{
		    		if ( matrixptr[j][k] == new_label)
		    		{
						tposition[pmarker][0] = i;
						tposition[pmarker][1] = k;
						pmarker++;
						
						one_cell = i;
						column = k;
						goto new_label_found_in_tlabel_matrices;
		    		}
				}
	    	}
		}
		new_label_found_in_tlabel_matrices:
		/************************************************************************				   
		Set those integers in 'tael' that appear in the column just found to be 
		negative, so we do not look for this vertex (column) again.    
		************************************************************************/

    	matrix<int>& matrixptr = *tlabels_for_1_cell[one_cell];
		for (unsigned int i=0; i< matrixptr.numrows(); i++)
		{
	    	int tlabel = matrixptr[i][column];

	    	for (int j = 0; j< tndtc; j++)
	    	{
				for (int k=0; k< 2; k++)
				{
					if (tael[j][k] == tlabel)
					{
				    	tael[j][k] *= -1;
				    	goto label_found_in_tael;
					}
				}
	    	}
	    	label_found_in_tael:;
		}

		/* We now see how we can extend the tree, ie find a new new_label value.  
		   Each vertex of V(t) corresponds to a column in one of the t-label matrices,
		   indicated by a row in tposition.  We execute a breath first search, looking down each 
		   column in turn (indexed by lpointer), the column to be searched being indexed through 
		   tposition by ppointer.  
		*/
    	unsigned int ppointer = 0;
	    int lpointer = 0;
		bool label_found = false;
	
		do
		{
			int one_cell = tposition[ppointer][0]; // a 1-cell in the rosette
			int column = tposition[ppointer][1];
			
	    	lpointer++;
	    	if ( lpointer > occurrences_of_1_cell[one_cell])
	    	/* ie if lpointer goes off the end of its current column. */
	    	{
				ppointer++;
				lpointer = 1;
	    	}

	    	if (ppointer >= tposition.numrows()) 
			{
				tree_complete = true;
			}
			else
			{
				one_cell = tposition[ppointer][0]; // a 1-cell in the rosette
				column = tposition[ppointer][1];
			}

	    	if (!tree_complete)
	    	{

				matrix<int>& tlabels = *(tlabels_for_1_cell[one_cell]);

				int il = tlabels[lpointer - 1][column];

				/*******************************************************
				"initial label" is the label referred to by 'lpointer' 
				in the current column.

				Look for il in tael.  When we find it, we check to see 
				if we already have the labels at the other end of this 
				arc.  If we have not, the tlabel at the other end of the 
				arc is our new label, and we have found another tree arc. 
				*****************************************************/
				for (int i=0; i< tndtc; i++)
				{
					for (int j=0; j< 2; j++)
					{
						if (abs(tael[i][j]) == il )
						{
							if (j == 0 && tael[i][1] > 0)
							{
					    		new_label = tael[i][1];
					    		label_found = true;
								ttrarcs[amarker][0] = il;
								ttrarcs[amarker][1] = new_label;
								amarker++;
								goto new_tree_arc_search_complete; 
							}
							else if (j == 1 && tael[i][0] > 0)
							{
					    		new_label = tael[i][0];
					    		label_found = true;
								ttrarcs[amarker][0] = il;
								ttrarcs[amarker][1] = new_label;
								amarker++;
								goto new_tree_arc_search_complete; 
							}
						}
					}
				}			
	    	}
			new_tree_arc_search_complete:;
		} while(!label_found && !tree_complete);	
    } while(!tree_complete);

    /*****************************************************************
     We note at this point the number of edges that have been assigned 
     to the t-tree.  The boolean trosette records whether the track is a
	 rosette (in which case there are no edges in the t-tree). 
     *****************************************************************/
    int nett = amarker; // number of edges in t-tree
	bool trosette = (nett? false: true);
	

if (decomp_control::DECOMPOSE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "decompose: ttrarcs" << endl;
	print (ttrarcs,debug,0, "decompose:   ");
	debug << "decompose: tposition" << endl;
	print(tposition,debug,0,"decompose:   ");
}


    /********************************************************************
    The b-graph, X.

    Recall that VX consists of intersections of the boundary of the band B
    with R(L) and the 0-cells of L.  The edges EX comprise the boundary arcs 
	of the band B on each relator disc together with edges in the boundary 
	of the relator discs themselves.  These latter edges are of two types, 
	those whose boundary labels include a label attached to an 0-cell of L 
	and those whose boundary labels are both attached to boundary arcs of 
	the band B.
	                                                                         
    Our strategy is as follows.  Choose b-label 1 as the initial point for
    a maximal tree in the component of X containing it.  Extend to a maximal 
    tree in the component of the boundary of B containing 1; we may use the 
    same algorithm as for the t-tree, with jbba[,2] taking the place of 
    tael.  (If the b-label 2 now lies in the tree, t must be twisted.)  Now 
    extend the tree to  as many 0-cells as possible starting from the vertices 
	already added;if t is twisted, all the 0-cells may be included and the 
	construction of a maximal tree in (the connected graph) X will be complete.  
	
	If t is untwisted but non separating, then the b-graph X is again connected and
	it follows that we may extend the tree from one of the 0-cells just added to a 
	0-cell that lies at the end of a boundary-arc of the band B that is not already 
	in the tree.  We therefore check whether this is indeed the case, if so, 
	we extend from  this new 1-cell b-label over the other component of the boundary 
	of B finally adding any remaining 0-cells to complete the construction.  If
    this is not the case the t must be separating and we shall have completed
    the construction of a maximal tree, T, in one component of X.  For the
    other component, notice that the b-label 2 cannot be contained in VT, so
    we may use this label as the initial point for a maximal tree in the
    second component of X.  Notice, also, that, if t separates, the boundary
    of B consists of two copies of t so by choosing 2 as the initial point for
    the second component, the two trees constructed in the boundary of B are 
    parallel in L, hence "parallel" on each disc.  This enables a quick 
    construction of this second "boundary tree" and then, finally, we add any 
    remaining 0-cells to complete the construction.

         We shall require two analogues, "btrarcs1","btrarcs2", of 'ttrarcs' 
    and one, "bposition", of 'tposition'.  As an upper bound for 'bposition'
    we may use '2*UPB tposition' and for 'btrarcs1(2)', '2*UPB ttrarcs+k',
    where k may be calculated from the number of 0-cells and their
    occurrences.

         Notice that the maximal tree(s) in the b-graph cannot consist of a 
    single vertex, since the first 0-cell b-label will be in the first tree,
    and if t intersects R(L) in exactly one point, it must be non
    separating.
    **********************************************************************/

    /* Evaluate k. */
    int dump = 2 * occ_of_cv;
    for (int i = 0; i< number_of_0_cells - 1 ; i++)
	dump += 2 * occ_of_order_2_gen[i];

    /* Declare space for records */
    matrix<int> btrarcs1(2*tndtc+dump,2);
    matrix<int> btrarcs2(2*tndtc+dump,2);

    matrix<int> bposition(2*tposition.numrows(),2);

    for (unsigned int i = 0; i <  bposition.numrows(); i++)
		bposition[i][0] = -1;
    /* The assignation is used as a control aid, 
       since the upper bound is not exact. */

    /********************************************************************
    Notice we have the objects 'new label', 'lpointer', 'ppointer', 'pm',
    'am', 'label not found' and 'tree complete' in range and available for
    use. 

    Initialize for "band extension" from b-label 1.  
    *********************************************************************/
    new_label = 1;
    unsigned int ppointer = 0;
    int lpointer = 0;
    pmarker = 0;
    amarker = 0;

    /* we shall also require the following control aids. */
    bool twisted = false;
    int nebt1 = 0; // number of edges in b-tree1
    int nebt2 = 0; 

   /************************************************************************
   "maximum band arc 1(2)" will denote the row of 'btrarcs1(2)' containing the
    last arc added by 'band extend' at the first (second) call, when dealing
    with the first component of the b-graph.  "maximum 0-cell connection arc"
    and "extension arc" will denote the last row of 'btrarcs1' added at the
    first call to 'add 0 cells' and the row corresponding to the arc extending
    the tree to the second component of the boundary of the band B
    (if it exists).
    ************************************************************************/
    int max_band_arc_2 = 0;
    int extension_arc = 0;

	band_extend(new_label, number_of_1_cells_in_rosette, occurrences_of_1_cell, blabels_for_1_cell, 
                lpointer, bposition, ppointer, pmarker, jbba, btrarcs1, amarker, twisted);
				 

    /* Assign 'max band arc 1'. */
    int max_band_arc_1 = amarker;

	add_0_cells(number_of_0_cells, blabels_for_0_cell, jbba, amarker, btrarcs1);
    

    /* Assign 'max 0 cell arc' */
    int max_0_cell_arc = amarker;

    /**********************************************************************
    If 'twisted=true' then 2 is in the current tree and all the 0-cells will
    have been added.  We assign the number of edges in the first b-tree
    (currently) to 'nebt1' to deal with the possibility that t separates or 
    is twisted.
    ***********************************************************************/
    nebt1 = amarker;

    /************************************************************************
    By default, we always evaluate the decomposition given by 2t if t is a 
    twisted track.  This procedure has been written to handle twisted tracks 
    in full but the following clause re-submits 2t in the case t is twisted 
    and returns control to the call.  Should the decomposition given by 
    twisted tracks be required, the switch 't' should be included at the 
    command line.
    ************************************************************************/

    if (twisted && !TWISTED_ALLOWED)
    {
		if (echo)
		{
		    cout << "\nThe track\n";
			for (unsigned int i=0; i< track.size(); i++)
				cout << track[i] << ' ';
		    cout << "\nis twisted. Decomposing 2t." << endl;
		}
		       
		vector<int> dbtrack(track); /* "doubled track" */
		for (unsigned int i = 0; i<dbtrack.size(); i++)
		    dbtrack[i] *= 2;

if (decomp_control::DECOMPOSE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
    debug << "decompose: track is twisted, decomposing 2t" << endl;

		return decompose (dbtrack, K, echo, full_decomp, txt);
    }

    if ( !twisted )
    {
		/*****************************************************************
		The tree is not necessarily complete, we can only say, at
		this stage, that t is untwisted.  We now see if we can extend the
		tree from one of the 0-cells just added to a 1-cell b-label not
		already in the tree.  That this is sufficient to continue
		the construction follows from remark (*) above. 
		******************************************************************/
		bool extension_found = false;

		for ( unsigned int i = max1cbl; i < jbba.numrows(); i++ )
		{ 
	    	if ( jbba[i][1] < 0 )
	    	{
				/* The ith b-label (a 0-cell b-label) is in the tree. */
				if ( jbba[i][0] <= max1cbl && jbba[jbba[i][0]-1][1] > 0
				   )
				   /*******************************************************
		    		Again the second half of the above condition is changed
		    		from the 1987 version, since the 1-cell b-label 
		    		*ptr_jbba->elt(i,0) is not in the tree iff the label at 
		    		the other end of the corresponding disc arc is not in the
		    		tree.
		    		*******************************************************/
				{
		    		/* i is joined to the 1-cell b-label 'jbba[i,1]', which
		    		   is not in the tree. */
		    		new_label = jbba[i][0];
					btrarcs1[amarker][0] = i+1;
					btrarcs1[amarker][1] = new_label;					
		    		extension_arc = ++amarker;
		    		extension_found = true;
		    		break;
				}
				else
				if ( jbba[i][2] <= max1cbl && jbba[jbba[i][2]-1][1] > 0
				   )
				{
		    		/* i is joined to the 1-cell b-label 'jbba[i,1]', which
		    		   is not in the tree. */
		    		new_label = jbba[i][2];
					btrarcs1[amarker][0] = i+1;
					btrarcs1[amarker][1] = new_label;									
		    		extension_arc = ++amarker;            
		    		extension_found = true;
		    		break;
				}
	    	}
		}
	
		if ( extension_found )
		{
	    	/*************************************************************
	    	t is untwisted but non separating.  We have adjusted
	    	'extension_arc' accordingly and now extend the tree over the
	    	(other) component of the boundary of B from the newly assigned
	    	'new label'.  'twisted' will be set true within the
	    	next call to 'band extend', so 'twisted' may be used to indicate
	    	whether or not t separates.

	    	Reinitialize for next call to 'band extend'.  We do not adjust
	    	'am' since X is connected, 'bposition' will have been at most
	    	half assigned to, so we continue to assign to the remainder,
	    	noting that 'ppointer' is correctly assigned for this
	    	continuation.
	    	**************************************************************/

	    	lpointer = 0;
	    	pmarker = ppointer;
			
			band_extend(new_label, number_of_1_cells_in_rosette, occurrences_of_1_cell, blabels_for_1_cell, 
                lpointer, bposition, ppointer, pmarker, jbba, btrarcs1, amarker, twisted);
						
	    	/* Adjust 'max band arc 2' to its 
	    	   correct value for this case. */
	    	max_band_arc_2 = amarker;

	    	/*Add any remaining o-cells and assign the number of edges in the
	    	  b-tree to 'nebt1',overwriting the previous assignation. */
			add_0_cells(number_of_0_cells, blabels_for_0_cell, jbba, amarker, btrarcs1);

	    	nebt1 = amarker;
		}
		else
		{
	    	/***********************************************************
	    	t separates, we have completed the maximal tree in the first
	    	component of X and now use the parallel properties mentioned
	    	above to write down a maximal tree in the other component of the
	    	boundary of B.
	    	***********************************************************/
	    	for ( int i =0; i< max_band_arc_1 ; i++)
	    	/* max_band_arc_1 is an upper bound for the part
	    	   of 'btrarcs1' assigned by 'band extend' */
	    	{
				if ( btrarcs1[i][0] % 2 == 1 )
				{
					btrarcs2[i][0] = btrarcs1[i][0]+1;
					btrarcs2[i][1] = btrarcs1[i][1]-1;
				}
				else
				{
					btrarcs2[i][0] = btrarcs1[i][0]-1;
					btrarcs2[i][1] = btrarcs1[i][1]+1;
				}
	    	}

	    	/* All the remaining 1-cell b-labels have been added to the second
	    	  tree, we adjust 'jbba[,2]' to signify this. */
	    	for ( int i =0; i< max1cbl ; i++)
	    		jbba[i][1] = -abs( jbba[i][1]);

	    	/* 'am' now becomes a marker for 'btrarcs2' as we add 
			any remaining 0-cells to this second tree. */
	    	amarker = max_band_arc_1;
			add_0_cells(number_of_0_cells, blabels_for_0_cell, jbba, amarker, btrarcs2);

	    	/* We, finally, assign the number of edges in the second 
	    	   b-tree to 'nebt2'. */
	    	nebt2 = amarker;
		}
    }

if (decomp_control::DECOMPOSE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
{

	debug << "decompose: btrarcs1" << endl;
	print(btrarcs1,debug,0,"decompose:   ");
	debug << "decompose: btrarcs2" << endl;
	print(btrarcs2,debug,0,"decompose:   ");
	debug << "decompose: bposition" << endl;
	print(bposition,debug,0,"decompose:   ");

	debug << "decompose: nett " << nett << " max_band_arc_1 " << max_band_arc_1
	      << " max_0_cell_arc " << max_0_cell_arc << " extension_arc "
	      << extension_arc << endl;
	debug << "decompose: max_band_arc_2 " << max_band_arc_2 << " nebt1 "
	      << nebt1 << " nebt2 " << nebt2 << endl;	  
}


	if (!full_decomp)
	{
		txt << "The track (";
		for(unsigned int i=0; i< track.size()-1; i++)
			txt << track[i] << ", ";       
				
		txt << track[track.size()-1] << ") ";
			
		if (echo)
		{
			cout << "The track (";
			for(unsigned int i=0; i< track.size()-1; i++)
				cout << track[i] << ", ";       
					
			cout << track[track.size()-1] << ") ";
		}

		/* the boolean twisted actually indicates whether the
		   track is separating at this stage, extension_arc
		   tells us whether it is twisted or not
		*/
		if ( twisted && extension_arc == 0)
		{
			txt << " is twisted" << endl;
			
			if (echo)
				cout << " is twisted" << endl;
			
			status |= DECOMPOSE_TWISTED;
			status |= DECOMPOSE_NON_SEPARATING;
		}
		else if ( twisted )
		{
			txt << " is untwisted and non-separating" << endl;

			if (echo)
				cout << " is untwisted and non-separating" << endl;
				
			status |= DECOMPOSE_NON_SEPARATING;
		}
		else
		{
			txt << " is untwisted and separating" << endl;
			
			if (echo)
				cout << " is untwisted and separating" << endl;
		}
		
		return status; // we don't want to calculate induced images
	}

    /* The maximal trees are now complete and we may begin calculating 
    induced images.  First, however, we create the flags on the 1-cells 
    in the rosette to enable this calculation. */

	vector<string> flag_on_1_cell(number_of_1_cells_in_rosette);
	
    for (int i =0; i< number_of_1_cells_in_rosette; i++)
	flag_on_1_cell[i] = string(1, generator_corresponding_to_1_cell[i]);
   
    /* We must remove any elements of 'flag on 1 cell' corresponding to the 
       first 1-cell "carried" by a 2-torsion generator in R(L).
       Since 2-torsion generators are entered first, this is not difficult. 
	*/
    for (int i=0; i< number_of_0_cells-1; i++)
    {
		if (FLAGS_ON_FIRST_2_TORSION_1_CELL)
		{
			flag_on_1_cell[2*i+1] = "";
		}
		else
			flag_on_1_cell[2*i] = "";
	}

if (decomp_control::DECOMPOSE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "decompose: flag_on_1_cell: ";

	for (int i=0; i< number_of_1_cells_in_rosette; i++)
	{
		if (flag_on_1_cell[i] == "" )
			debug << ". ";
		else
			debug << flag_on_1_cell[i] << ' ';
	}
	debug << endl;
}

    /**********************************************************************
    We work towards a procedure 'image', for evaluating paths (and their
    corresponding words) between arbitrary t- or b-labels in some tree and a
    chosen basepoint in that tree.  To this end, we create two []INT objects,
    "ttree" and "btree".  We shall regard the trees as being oriented from 
    our chosen basepoint.  In 'ttree[i]' ('btree[i]') we shall store the 
    unique vertex of the tree that is joined by an edge, e, to the ith 
    vertex, such that the terminal point of e is the ith vertex 
	(2006: i.e. the parent vertex of i (i+1 since we label vertices from 1) in the tree).  

	We also record, for the ith vertex, the word corresponding to moving along e 
    (e as above) from i (towards the base point).
    ************************************************************************/
    int bsize = nebt1 + nebt2 + (twisted ? 1 : 2);
    /* 'bsize' is the total number of b-tree vertices. */

     /* Declare space for the above records. */
    int ttree[nett+1];
	for (int i=0; i< nett+1; i++)
		ttree[i] = 0;
	
    int btree[bsize];
	for (int i=0; i< bsize; i++)
		btree[i] = 0;

    string tvw[nett+1];
    string bvw[bsize];

    /****************************************************************
    We consider the vertices of the various trees labelled as follows.
    Suppose l is a label attached to some vertex, v.  If l occurrs in the ith
    column of the corresponding label matrices (taken in order) then v is
    labelled by i.  If l is attached to the jth 0-cell and there are n b-label
    columns in all, v is labelled by n+j.  We now create copies of the tree
    arc matrices with the new labelling.  
    ******************************************************************/
    matrix<int> rttrarcs(nett,2); /* "relabelled t-tree arcs" */
    matrix<int> rbtrarcs1(nebt1,2);
    matrix<int> rbtrarcs2(nebt2,2);

    int tccc[number_of_1_cells_in_rosette];
    int bccc[number_of_1_cells_in_rosette];
    /* "t(b)-label cumulative column count" */

    /* Assign 'tccc' and 'bccc'. */
    tccc[0] = 0;
    bccc[0] = 0;

    for ( int i=1; i < number_of_1_cells_in_rosette; i++)
    {
		tccc[i] = tccc[i-1] + n_tpts_on_1_cell[i-1];
		bccc[i] = bccc[i-1] + 2 * n_tpts_on_1_cell[i-1];
    }

    /* assign the number of t- and b-label columns for future use */
    int numbcols = bccc[number_of_1_cells_in_rosette-1] + 2 * n_tpts_on_1_cell[number_of_1_cells_in_rosette-1];

    /* Assign relabelled arcs. */
    vector<int> position(3);

    for (int i = 0; i< nett; i++)
	for (int j = 0; j < 2; j++)
	{

	    find(ttrarcs[i][j],number_of_1_cells_in_rosette,tlabels_for_1_cell, position);
	    rttrarcs[i][j] = tccc[position[0]] + position[2] + 1;
	}

    for (int i = 0; i< nebt1; i++)
	for (int j = 0; j < 2; j++)
	{
	    if (btrarcs1[i][j] <= max1cbl )
	    {
			find(btrarcs1[i][j], number_of_1_cells_in_rosette, blabels_for_1_cell, position);
			rbtrarcs1[i][j] =  bccc[position[0]] + position[2] + 1;
	    }
	    else
			rbtrarcs1[i][j] = bvertex(btrarcs1[i][j], number_of_0_cells, numbcols, blabels_for_0_cell);
	}

    for (int i = 0; i< nebt2; i++)
	for (int j = 0; j < 2; j++)
	{
	    if (btrarcs2[i][j] <= max1cbl )
	    {
			find(btrarcs2[i][j], number_of_1_cells_in_rosette, blabels_for_1_cell, position);
			rbtrarcs2[i][j] =  bccc[position[0]] + position[2] + 1;
	    }
	    else
			rbtrarcs2[i][j] = bvertex(btrarcs2[i][j], number_of_0_cells, numbcols, blabels_for_0_cell);
	}

if (decomp_control::DECOMPOSE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "decompose: bsize = " << bsize << endl; //" max_disc_word_size = "   << max_disc_word_size;

	debug << "decompose: rttrarcs" << endl;
	print (rttrarcs,debug,0,"decompose:   ");
	debug << "decompose: rbtrarcs1" << endl;
	print(rbtrarcs1,debug,0,"decompose:   ");
	debug << "decompose: rbtrarcs2" << endl;
	print (rbtrarcs2,debug,0,"decompose:   ");

	debug << "decompose: cumulative t-label column count";
	for (int i =0; i < number_of_1_cells_in_rosette; i++)	
		debug << setw(4) << tccc[i] << ' ';
	debug << endl;

	debug << "decompose: cumulative b-label column count";
	for (int i =0; i < number_of_1_cells_in_rosette; i++)	
		debug << setw(4) << bccc[i] << ' ';
	debug << endl;
}
    /*********************************************************************
    We now work with this new labelling.  We shall choose the t-label 1 as 
    basepoint for the t-graph, b-label 1 as basepoint for the
    first component of the b-graph and, if the b-graph is not connected, 
    b-label 2 as basepoint for the second component of the b-graph.  This 
    ensures we calculate the correct conjugates of the induced images in G.  
    Suppose v0 is one such basepoint.  By construction, the natural
    orientation of the tree arcs (column1-->column2) is "away" from v0, for
    those arcs in t or the boundary of B (the construction was
    "breadth-first") or the extension arc (if it exists).  For arcs adding
    0-cells of L to the b-tree(s), the orientation is "towards" v0.  Notice
    we cannot be sure where the 0-cells of L are added to the b-tree(s)
    (ie, at which level).
    *******************************************************************/

    /* Assign 'ttree' and 'tvw'. */
    for (int i= 0; i< nett; i++)
    {
		ttree[rttrarcs[i][1]-1] = rttrarcs[i][0];
		tvw[rttrarcs[i][1]-1] = disc_word(ttrarcs[i][1], ttrarcs[i][0], true, true, mbtl, tlabels_for_1_cell, 
										  flag_on_1_cell, n_tpts_on_1_cell, ext0cbl, max1cbl,
										  number_of_1_cells_in_rosette, length_of_attaching_map, attaching_map);
    }

    /* Assign 'btree' and 'bvw'.*/
    for (int i=0; i < nebt1; i++ )
    {
		if ( i< max_band_arc_1 || (i >=extension_arc-1 && i < max_band_arc_2))
		{
		    btree[rbtrarcs1[i][1]-1] = rbtrarcs1[i][0];
		    bvw[rbtrarcs1[i][1]-1] = disc_word(btrarcs1[i][1], btrarcs1[i][0],false, true, m1cbl, blabels_for_1_cell, 
										  flag_on_1_cell, n_tpts_on_1_cell, ext0cbl, max1cbl,
										  number_of_1_cells_in_rosette, length_of_attaching_map, attaching_map);
		}
		else
		{
		    btree[rbtrarcs1[i][0]-1] = rbtrarcs1[i][1];
		    bvw[rbtrarcs1[i][0]-1] = disc_word(btrarcs1[i][0],btrarcs1[i][1],false, true, m1cbl, blabels_for_1_cell, 
										  flag_on_1_cell, n_tpts_on_1_cell, ext0cbl, max1cbl,
										  number_of_1_cells_in_rosette, length_of_attaching_map, attaching_map);
		}
    }

    if (!twisted)
    {
		for(int i=0; i < max_band_arc_1; i++)
		{
		    btree[rbtrarcs2[i][1]-1] = rbtrarcs2[i][0];
		    bvw[rbtrarcs2[i][1]-1] = disc_word(btrarcs2[i][1],btrarcs2[i][0],false,true, m1cbl, blabels_for_1_cell, 
										  flag_on_1_cell, n_tpts_on_1_cell, ext0cbl, max1cbl,
										  number_of_1_cells_in_rosette, length_of_attaching_map, attaching_map);
		}

		for (int i = max_band_arc_1; i< nebt2; i++)
		{
		    btree[rbtrarcs2[i][0]-1] = rbtrarcs2[i][1];
		    bvw[rbtrarcs2[i][0]-1] = disc_word(btrarcs2[i][0],btrarcs2[i][1],false,true, m1cbl, blabels_for_1_cell, 
										  flag_on_1_cell, n_tpts_on_1_cell, ext0cbl, max1cbl,
										  number_of_1_cells_in_rosette, length_of_attaching_map, attaching_map);
		}
    }

    /************************************************************* 
     *  We now note the base point vertices, "t-graph base point", 
     *  "b-graph base point1(2)". 
     ***************************************************************/
	int tbp;
	int bbp1;
	int bbp2;
	
    if (trosette)
    {
	    find(1,number_of_1_cells_in_rosette,tlabels_for_1_cell, position);
		tbp = tccc[position[0]] + position[2] + 1;
    }
    else
		tbp = rttrarcs[0][0];

    find(1,number_of_1_cells_in_rosette,blabels_for_1_cell, position);
    bbp1 = bccc[position[0]] + position[2] + 1;
    
    if (twisted)
		bbp2 = 0;
    else
    {
	    find(2,number_of_1_cells_in_rosette,blabels_for_1_cell, position);
		bbp2 = bccc[position[0]] + position[2] + 1;
    }

if (decomp_control::DECOMPOSE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "decompose: relabelled t-tree: "; 
	for (int i=0; i< nett+1; i++)
		debug << ttree[i] << ' ';
	debug << endl;

	debug << "decompose: t-tree vertex words" << endl;
	for (int i=0; i<nett; i++)
	    debug << "decompose:   vertex " << rttrarcs[i][1] << ' ' << tvw[rttrarcs[i][1]-1] << endl;

	debug << "decompose: relabelled b-tree ";
	for (int i=0; i< bsize; i++)
		debug << btree[i] << ' ';
	debug << endl;

	debug << "decompose: b-tree vertex words" << endl;
	for (int i=0; i<nebt1; i++)
	{
	    if (i<max_band_arc_1 || (i >=extension_arc-1 && i < max_band_arc_2))
			debug << "decompose:   vertex " << rbtrarcs1[i][1] << ' ' << bvw[rbtrarcs1[i][1]-1] << endl;
	    else
			debug << "decompose:   vertex " << rbtrarcs1[i][0] << ' ' << bvw[rbtrarcs1[i][0]-1] << endl;
	}
    
	if (!twisted)
	{
	    for(int i=0; i < max_band_arc_1; i++)
			debug << "decompose: vertex " << rbtrarcs2[i][1]  << ' ' << bvw[rbtrarcs2[i][1]-1] << endl;

	    for (int i = max_band_arc_1; i< nebt2; i++)
			debug << "decompose: vertex " << rbtrarcs2[i][0]  << ' ' << bvw[rbtrarcs2[i][0]-1] << endl;
	}

	debug << "decompose: tbp = " << tbp << " bbp1 = " << bbp1 << " bbp2 = " << bbp2 << endl;
}

    /* now we deal with the t-graph to calculate the edge stabilizer */
    int nesg = tndtc-nett; /* "number of edge stabilizer generators" */
	
	string estgen[nesg];
    int gentr = 0; /* "generator" will act as a marker */

    /* We assign 'tael[i,]' to be negative iff the corresponding arc 
       lies in the t-tree. The elements of tael are currently all
       negative, so change the sign first.
    */
    for (int i=0; i< tndtc; i++)
	for (int j=0; j<2; j++)
		tael[i][j] *= -1;

    for ( int i=0; i< nett; i++)
    {
		/* look for ttrarcs[i,0] in tael */
		int tlabel = ttrarcs[i][0];
		for (int j=0; j< tndtc; j++)
		for (int k=0; k<2; k++)
		{
	    	if ( abs(tael[j][k]) == tlabel )
	    	{
				/* set the arc elements negative */
				tael[j][0] *= -1;
				tael[j][1] *= -1;
				goto t_arc_set_negative;
	    	}
		}
		t_arc_set_negative:;
	}

   /* Evaluate the edge stabilizer. */

    for (int i=0; i< tndtc; i++)
	{
		if (tael[i][0] > 0)
		{		
	    	estgen[gentr++] = image(tael[i][0], tael[i][1], true, ttree, tvw, mbtl, tlabels_for_1_cell,
									twisted, tbp, bbp1, bbp2, tccc, bccc, number_of_0_cells, numbcols, 
									blabels_for_0_cell, flag_on_1_cell, n_tpts_on_1_cell, ext0cbl, max1cbl, 
									number_of_1_cells_in_rosette, length_of_attaching_map, attaching_map );
		}
	}
	
if (decomp_control::DECOMPOSE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "decompose: edge stabilizer generator" << endl;
	for (int i=0; i<nesg; i++)
	    debug << "decompose:   generator " << i+1   << " " << estgen[i] << endl;
}

    /* For the b-graph, we store all the generators for the vertex stabilizer(s) 
	   in one array of strings. We first evaluate the required upper bound for such an object. */

    int nvsg = 2*tndtc;

    for (int i = 0 ; i< number_of_1_cells_in_rosette; i++)
		nvsg += n_tpts_on_1_cell[i]+1;
    /* 'nvsg' now contains #(EX). */

    nvsg -= (nebt1+nebt2);
	string vstgen[nvsg];
    gentr=0; /* Return 'gentr' for use as a marker. */

    /* Return 'jbba[,2]' to its absolute value. */
    for (unsigned int i =0 ; i< jbba.numrows() ; i++)
		jbba[i][1] *= -1;

    /* Remove any repetitions amongst the entries in 'jbba': the entry jbba[i][j] represents
	   an arc from (i+1)-->jbba[i][j] and the reverse of this arc jbba[i][j]-->(i+1) appears
	   in row jbba[i][j]-1.  We set the second of these arc representations to be negative.
    */

    for (int i =0; i< max1cbl; i++)
	for (int j=0; j < 3 ; j++)
	{
	    if (jbba[i][j]>0)
			jbba[jbba[i][j]-1][2-j] *= -1;
	}
					     
    for (unsigned int i=max1cbl; i< jbba.numrows(); i++)
    {
		if(jbba[i][0] > 0)
		    jbba[jbba[i][0]-1][2] *= -1;
		if(jbba[i][2] > 0)
		    jbba[jbba[i][2]-1][0] *= -1;
    }

    /* Set 'jbba[i][1]' to be negative if (i+1)-->jbba[i][1] is a "band arc"
       in a b-tree. */

    for (int i=0; i<nebt1; i++)
	{
		if (i < max_band_arc_1 || (i >=extension_arc && i < max_band_arc_2))
		{
		    jbba[btrarcs1[i][0]-1][1]= -abs(jbba[btrarcs1[i][0]-1][1]);
		    jbba[btrarcs1[i][1]-1][1]= -abs(jbba[btrarcs1[i][1]-1][1]);
		}
	}
	
    if (!twisted)
	{
		for (int i =0; i< max_band_arc_1; i++)
		{
	    	jbba[btrarcs2[i][0]-1][1]= -abs(jbba[btrarcs2[i][0]-1][1]);
		    jbba[btrarcs2[i][1]-1][1]= -abs(jbba[btrarcs2[i][1]-1][1]);
		}
	}
	
    /********************************************************************
    Set negative 'jbba[i][j]' if (i+1)-->jbba[i][j] is a "boundary arc" in a 
    b-tree.  Suppose the labels in 'btrarcs1(2)[i][]' are l1,l2 respectively 
    and indicate a boundary arc (this includes the extension arc, if it 
    exists); l1 will always be a 0-cell b-label but l2 may be a 0-cell or 
    1-cell b-label.  If l2 is a 1-cell b-label, we isolate the row of 
    b-labels containing l1 and the column of b-labels containing l2 then 
    check for all occurrences of an arc identified with l1-->l2 in R(L).  
    By the way we removed repetitions, all occurrences of the arcs will 
    occur in the rows of 'jbba' determined by the column of b-labels 
    containing l2.   Here we have to look "both sides" of l2 for possible 
    occurrences.  If l2 is a 0-cell b-label the arc l1-->l2 must correspond 
    to a 1-cell in the rosette, "carried" by a generator of order 2, and l2 
    must be attached to the first 0-cell in L. In this case we utilize 
    'b0cbl for 1 cell' to set the required occurrences of the arc l1-->l2 
    to be negative.  Notice the latter case can only occur in the first 
    component of X, since the first 0-cell b-label lies in the same 
    component of L-t as the b-label 1.
    *******************************************************************/
if (decomp_control::DECOMPOSE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << "decompose: setting negative boundary arcs in jbba" << endl;

    for (int i = max_band_arc_1; i< nebt1; i++)
    {
		if ( extension_arc == 0 || (extension_arc != 0 && ( i < extension_arc || i >= max_band_arc_2 ) ) )
		{
	    	if ( btrarcs1[i][1] <= max1cbl )
	    	{
				/* second label in arc is a 1-cell b-label */
			    find(btrarcs1[i][1],number_of_1_cells_in_rosette,blabels_for_1_cell, position);
				int one_cell = position[0]; /* carrying l2 */
				int column = position[2]; /* carrying l2 */
				matrix<int>& matrixptr = *blabels_for_1_cell[one_cell];

				vector<int>& Ocblptr = *find_row (btrarcs1[i][0], number_of_0_cells, blabels_for_0_cell);

				for (int j = 0; j < occurrences_of_1_cell[one_cell]; j++ ) /* down column containing l2 */
				{
		    		for (unsigned int k=0; k < Ocblptr.size(); k++) /* along row containing l1 */
		    		{                         
						if (jbba[matrixptr[j][column]-1][0] == Ocblptr[k])
						{
			    			jbba[matrixptr[j][column]-1][0] *= -1;
			    			break;
						}
						else 
						if (jbba[matrixptr[j][column]-1][2] == Ocblptr[k])
						{
			    			jbba[matrixptr[j][column]-1][2] *= -1;
			    			break;
						}
		    		}
				}
	    	}
	    	else
	    	{
if (decomp_control::DECOMPOSE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
	debug << endl << "decompose: label lies on added 0-cell" << endl;
				
				/* l1 lies on an added o-cell and so the 1-cell carrying
				   l1-->l2 may be determined from l1, since 2-torsion
				   generators were entered first and their corresponding
				   relators were entered first, in the positive.
				*/
				int one_cell = 2*((btrarcs1[i][0]-max1cbl)/4+1)-1;

				/* We now use the way we removed repetitions and the
				  construction of 'b0cbl_for_1_cell' to set negative the
				  required occurrences. */
				matrix<int>& matrixptr = *b0cbl_for_1_cell[one_cell-1];

				for (unsigned int j=0; j < matrixptr.numrows(); j++)
				{
		    		if (jbba[matrixptr[j][0]-1][0] == matrixptr[j][1])
						jbba[matrixptr[j][0]-1][0] *= -1;
		    		else
						jbba[matrixptr[j][0]-1][2] *= -1;
				}
	    	}
		}
    }

    if (!twisted)
    {
		for (int i = max_band_arc_1 ; i< nebt2; i++)
		{
	    	/* the second label in a tree arc must be a 1-cell b-label 
	    	   in the second component */
	    	find (btrarcs2[i][1],number_of_1_cells_in_rosette,blabels_for_1_cell,position); /* l2 */
	    	int one_cell = position[0]; /* carrying l2 */
	    	int column = position[2]; /* carrying l2 */
	    	matrix<int>& matrixptr = *blabels_for_1_cell[one_cell];

			vector<int>& Ocblptr = *find_row (btrarcs2[i][0], number_of_0_cells, blabels_for_0_cell);

	    	for (int j = 0; j < occurrences_of_1_cell[one_cell]; j++) /* down column containing l2 */
	    	{
				for (unsigned int k=0; k < Ocblptr.size(); k++) 
				{                         /* along row containing l1 */
		    		if (jbba[matrixptr[j][column]-1][0] == Ocblptr[k])
		    		{
						jbba[matrixptr[j][column]-1][0] *= -1;
						break;
		    		}
		    		else 
		    		if (jbba[matrixptr[j][column]-1][2] == Ocblptr[k])
		    		{
						jbba[matrixptr[j][column]-1][2] *= -1;
						break;
		    		}
				}
	    	}
		}
    }

    /**********************************************************************
    We are now left with the positive entries of 'jbba' denoting edges in X
    representing generators for one (or other) vertex stabilizer; we still
    have all the occurrences of the "boundary" generators but no inverses.  
    If the arc l1-->l2 is present, it occurrs in the ith row of 'jbba', where
    (i+1)=min(l1,l2).

    We use 'bposition' and 'btrarcs1' to deal with the first component.  All
    the b-labels attached to elements of VX in this component are stored in
    the columns of the b-label matrices referenced by 'bposition', or the
    rows of 0-cell b-labels added to the first "band forest" (a tree iff t
    separates).
    **********************************************************************/

    for (unsigned int i=0; i < bposition.numrows(); i++)
    {
		if ( bposition[i][0] == -1 )
	    	break;

		/* We check whether any of the arcs from the labels in the
		   column indicated by the ith row of bposition give generators 
		   for the first component.
		*/

		int one_cell = bposition[i][0];
		int column = bposition[i][1];
		matrix<int>& matrixptr = *blabels_for_1_cell[one_cell];
		for (int j = 0; j < occurrences_of_1_cell[one_cell]; j++) /* down column containing l2 */
		{
	    	for ( int k=0; k < 3 ; k++)
	    	{
				if ( jbba[matrixptr[j][column]-1][k] > 0 )
				{
			    	vstgen[gentr++] = image(matrixptr[j][column], jbba[matrixptr[j][column]-1][k], false, btree, 
											bvw, m1cbl, blabels_for_1_cell,
											twisted, tbp, bbp1, bbp2, tccc, bccc, number_of_0_cells, numbcols, 
											blabels_for_0_cell, flag_on_1_cell, n_tpts_on_1_cell, ext0cbl, max1cbl, 
											number_of_1_cells_in_rosette, length_of_attaching_map, attaching_map );

		    		if ( k == 1 )
						jbba[matrixptr[j][column]-1][k] *= -1;
		    		else /*1*/
		    		{
						set_negative(matrixptr[j][column],jbba[matrixptr[j][column]-1][k], jbba, max1cbl,
									 number_of_1_cells_in_rosette, blabels_for_1_cell, b0cbl_for_1_cell,
									 occurrences_of_1_cell, length_of_attaching_map, attaching_map, ext0cbl);
		    		}
				}
	    	}
		}
    }

 /**************************************************************************
 We set the first component arcs negative above to prevent us taking this 
 generator again, and to assist the calculation of the second vertex 
 stabilizer (if it exists).  Notice, at *1*, 'jbba[col[j],2]' was the only 
 occurrence of this arc left positive in 'jbba' after the repetitions were 
 removed.

 Notice that we may ignore the labels of the extension arc (if it exists)
 since the first will be an added 0-cell and the second will have its
 position in the b-label matrices referenced by a row of 'bposition'.
 *************************************************************************/

    for (int i = max_band_arc_1; i< nebt1; i++)
    {
		if ( extension_arc == 0 || (extension_arc != 0 && (i < extension_arc || i >= max_band_arc_2)))
		{
			vector<int>& Ocblptr = *find_row (btrarcs1[i][0], number_of_0_cells, blabels_for_0_cell);
	    	for (unsigned int j =0; j<Ocblptr.size(); j++ )
	    	{
				if (jbba[Ocblptr[j]-1][0] > 0)
				{
			    	vstgen[gentr++] = image(Ocblptr[j],jbba[Ocblptr[j]-1][0], false, btree, 
											bvw, m1cbl, blabels_for_1_cell,
											twisted, tbp, bbp1, bbp2, tccc, bccc, number_of_0_cells, numbcols, 
											blabels_for_0_cell, flag_on_1_cell, n_tpts_on_1_cell, ext0cbl, max1cbl, 
											number_of_1_cells_in_rosette, length_of_attaching_map, attaching_map );
					   
						set_negative(Ocblptr[j], jbba[Ocblptr[j]-1][0], jbba, max1cbl,
									 number_of_1_cells_in_rosette, blabels_for_1_cell, b0cbl_for_1_cell,
									 occurrences_of_1_cell, length_of_attaching_map, attaching_map, ext0cbl);
				}
				if (jbba[Ocblptr[j]-1][2] > 0 )
				{
			    	vstgen[gentr++] = image(Ocblptr[j],jbba[Ocblptr[j]-1][2], false, btree, 
											bvw, m1cbl, blabels_for_1_cell,
											twisted, tbp, bbp1, bbp2, tccc, bccc, number_of_0_cells, numbcols, 
											blabels_for_0_cell, flag_on_1_cell, n_tpts_on_1_cell, ext0cbl, max1cbl, 
											number_of_1_cells_in_rosette, length_of_attaching_map, attaching_map );
					   
						set_negative(Ocblptr[j], jbba[Ocblptr[j]-1][2], jbba, max1cbl,
									 number_of_1_cells_in_rosette, blabels_for_1_cell, b0cbl_for_1_cell,
									 occurrences_of_1_cell, length_of_attaching_map, attaching_map, ext0cbl);
				}
	    	}
		}
    }

if (decomp_control::DECOMPOSE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "decompose: first vertex stabilizer generator" << endl;
	for (int i=0; i<gentr; i++)
	    debug << "decompose:   generator " << i+1  << " " << vstgen[i] << endl;
}

    /*********************************************************************
    This completes the first vertex stabilizer.  We have left, positive, in
    'jbba' only those arcs generating the second vertex stabilizer (if it
    exists), so we may simply run through 'jbba' to find them.  We note at
    this stage how many generators we already have in 'vstgen'.
    **********************************************************************/

    int n_gen_for_1st_vstab = gentr; // "number of generators for first stabilizer"

    for (int i = 0; i< max1cbl; i++)
	for (int j = 0; j < 3; j++)
	{
	    if ( jbba[i][j] > 0 )
	    {
	    	vstgen[gentr++] = image(i+1, jbba[i][j], false, btree, bvw, m1cbl, blabels_for_1_cell,
									twisted, tbp, bbp1, bbp2, tccc, bccc, number_of_0_cells, numbcols, 
									blabels_for_0_cell, flag_on_1_cell, n_tpts_on_1_cell, ext0cbl, max1cbl, 
									number_of_1_cells_in_rosette, length_of_attaching_map, attaching_map );

			if ( j == 1 )
		    	jbba[i][1] *= -1;
			else 
			{
				set_negative(i+1,jbba[i][j], jbba, max1cbl,
							 number_of_1_cells_in_rosette, blabels_for_1_cell, b0cbl_for_1_cell,
							 occurrences_of_1_cell, length_of_attaching_map, attaching_map, ext0cbl);
			}
	    }
	}
	
    for (unsigned int i = max1cbl; i< jbba.numrows(); i++)
    {
		if ( jbba[i][0] > 0 )
		{
	    	vstgen[gentr++] = image(i+1, jbba[i][0], false, btree, bvw, m1cbl, blabels_for_1_cell,
									twisted, tbp, bbp1, bbp2, tccc, bccc, number_of_0_cells, numbcols, 
									blabels_for_0_cell, flag_on_1_cell, n_tpts_on_1_cell, ext0cbl, max1cbl, 
									number_of_1_cells_in_rosette, length_of_attaching_map, attaching_map );

				set_negative(i+1,jbba[i][0], jbba, max1cbl,
							 number_of_1_cells_in_rosette, blabels_for_1_cell, b0cbl_for_1_cell,
							 occurrences_of_1_cell, length_of_attaching_map, attaching_map, ext0cbl);
		}
		if ( jbba[i][2] > 0 )
		{
	    	vstgen[gentr++] = image(i+1, jbba[i][2], false, btree, bvw, m1cbl, blabels_for_1_cell,
									twisted, tbp, bbp1, bbp2, tccc, bccc, number_of_0_cells, numbcols, 
									blabels_for_0_cell, flag_on_1_cell, n_tpts_on_1_cell, ext0cbl, max1cbl, 
									number_of_1_cells_in_rosette, length_of_attaching_map, attaching_map );

				set_negative(i+1,jbba[i][2], jbba, max1cbl,
							 number_of_1_cells_in_rosette, blabels_for_1_cell, b0cbl_for_1_cell,
							 occurrences_of_1_cell, length_of_attaching_map, attaching_map, ext0cbl);
		}
    }

if (decomp_control::DECOMPOSE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "decompose: second vertex stabilizer generator" << endl;
	for (int i=n_gen_for_1st_vstab; i<gentr; i++)
	    debug << "decompose:   generator " << i-n_gen_for_1st_vstab+1 << " " << vstgen[i] << endl;
}

    /***********************************************************************
    This completes the evaluation of the stabilizers.  We now reassign 'nvsg'
    to be the actual number of vertex stabilizer generators, rather than our
    initial upper bound.
    **********************************************************************/
    nvsg = gentr;

    /* We now reduce the generators of the stabilizers 
       and write out the resulting generators into the file 'txt'. */
    reduce(estgen, 0, nesg, relator, number_of_relators);
    reduce(vstgen, 0, n_gen_for_1st_vstab, relator, number_of_relators);
    reduce(vstgen,n_gen_for_1st_vstab, nvsg - n_gen_for_1st_vstab, relator, number_of_relators);

if (decomp_control::DECOMPOSE_DEBUG && decomp_control::DEBUG >= decomp_control::DETAIL)
{
	debug << "decompose: reduced edge stabilizer generators" << endl;
	for (int i=0; i< nesg; i++)
	    debug << "decompose:   generator " << i+1  << " <" << estgen[i] << ">" << endl;

	debug << "decompose: reduced first vertex stabilizer generator" << endl;
	for (int i=0; i<n_gen_for_1st_vstab; i++)
	    debug << "decompose:   generator " << i+1  << " <" << vstgen[i] << ">"<< endl;

	debug << "decompose: reduced second vertex stabilizer generator" << endl;
	for (int i=n_gen_for_1st_vstab; i<gentr; i++)
	    debug << "decompose:   generator " << i-n_gen_for_1st_vstab+1 << " <" << vstgen[i] << ">"<< endl;
}

    if ( echo )
    {
		if ( twisted && extension_arc == 0)
		    cout << "\nThe twisted track\n\n";
		else if ( twisted )
		    cout << "\nThe non-separating and untwisted track\n\n";
		else
		    cout << "\nThe separating track\n\n";

		for(unsigned int i=0; i< track.size()-1; i++)
			cout << track[i] << ", ";       
		cout << track[track.size()-1] << endl;
    }

    /* check for an obvious trivial decomposition */
    bool trivial_decomp = false;
    int single_gen = 0;
    int non_trivial_generators=0; // used to count number of non-trivial generators
	
    for (int i=0; i<nesg; i++)
	{		
		non_trivial_generators += estgen[i].length();
		
		if (estgen[i].length() == 1)
		    single_gen++;
	}
    
	if (single_gen == number_of_generators)
	{
		trivial_decomp = true;
		estgen[0] = "G";
		for (int i=1; i< nesg; i++)
			estgen[i] = "";
	}

	if (non_trivial_generators == 0)
	    estgen[0] = "1";
	    
	single_gen = 0;
	non_trivial_generators=0;

	for (int i=0; i<n_gen_for_1st_vstab; i++)
	{
		non_trivial_generators += vstgen[i].length();
		
	    if (vstgen[i].length() == 1)
			single_gen++;
	}
		
	if (single_gen == number_of_generators)
	{
	    trivial_decomp = true;
		vstgen[0] = "G";
		for (int i=1; i< n_gen_for_1st_vstab; i++)
			vstgen[i] = "";
	}

	if (non_trivial_generators == 0)
	{
		vstgen[0] = "1";
	    trivial_decomp = true; // vertex stabilizer is generated by identity
	}

	single_gen = 0;
	non_trivial_generators=0;

	for (int i=n_gen_for_1st_vstab; i<nvsg ; i++)
	{
		non_trivial_generators += vstgen[i].length();
		
	    if (vstgen[i].length() == 1)
			single_gen++;
	}

	if (single_gen == number_of_generators)
	{
	    trivial_decomp = true;
		vstgen[n_gen_for_1st_vstab] = "G";
		for (int i=n_gen_for_1st_vstab+1; i< nvsg; i++)
			vstgen[i] = "";
	}

	if (non_trivial_generators == 0)
	{
		vstgen[n_gen_for_1st_vstab] = "1";
	    trivial_decomp = true; // vertex stabilizer is generated by identity
	}

	if ( twisted && extension_arc == 0)
		txt << "\nThe twisted track\n\n";
	else if ( twisted )
		txt << "\nThe non-separating and untwisted track\n\n";
	else
		txt << "\nThe separating track\n\n";

	txt << '(';
	for(unsigned int i=0; i< track.size()-1; i++)
		txt << track[i] << ", ";       
	txt << track[track.size()-1] << ')' << endl;

    if (trivial_decomp)
    {
	    txt << "\nGives a trivial decomposition.\n";

		if (echo)
		    cout << "\n\ngives a trivial decomposition.\n";
    }

	txt << "\nEdge stabilizer";
	if (estgen[0] != "1" && estgen[0] != "G")
		txt << " generators.\n";
	else
		txt << ".\n";
		
	for (int i=0; i<nesg; i++)
	{
    	if ( estgen[i] != "" )
			txt << estgen[i] << "\n";
	}

	if (echo)
	{
    	cout << "\n\nEdge stabilizer";

		if (estgen[0] != "1" && estgen[0] != "G")
			cout << " generators.\n";
		else
			cout << ".\n";

    	for (int i=0; i<nesg; i++)
		{
			if ( estgen[i] != "" )
	    		cout << estgen[i] << "\n";
		}
	}
	
	txt << "\n\nFirst vertex stabilizer";
	if (vstgen[0] != "1" && vstgen[0] != "G")
		txt << " generators.\n";
	else
		txt << ".\n";

	for (int i=0; i < n_gen_for_1st_vstab; i++)
	{
    	if ( vstgen[i] != "" )
			txt << vstgen[i] << "\n";
	}

	if (echo)
	{
    	cout << "\n\nFirst vertex stabilizer";
		if (vstgen[0] != "1" && vstgen[0] != "G")
			cout << " generators.\n";
		else
			cout << ".\n";

    	for (int i=0; i < n_gen_for_1st_vstab; i++)
		if ( vstgen[i] != "" )
	    	cout << vstgen[i] << "\n";
	}

	if (!twisted)
	{
    	txt << "\n\nSecond vertex stabilizer";
		if (vstgen[n_gen_for_1st_vstab] != "1" && vstgen[n_gen_for_1st_vstab] != "G")
			txt << " generators.\n";
		else
			txt << ".\n";

    	for (int i=n_gen_for_1st_vstab; i< nvsg; i++)
		{
			if ( vstgen[i] != "" )
	    		txt << vstgen[i] << "\n";
		}

    	if (echo)
		{
			cout << "\n\nSecond vertex stabilizer";
			if (vstgen[n_gen_for_1st_vstab] != "1" && vstgen[n_gen_for_1st_vstab] != "G")
				cout << " generators.\n";
			else
				cout << ".\n";

			for (int i=n_gen_for_1st_vstab; i<nvsg; i++)
		    	if ( vstgen[i] != "" )
				cout << vstgen[i] << "\n";
		}
	}


	/* we can now do some clearing up of memory */

	for (int i=0; i< number_of_1_cells; i++)
		delete tlabels_for_1_cell[i];

	for (int i=0; i< number_of_0_cells; i++)
		delete blabels_for_0_cell[i];

	for (int i=0; i< number_of_1_cells_in_rosette; i++)
	{
		delete b0cbl_for_1_cell[i];
		delete blabels_for_1_cell[i];
	}
	
	for (int i=0; i< number_of_relators; i++)
	{
		delete pair[i];
		delete dperm[i];
		delete darcs[i];
	}
	
	if (!trivial_decomp)
		status |= DECOMPOSE_NON_TRIVIAL;
		
	return status;
}

